from chempy.equilibria import Equilibrium, EqSystem
from chempy import Substance
import copy
import networkx as nx
import warnings
import numpy as np
import tqdm_pathos
import itertools as it
from arcs.generate import parse_molecule, cost_function
from scipy.optimize import fsolve
import pandas as pd
from collections import defaultdict
from typing import Union, Optional
from tqdm.auto import tqdm


class TableTraversal:
    def __init__(self, table: pd.DataFrame, temperature: float = 298, **kws):
        '''
        TableTraversal contains the algorithm that sorts and filters the reaction table currently:
        1. Filtering:
             * a) filter reactions to show only those with species present in `initial_concentrations` (basic filter)

        2. Sorting:

             * b) filter based on availability (i.e. prioritise larger concentrations first)


             * c) filter based on stoichiometry (i.e. if you have `{O2:50,H2:100}` then 2H2 + O2 = 2H2O might be ranked higher)


             * d) filter based on simplicity (reactions with lots of molecules will be ranked lower - occams razor)


             * e) filter based on gibbs free energy (this could be moved forward)

        '''
        self.table = table
        self.temperature = temperature
        self.meta_cols = (
            'K',
            'G',
            "K_rev",
            "G_rev",
            "availability_score",
            "simplicity_score",
            "gibbs_score",
            "combined_score",
            "quotient_score",
            "probabilities"
        )

        self.__dict__.update(kws)

    def filter_table(
        self,
        initial_concentrations: dict
    ) -> Union[pd.DataFrame, pd.Series]:
        '''
        1. find species with a concentration
        2. only find reactions which fit those species

        i.e. for {H2:10, O2: 10, H2O:10}

        it will be able to find:
        a) 2 H2 + O2 = 2H2O
        b) H2 + O2 = H2O2
        c) O3 = O2
        etc.
        '''

        # get species present in initial concententrations
        present_species = [compound for compound,
                           concentration in initial_concentrations.items() if concentration]

        # filter the reaction table to only show reactions that are possible
        reaction_data = ["K", "G", "K_rev", "G_rev"]
        species_cols = [
            c for c in self.table.columns if c not in reaction_data]

        disallowed = [c for c in species_cols if c not in present_species]

        # disallowed never a          reactant
        forward = (self.table[disallowed] >= 0).all(axis=1)
        # disallowed never a          product
        reverse = (self.table[disallowed] <= 0).all(axis=1)
        filtered_table = self.table[forward | reverse]

        return filtered_table

    def stoichiometry_scores(
            self,
            filtered_table: pd.DataFrame,
            concentrations: dict,
            c_half=1.0
    ) -> pd.Series:

        species_cols = [
            c for c in filtered_table.columns if c not in self.meta_cols]

        def side_score(side):
            if not side:
                return 0.0
            extents = [concentrations[s] / c for s, c in side.items()]
            if max(extents) == 0:
                return 0.0
            ratio_fit = min(extents) / max(extents)
            throughput = min(extents)
            magnitude = throughput / (throughput + c_half)
            return ratio_fit * magnitude

        def score_row(row):
            coeffs = {s: row[s]
                      for s in species_cols if pd.notna(row[s]) and row[s] != 0}
            if not coeffs:
                return 0.0
            reactants = {s: -c for s, c in coeffs.items() if c < 0}
            products = {s:  c for s, c in coeffs.items() if c > 0}

            #  this might be problematic...
            best = max(side_score(reactants), side_score(products))

            # concentration-weighted coverage: trace species count for little
            weights = [concentrations[s] / (concentrations[s] + c_half)
                       for s in coeffs]
            coverage = sum(weights) / len(weights)
            return best * coverage

        return filtered_table.apply(score_row, axis=1)

    def simplicity_scores(
            self,
            filtered_table: pd.DataFrame,
            w_species: float = 1.0,
            w_coeff: float = 0.5,
            direction_scaling: float = 2
    ) -> pd.Series:
        '''
        Rank the reactions in `filtered_table` based on its simplicity. 

        A reaction, i.e. 2 NO + O2 = 2 NO2 (Delta G_reac = -1.36 eV), is ranked based on: 
            a) sum of the coefficients (i.e. 3 reactants, 2 products)
                * coefficients are weighted by `w_coeff` = 0.5 * 
            b) number of species ( i.e. 2 reactants, 1 product)
                * species are weighted by `w_species` = 1.0 * 

        The score determines the direction of the reaction based on the Gibbs Free energy of reaction.
        Therefore the simplicity score of the reactants are weighted above those of the products by `direction_scaling` = 2 

        The score is 1 / exp(reactants_score + products_score)

        RETURNS: 
        pd.Series

        '''

        species_cols = [
            c for c in filtered_table.columns if c not in self.meta_cols]

        def score_row(row):

            coeffs = {s: row[s] for s in species_cols
                      if pd.notna(row[s]) and row[s] != 0}

            reactants = {s: -c for s, c in coeffs.items() if c < 0}

            products = {s:  c for s, c in coeffs.items() if c > 0}

            reactants_score = len(reactants) * w_species + \
                sum(reactants.values()) * w_coeff
            products_score = (len(products) * w_species +
                              sum(products.values()) * w_coeff) / direction_scaling

            #  exponential function
            return 1 / (1.0 + np.exp(reactants_score+products_score))

        simplicity_score = filtered_table.apply(score_row, axis=1)

        return simplicity_score / simplicity_score.max()  # normalise it

    def reaction_quotient_scores(
            self,
            filtered_table: pd.DataFrame,
            concentrations: dict,
            scale: float = 2.0,
            floor: float = 1e-30
    ) -> pd.Series:
        """Score reactions by distance-from-equilibrium using Q vs K.    

        Q = product([conc]^coeff) over all species (products positive exponent,
        reactants negative). Driving force = ln(K / Q); forward-favorable > 0.
        Score = sigmoid(ln(K/Q) / scale) → 1 strongly forward, 0 strongly reverse.
        """
        species_cols = [
            c for c in filtered_table.columns if c not in self.meta_cols]

        def score_row(row):
            K = row["K"]
            if pd.isna(K) or K <= 0:
                return 0.0

            ln_Q = 0.0
            for s in species_cols:
                coeff = row[s]
                if pd.isna(coeff) or coeff == 0:
                    continue
                c = concentrations.get(s, 0.0)
                # avoid log(0); trace ≠ absent
                c = max(c, floor)
                ln_Q += coeff * np.log(c)         # Σ coeff · ln(conc)

            # ln(K/Q); >0 ⇒ forward favorable
            driving_force = np.log(K) - ln_Q
            return 1.0 / (1.0 + np.exp(-driving_force / scale))
        quotient_scores = filtered_table.apply(score_row, axis=1)

        return quotient_scores / quotient_scores.max()

    def gibbs_scores(
            self,
            filtered_table: pd.DataFrame,
            concentrations: dict,
            per_reactant_atom: bool = True,
    ) -> pd.Series:
        '''
        Rank the reactions in `filtered_table` based on their Gibbs Free Energy of Reaction. 

        The score is direction aware and concentration aware. i.e. if you have: 

        2 NO + O2 = 2 NO2 (Delta G_reac = -1.36eV) 

        and the concentrations: {"NO2":20} 

        It will assume that the products in this instance are reactants and use Delta G_reac = 1.36 eV 

        Gibbs Free energies of reaction are normalised per reactant atom, based on arcs.generate.cost_function which takes into account self.temperature.

        Returns 
        pd.Series
        '''

        species_cols = [
            c for c in filtered_table.columns if c not in self.meta_cols]

        def all_present(side):
            return bool(side) and all(concentrations[s] > 0 for s in side)

        def score_row(row):
            coeffs = {
                s: row[s] for s in species_cols if pd.notna(row[s]) and row[s] != 0
            }
            reactants = [s for s, c in coeffs.items() if c < 0]
            products = [s for s, c in coeffs.items() if c > 0]

            candidates = []
            if all_present(reactants):
                candidates.append(
                    cost_function(
                        gibbs_free_energy=row["G"],
                        temperature=self.temperature,
                        reactants={s: int(-c)
                                   for s, c in coeffs.items() if c < 0},
                        normalise_by_reactant_atoms=per_reactant_atom
                    )
                )
            if all_present(products):
                candidates.append(
                    cost_function(
                        gibbs_free_energy=row["G_rev"],
                        temperature=self.temperature,
                        reactants={s: int(c)
                                   for s, c in coeffs.items() if c > 0},
                        normalise_by_reactant_atoms=per_reactant_atom
                    )
                )
            if not candidates:
                # can't feed either direction, however need additional any clause if reaction can indeed proceed...
                return 0

            #  1 / in order to make it the same as the other scores...
            return 1 / min(candidates)

        gibbs_score = filtered_table.apply(score_row, axis=1)

        return gibbs_score / gibbs_score.max()

    def combined_scores(self,
                        filtered_table: pd.DataFrame,
                        score_weighting={'availability_score': 1.0,
                                         'simplicity_score': 1.0, 'gibbs_score': 1.0, "quotient_score": 1.0},
                        method='sum',
                        normalise=True,
                        ):

        # needs all the scores present in the dataframe
        """Combine sub-score columns into one final score.
        weights = {'availability_score': 1.0, 'simplicity_score': 1.0, 'gibbs_score': 1.0} (ratio of importance)
        method='geometric': weighted geometric mean — a near-zero on ANY factor
            tanks the total (use when all factors are necessary).
        method='sum': weighted arithmetic mean — factors compensate for each other.
        normalize: rescale each factor to [0,1] across the set first, so no factor
            dominates just because of its raw range.
        """
        score_cols = ('availability_score', 'simplicity_score',
                      'gibbs_score', "quotient_score")
        w = np.array([score_weighting[c] for c in score_cols], dtype=float)
        w = w / w.sum()                                  # weights sum to 1

        S = filtered_table[list(score_cols)].astype(float).copy()

        # per-column min-max to [0,1]
        if normalise:  #  can remove this eventially
            rng = S.max() - S.min()
            # avoid div-by-zero on flat cols
            rng = rng.replace(0, 1)
            S = (S - S.min()) / rng

        if method == 'geometric':  #  check this
            eps = 1e-9                                    # keep log finite at 0
            final = np.exp((np.log(S + eps) * w).sum(axis=1))
        elif method == 'sum':
            final = (S * w).sum(axis=1)
        else:
            raise ValueError("method must be 'geometric' or 'sum'")

        return final

    def filter_and_sort_reactions_table(
            self,
            concentrations: dict,
            #  experimental - may not be needed
            method: str = "sum",
            score_weighting: dict = {
                'availability_score': 2.0,
                'simplicity_score': 1.0,
                'gibbs_score': 2.0,
                "quotient_score": 1.0
            },
            w_species: float = 1.0,
            w_coeff: float = 0.5,
            normalise: bool = False
            # kws
    ) -> pd.DataFrame:
        # filter table
        filtered_table = self.filter_table(
            initial_concentrations=concentrations
        )

        # filter based on availability
        filtered_table["availability_score"] = availability_score = self.stoichiometry_scores(
            filtered_table=filtered_table,
            concentrations=concentrations
        )

        #  filter based on simplicity:
        filtered_table["simplicity_score"] = self.simplicity_scores(
            filtered_table=filtered_table, w_species=w_species, w_coeff=w_coeff)

        #  filter based on reaction quotient
        filtered_table["quotient_score"] = self.reaction_quotient_scores(
            filtered_table=filtered_table, concentrations=concentrations
        )

        # filter based on Gibbs Free Energy:
        gibbs_score = self.gibbs_scores(
            filtered_table=filtered_table, concentrations=concentrations
        )
        filtered_table["gibbs_score"] = gibbs_score

        combined_score = self.combined_scores(
            filtered_table=filtered_table,
            method=method,
            score_weighting=score_weighting,
            normalise=normalise
        )

        filtered_table["combined_score"] = combined_score

        final_table = filtered_table.sort_values(
            by="combined_score", ascending=False
        )

        return final_table

    def row_to_reaction(self, row):  #  make this take a whole table instead
        """Turn a DataFrame row into a reaction dict.

        Species columns hold stoichiometric coefficients:
        negative = reactant, positive = product, NaN/0 = not involved.
        """

        reactants, products = {}, {}
        for species, coeff in row.items():
            if species in self.meta_cols or pd.isna(coeff) or coeff == 0:
                continue
            (reactants if coeff < 0 else products)[species] = abs(coeff)

        def fmt(c):                       # drop trailing .0 on whole numbers
            return int(c) if float(c).is_integer() else c

        def side(d):
            return " + ".join(f"{fmt(c)} {s}" for s, c in d.items())

        result = {"reaction_string": f"{side(reactants)} = {side(products)}"}
        for m in self.meta_cols:
            if m in row.index:
                result[m] = row[m]
        result["reactants"] = {s: fmt(c) for s, c in reactants.items()}
        result["products"] = {s: fmt(c) for s, c in products.items()}
        return result

    def choose_reaction(
            self,
            filtered_table: pd.DataFrame,
            score_col='combined_score',
            choice_temp=0.8,
            mode='softmax',
            return_probabilities: bool = False
    ):
        """Probabilistic pick that sharpens a flat score landscape.    

        mode='softmax': p ∝ exp(score / T). Low T → nearly greedy, high T → uniform.
        mode='rank':    sample by rank, discarding the (flat) magnitudes entirely.
        """

        if mode == 'rank':
            s = filtered_table[score_col].rank(
                ascending=True).to_numpy(dtype=float)
        else:
            s = filtered_table[score_col].to_numpy(dtype=float)
            # standardise so T is meaningful
            s = (s - s.mean()) / (s.std() + 1e-12)

        z = (s - s.max()) / choice_temp
        p = np.exp(z)
        p /= p.sum()

        idx = np.random.default_rng().choice(filtered_table.index, p=p)
        if return_probabilities:
            return p
        else:
            return filtered_table.loc[idx]

    def generate_chempy_eqsystem(
            self,
            reaction_dict: dict,
    ) -> EqSystem:
        """
        given a reaction index form a chempy.equilibria.EqSystem

        eventually this will be deprecated as it is a speed bottleneck

        needs to involve charged species

        """
        reactants = reaction_dict["reactants"]
        products = reaction_dict["products"]
        k = reaction_dict["K"]

        substances = {}
        for compound in list(it.chain(*[list(reactants) + list(products)])):
            substances[compound] = Substance.from_formula(
                compound, **{"charge": 0})

        equation = Equilibrium(reac=reactants, prod=products, param=k)
        try:
            return EqSystem([equation], substances=substances)
        except Exception:
            return None

    def chempy_equilibrium_concentrations(
        self,
        concentrations: dict,
        equilibrium_reaction: EqSystem,
        chempy_sane=True,
    ) -> dict:
        """
        generate equilibrium concentrations
        if the reaction is "sane" and a "success" then it returns the equilibrium concentrations as a dict     of concentrations
        elsewise
        return None
        """

        warnings.simplefilter("ignore")
        _concs = copy.deepcopy(concentrations)
        try:
            result = equilibrium_reaction.solve(init_concs=_concs)
            # assert result.success and result.sane
            if chempy_sane:
                assert result.success and result.sane
            else:
                assert result.success
            for compound, concentration in enumerate(result.conc):
                _concs[equilibrium_reaction.substance_names()[compound]
                       ] = concentration
            return _concs
        except Exception:
            return None

    def inner_loop(
            self,
            final_table: pd.DataFrame,
            concentrations: dict,
            max_attempts: int = 10,
            choice_temp: float = 0.3,
            mode="softmax"
    ):
        '''
        1. weighted random picking of a reaction
        2. calculate equilibrium concentrations
        3. save {reaction_stats:reaction,concentrations:concentrations}

        to avoid an infinite loop, breaks after max_attempts.
        '''
        new_concs = None
        chosen_reaction = None
        attempts = 1

        while not new_concs and not attempts == max_attempts:
            # 1.
            chosen_reaction = self.choose_reaction(
                final_table,
                choice_temp=choice_temp,
                mode=mode
            ).name
            # 2.
            reaction_dict = self.row_to_reaction(
                final_table.T[chosen_reaction])
            # 3.
            eqsystem = self.generate_chempy_eqsystem(reaction_dict)
            # 4.
            new_concs = self.chempy_equilibrium_concentrations(
                copy.deepcopy(concentrations),
                eqsystem,
                chempy_sane=True
            )
            chosen_reaction = chosen_reaction if new_concs else None
            attempts += 1

        return new_concs, chosen_reaction

    @staticmethod
    def check_convergence(concentration_stats, tol=0.5, frac=0.5, how_far_back=5):
        """
        Returns True if the majority of rows have all values within `tol` of each other.

        tol  : max allowed spread (max - min) within a row
        frac : fraction of rows that must converge (0.5 = strict majority)
        """

        how_far_back = 5

        dfs = [
            pd.DataFrame(
                [
                    x for x in concentration_stats[i][1:] if x
                ]
            ) for i in range(
                len(concentration_stats)-how_far_back, len(concentration_stats)
            )
        ]
        dfs = pd.DataFrame(
            [df[df > 0].dropna(axis=1, how="all").mean() for df in dfs]
        ).T.round(2)

        row_spread = dfs.max(axis=1) - dfs.min(axis=1)
        converged = row_spread <= tol
        return converged.mean() > frac

    def outer_loop(self,
                   initial_concentrations: Union[dict, pd.Series],
                   path_length: int = 50,
                   inner_loop_runs: int = 50,
                   method: str = "sum",
                   score_weighting: dict = {
                       'availability_score': 2.0, 'simplicity_score': 1.0, 'gibbs_score': 2.0, "quotient_score": 1.0
                   },
                   normalise=False,
                   w_species: float = 1.0,
                   w_coeff: float = 0.5,
                   choice_temp: float = 0.8,
                   mode="softmax",
                   convergence_how_far_back=5
                   ):
        reaction_stats = defaultdict(list)
        concentration_stats = defaultdict(list)
        scoring_stats = defaultdict(list)
        updated_concs = copy.deepcopy(initial_concentrations)

        for i in tqdm(range(path_length)):
            #  first step is the initial concentrations
            reaction_stats[i] = [None]
            concentration_stats[i] = [updated_concs]
            scoring_stats[i] = [None]
            if i > 0:
                #  takes the mean from the previous step
                updated_concs = copy.deepcopy(
                    pd.DataFrame([x for x in concentration_stats[i-1]
                                 [1:] if not x == None]).mean().to_dict()
                )
                final_table = self.filter_and_sort_reactions_table(
                    # self.table,
                    concentrations=updated_concs,
                    method=method,
                    score_weighting=score_weighting,
                    normalise=normalise,
                    w_species=w_species,
                    w_coeff=w_coeff
                )
                scoring_stats[i] = final_table.filter(
                    ["availability_score", "gibbs_score", "simplicity_score", "quotient_score", "combined_score"])

            #  the first entry is always the previous concentrations (mean)
            if i == 0:
                reaction_stats[i].append(None)
                concentration_stats[i].append(updated_concs)
            else:
                for n in range(1, inner_loop_runs+1):
                    _conc, _reaction = self.inner_loop(
                        final_table,
                        updated_concs,
                        choice_temp=choice_temp,
                        mode=mode
                    )
                    reaction_stats[i].append(_reaction)
                    concentration_stats[i].append(_conc)
                if i > convergence_how_far_back:
                    converged = self.check_convergence(
                        concentration_stats=concentration_stats, how_far_back=convergence_how_far_back)
                    if converged:
                        print("convergence reached.")
                        break

        return concentration_stats, reaction_stats, scoring_stats

####################################################################################################


class GraphTraversal:
    def __init__(self, graph, max_reaction_length=5, **kws):
        '''
        Algorithm that traverses the graph - legacy version of ARCS.
        '''
        self.graph = graph

        # default values:
        self.exclude_co2 = True
        self.max_compounds = 5
        self.discovery_threshold = 5  # % percent
        self.maximum_reaction_number = 10
        self.max_steps = 5
        self.ncpus = 4
        self.ceiling = 2000
        self.scale_largest = 10
        self.rank_small_reactions_higher = True
        self.rank_by_number_of_atoms = True
        self.shortest_path_method = "Djikstra"
        self.__dict__.update(kws)

    def length_multiplier(self, candidate_reaction: int, **kws) -> float:
        """
        given a candidate reaction,

        if self.rank_small_reactions_higher == True, then return the sum of the coefficients as a multiplier

        i.e. H2 + 1/2 O2 = H2O has a length multiplier of 2.5

        if self.rank_by_number_of_atoms = True, then rank by the number of atoms in the reactants (which = num_atoms_in_products)

        i.e. H2 + 1/2 O2 = H2O has a length multiplier of 3

        """
        self.__dict__.update(kws)

        if self.rank_small_reactions_higher:
            if self.rank_by_number_of_atoms:
                reaction_dict = self.graph.nodes[candidate_reaction]["reaction"]
                reactants = reaction_dict["reactants"]
                num_atoms = []
                for species, coefficient in reactants.items():
                    num_atoms.append(
                        np.sum(list(parse_molecule(species).values())) *
                        coefficient
                    )
                return np.sum(num_atoms)
            else:
                reaction_dict = self.graph.nodes[candidate_reaction]["reaction"]
                num_reactants = np.sum(
                    list(reaction_dict["reactants"].values()))
                num_products = np.sum(list(reaction_dict["products"].values()))
                return num_reactants + num_products
        else:
            return 1

    def check_reactant_atoms(
        self, reaction_index: int, weighted_random_compounds: list, **kws
    ) -> bool:
        """
        takes a reaction index, chosen random_compounds, and checks for atom balance
        returns bool
        """
        self.__dict__.update(kws)

        reac_atoms = list(
            dict.fromkeys(
                list(
                    it.chain(*[
                        list(parse_molecule(x))
                        for x in self.graph.nodes()[reaction_index]["reaction"][
                            "reactants"
                        ]
                    ])
                )
            )
        )
        random_compounds_atoms = list(
            dict.fromkeys(
                list(
                    it.chain(*[
                        list(parse_molecule(x)) for x in weighted_random_compounds
                    ])
                )
            )
        )
        if sorted(reac_atoms) == sorted(random_compounds_atoms):
            return True
        else:
            return False

    def scale_large_concentrations(self, concentrations: dict, **kws) -> dict:
        """
        function that takes a dict of concentrations and scales abnormally large concentrations (above ceiling DEFAULT = 3000%) and scales them by scale_highest (DEFAULT = 10% of original value)

        this is to be used with self.get_weighted_random_compounds
        """
        self.__dict__.update(kws)

        median_conc = np.median([v for v in concentrations.values() if v > 0])
        species_above_ceiling = {
            k: v
            for k, v in concentrations.items()
            if v > (median_conc * (1 + (self.ceiling / 100)))
        }
        # modify the ceiling by scaling it down to a suitable value
        # should still max out if concentrations become way to high
        for k, v in species_above_ceiling.items():
            concentrations[k] = v * 1 / self.scale_largest

        return concentrations

    def get_weighted_random_compounds(self, concentrations: dict, **kws) -> list:
        """
        given a dictionary of concentrations e.g. {'H2O':100,'NO2':50} a weighted ranking can be returned with probabilities given a discovery threshold DEFAULT = 5%.

        exceedingly large concentrations (up to ceiling % DEFAULT = 1000% above the median concentration) that may occur are scaled using self.scale_large_concentrations (scaled with scale_largest DEFAULT = 10% of original value) such that reactions may continue even with very large concentrations of species up to a point.

        returns a list with length up to max_compounds depending on the discovery_threshold and scale_largest factors.

        CO2 is by default excluded (exclude_co2 = True) as it is considered background, however this can be turned on if you want to test CO2 containing reactions.
        """
        self.__dict__.update(kws)

        concs = copy.deepcopy(concentrations)

        if self.exclude_co2 and "CO2" in concs:
            # CO2 will always be too large as it is the background
            del concs["CO2"]

        # scale potential large concentrations
        concs = self.scale_large_concentrations(
            concentrations=concs, scale_largest=self.scale_largest, ceiling=self.ceiling
        )

        # get the probabilities based upon relative concentrations:
        p_1 = {k: v / sum(concs.values()) for k, v in concs.items()}
        # now filter based upon the probability threshold: (discovery)
        p_2 = {k: v for k, v in p_1.items() if v >=
               self.discovery_threshold / 100}
        if not p_2:
            return []
        # remake the probabilities
        p_3 = {k: v / sum(p_2.values()) for k, v in p_2.items()}
        # make a list of choices based upon the probabilities
        # orig:
        available = list(
            np.random.choice(a=list(p_3), size=len(
                concs) * 10, p=list(p_3.values()))
        )
        # experimental:
        # available = list(
        #    np.random.choice(a=list(p_3), size=len(
        #        p_3), replace=False, p=list(p_3.values()))
        # )
        #  now make a list max_compounds long of random choices based on available
        choices = {}
        for i in range(self.max_compounds):
            try:
                compound = np.random.choice(available)
                choices[compound] = p_3[compound]

                available = list(
                    filter(lambda a: a != list(choices)[i - 1], available))
            except ValueError:
                pass
        return list(choices)[0: np.random.randint(2, self.max_compounds)]

    def get_weighted_reaction_rankings(self, weighted_random_compounds: list, **kws) -> dict:
        """
        returns a dictionary of {<reaction_index>:<weighting>} given a list of weighted_random_compounds from self.get_weighted_random_compounds (needs at least 2 to give a result)

        algorithm follows:
        given weighted_random_compounds = ['NO2','H2O','O2']
        1. generates combinations of useable compounds
            i.e. [['NO2','H2O'],['NO2','O2'],['H2O','O2']]
        2. for each combination, generate a list of shortest paths using networkx.shortest_paths.all_shortest_paths
        3. check that the reactant

        """
        self.__dict__.update(kws)

        # return None if there isn't enough to make a reaction
        if len(weighted_random_compounds) <= 1:
            return None
        # 1. generate possible combinations
        combinations = list(it.combinations(weighted_random_compounds, 2))
        # 2. generate shortest path possibilities from the combinations
        possibilities = []
        for compounds in combinations:
            possibilities.extend([
                x[1]
                for x in list(
                    nx.shortest_paths.all_shortest_paths(
                        G=self.graph,
                        source=compounds[0],
                        target=compounds[1],
                        method=self.shortest_path_method,
                    )
                )
            ])

        # 3. check that all reactant atoms (=product atoms) are accounted for.
        # this is so that self.get_chempy_equilibrium_concentrations gives a result.
        possibilities = [
            i
            for i in possibilities
            if self.check_reactant_atoms(
                reaction_index=i,
                weighted_random_compounds=weighted_random_compounds,
                **kws,
            )
        ]
        # 4. rank the possibilities based on edge weight and a length_multiplier
        # the length_multiplier is based on number of reaction coefficients
        # idea: perhaps by number of atoms ?
        rankings = {}
        for i, reaction in enumerate(possibilities):
            for compound in weighted_random_compounds:
                try:
                    weight = self.graph.get_edge_data(u=compound, v=reaction)[0][
                        "weight"
                    ] * self.length_multiplier(reaction, **kws)
                except TypeError:
                    pass
            rankings[reaction] = (
                weight  # this should be ammended to return None as well
            )

        rankings = dict(
            sorted(rankings.items(), key=lambda item: item[1])[
                0: self.maximum_reaction_number
            ]
        )

        return rankings

    @staticmethod
    def choose_reaction(ranked_reactions: dict) -> int:
        """
        given a dictionary of ranked reactions from self.get_weighted_reaction_rankings
        chose a reaction based on weights and probabilities
        """
        weights = {
            k: 1 / v**2 for k, v in ranked_reactions.items()
        }  # here higher is better
        # added a square multiplier to force more the larger coefficients
        probabilities = {k: v / sum(weights.values())
                         for k, v in weights.items()}
        chosen_reaction = np.random.choice(
            [
                np.random.choice(
                    a=list(probabilities.keys()),
                    size=len(probabilities) * 10,
                    p=list(probabilities.values()),
                )
            ][0]
        )
        return chosen_reaction

    def generate_chempy_eqsystem(self, index: int) -> EqSystem:
        """
        given a reaction index form a chempy.equilibria.EqSystem

        eventually this will be deprecated as it is a speed bottleneck

        needs to involve charged species

        """
        node_dict = self.graph.nodes[index]
        reactants = node_dict["reaction"]["reactants"]
        products = node_dict["reaction"]["products"]
        k = node_dict["equilibrium_constant"]

        substances = {}
        for compound in list(it.chain(*[list(reactants) + list(products)])):
            substances[compound] = Substance.from_formula(
                compound, **{"charge": 0})

        equation = Equilibrium(reac=reactants, prod=products, param=k)
        try:
            return EqSystem([equation], substances=substances)
        except Exception:
            return None

    @staticmethod
    def chempy_equilibrium_concentrations(
        concentrations: dict,
        equilibrium_reaction: EqSystem,
        chempy_sane=True,
    ) -> dict:
        """
        generate equilibrium concentrations
        if the reaction is "sane" and a "success" then it returns the equilibrium concentrations as a dict of concentrations
        elsewise
        return None
        """

        warnings.simplefilter("ignore")
        _concs = copy.deepcopy(concentrations)
        try:
            result = equilibrium_reaction.solve(init_concs=_concs)
            # assert result.success and result.sane
            if chempy_sane:
                assert result.success and result.sane
            else:
                assert result.success
            for compound, concentration in enumerate(result.conc):
                _concs[equilibrium_reaction.substance_names()[compound]
                       ] = concentration
            return _concs
        except Exception:
            return None

    def random_walk(
        self,
        initial_concentrations: dict,
        chempy_sane=True,  # typically for very large
        ** kws,
    ) -> dict:
        """
        does a random sampling of the reaction network with max_steps  DEFAULT = 10.

        for each sample step:
        1. get weighted random compounds
        2. get ranked reactions
        3. choose a reaction
        3. generate a chempy eqsystem
        4. calculate the equilibrium concentrations
        5. update the concentrations and reaction statistics
        """
        self.__dict__.update(kws)

        concentrations = {0: initial_concentrations}
        reactionstats = {0: None}
        i = 0
        for step in range(1, self.max_steps + 1):
            _concentrations = copy.deepcopy(concentrations[i])
            # 1 grab weighted_random_compounds
            weighted_random_compounds = self.get_weighted_random_compounds(
                concentrations=_concentrations,
                exclude_co2=self.exclude_co2,
                max_compounds=self.max_compounds,
                discovery_threshold=self.discovery_threshold,
                scale_largest=self.scale_largest,
                ceiling=self.ceiling,
                **kws,
            )
            # 2 grab reaction rankings
            ranked_reactions = self.get_weighted_reaction_rankings(
                weighted_random_compounds=weighted_random_compounds, **kws
            )
            # 3 if no reactions found then break the for loop
            if not ranked_reactions:
                break

            # 4 choose a reaction
            chosen_reaction_index = self.choose_reaction(
                ranked_reactions=ranked_reactions
            )
            # 5 generate a chempy eqsystem
            eqsystem = self.generate_chempy_eqsystem(
                index=chosen_reaction_index)

            # 6 get equilibrium_concentrations and update relevant dictionaries.
            final_concentrations = self.chempy_equilibrium_concentrations(
                concentrations=_concentrations,
                equilibrium_reaction=eqsystem,
                chempy_sane=chempy_sane,
            )

            if final_concentrations:
                i += 1
                concentrations[i] = final_concentrations
                try:
                    reactionstats[i] = {
                        "reaction": self.graph.nodes[chosen_reaction_index]["reaction"],
                        "equilibrium_constant": self.graph.nodes[chosen_reaction_index][
                            "equilibrium_constant"
                        ],
                    }
                except KeyError:
                    reactionstats[i] = {
                        "reaction": self.graph.nodes[chosen_reaction_index]["reaction"],
                        "log_equilibrium_constant": self.graph.nodes[chosen_reaction_index][
                            "log_equilibrium_constant"
                        ],
                    }

        return {"concentrations": concentrations, "reaction_statistics": reactionstats}

    def sampling_function(self, iterable, **kws):
        """
        sampling function to be multiprocessed - runs one random walk
        """
        self.__dict__.update(kws)
        initial_concentrations = self.initial_concentrations
        return self.random_walk(initial_concentrations=initial_concentrations, **kws)

    def sample(
        self,
        initial_concentrations: dict,
        nsamples: int,
        ncpus: int,
        tqdm_kws: dict = {},
        **kws,
    ) -> dict:
        """
        samples the graph network nsamples DEFAULT = 1000
        multiprocessed with ncpus DEFAULT = 4
        """
        self.__dict__.update(kws)
        self.initial_concentrations = initial_concentrations

        data = tqdm_pathos.map(
            self.sampling_function,
            list(range(nsamples)),
            n_cpus=ncpus,
            tqdm_kwargs=tqdm_kws,
            **kws,
        )

        return data
