from thermochem.janaf import Janafdb
from thermochem.burcat import Elementdb
from collections import defaultdict
import numpy as np
from scipy.constants import physical_constants
import pandas as pd
import matplotlib.pyplot as plt
from monty.serialization import loadfn, dumpfn
import tqdm
from arcs.generate import parse_molecule
from collections import Counter


class RMSEErrors:

    def __init__(self, reaction_dict):
        self.burcat_species = {'H2': 'H2  REF ELEMENT',
                               'N2': 'N2  REF ELEMENT',
                               'O2': 'O2 REF ELEMENT',
                               'S8': 'S8',
                               'H2O': 'H2O',
                               'SO2': 'SO2',
                               'SO3': 'SO3',
                               'H2S': 'H2S',
                               'CS2': 'CS2',
                               'CO': 'CO',
                               'COS': 'COS',
                               'CO2': 'CO2',
                               'HCN': 'HCN',
                               'NO2': 'NO2',
                               'NO': 'NO',
                               'NO3': 'NO3',
                               'N2O': 'N2O',
                               'N2O4': 'N2O4',
                               'O3': 'O3',
                               'H2SO4': 'H2SO4',
                               'HNO3': 'HNO3',
                               'HNO2': 'HNO2',
                               'CH3COOH': 'CH3COOH',
                               'C2H5OH': 'C2H5OH',
                               'C2H6O2': 'C2H6O2',
                               'CH2O': 'CH2O',
                               'CH3CHO': 'CH3CHO',
                               'C6H6': 'C6H6',
                               'C2H6': 'C2H6',
                               'C3H8': 'C3H8',
                               'C2H4': 'C2H4',
                               'C4H8': 'C4H8',
                               'H2O2': 'H2O2(L)',
                               'N2O2': None,
                               'CH4': 'CH4   ANHARMONIC',
                               'NH3': 'NH3 Anharmonic',
                               'H2SO3': None,
                               'H2CO3': None,
                               'HCOOH': 'HCOOH FORMIC ACID',
                               'CH3OH': 'CH3OH(L)',
                               'C6H14O4': None,
                               'CH3NH2': 'CH5N',
                               'C2H5NH2': None,
                               'C3H7NH2': None,
                               'CH3_2_CHCH2NH2': None,
                               'C4H9NH2': None,
                               'CH3(CO)CH3': 'C3H6O Acetone',
                               'CH3_C6H5': 'C7H8  TOLUENE',
                               'C2H5_C6H5': 'C8H10  C6H5C2H5',
                               'C4H10': 'C4H10 n-butane',
                               'C5H12': 'C5H12,i-pentane',
                               'C6H14': 'C6H14,n-hexane',
                               'C7H16': 'C7H16 n-heptane',
                               'C8H18': 'C8H18,n-octane',
                               'C3H6': 'C3H6 propylene'}

        self.ev_conversion = physical_constants["Faraday constant"][0] / 1000

        self.burcatdb = Elementdb()

        self.reaction_dict = reaction_dict

    def get_burcat_gibbs_free_energy_of_reaction(self, temperature):

        reactant_dict = self.reaction_dict["reactants"]
        sum_reactants = 0
        for k, v in reactant_dict.items():
            sum_reactants += (
                self.burcatdb.getelementdata(
                    self.burcat_species[k]
                ).go(T=temperature) / 1000 / self.ev_conversion
            ) * v

        product_dict = self.reaction_dict["products"]
        sum_products = 0
        for k, v in product_dict.items():
            sum_products += (
                self.burcatdb.getelementdata(
                    self.burcat_species[k]
                ).go(T=temperature) / 1000 / self.ev_conversion
            ) * v

        return sum_products - sum_reactants

    def parse_reaction(self, reaction_dict):
        "probably should go in arcs.generate"
        species = Counter()
        for d in list(reaction_dict["reactants"]) + list(reaction_dict["products"]):
            species.update(parse_molecule(d))

        return species

    def rmse_reaction(self, parsed_reaction, rmse_element_wise):
        """probably should go in arcs.generate"""
        summed = []
        for k, v in parsed_reaction.items():
            summed.append(np.square(rmse_element_wise[k] * v))

        return np.sqrt(np.mean(summed))

    def get_element_rmse(self, quantum_reactions, burcat_reactions):
        "root median square error"

        errors = defaultdict(dict)
        for (q, qq), (b, bb), (r, rr) in zip(quantum_reactions.items(),     burcat_reactions.items(), self.reaction_dict.items()):
            species = Counter()
            # list(rr["reactants"]) + list(rr["products"]):
            for d in rr["reactants"]:  # only total species involved?
                species.update(parse_molecule(d))

            try:
                errors[q]["diff"] = qq - bb
                errors[q]["species"] = species
            except Exception:
                pass
        # generate a dataframe
        df = pd.DataFrame(errors).T
        tmp = (
            df['species']
            .apply(pd.Series)
            .stack()
            .reset_index(level=1)
            .rename(columns={'level_1': 'element', 0: 'count'})
        )

        tmp['diff'] = df.loc[tmp.index, 'diff'].values
        tmp["diff_per_atom"] = tmp["diff"] / tmp["count"]

        tmp = tmp.dropna()

        rmse = tmp.groupby('element').apply(
            lambda g: np.sqrt(
                np.median(
                    np.square(
                        g['diff_per_atom'])
                )
            )
        )

        return rmse
