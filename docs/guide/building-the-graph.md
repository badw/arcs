# Building the reaction graph
<!--
[`GraphGenerator`][arcs.generate.GraphGenerator] is the entry point for turning
quantum data into a weighted reaction graph. It has two public builders:

- [`from_file`][arcs.generate.GraphGenerator.from_file] — load from a `.json` /
  `.json.gz` file on disk.
- [`from_dict`][arcs.generate.GraphGenerator.from_dict] — build from an
  already-loaded dictionary.

```python
from arcs.generate import GraphGenerator

graph = GraphGenerator().from_file(
    filename="src/arcs/data/quantum_data.json.gz",
    temperature=248,
    pressure=20,
    max_reaction_length=4,
)
```

## Parameters that shape the graph

| Parameter | Meaning |
| --- | --- |
| `temperature` | Temperature in **K** at which Gibbs energies are evaluated. |
| `pressure` | Pressure in **bar** (converted to Pa internally). |
| `max_reaction_length` | Maximum number of reactants + products in a reaction. Larger values admit more complex reactions. |
| `log_K` | If `True`, store the logarithm of the equilibrium constant instead of `K` itself — useful when `K` spans many orders of magnitude. |
| `filter_large_gibbs` | If `True`, drop reactions whose \|ΔG\| exceeds the internal `gibbs_filter` (709 eV) to avoid numerical overflow in `K`. |

## The input data format

Both builders expect a dictionary (or a file that deserialises to one) with a
per-compound entry plus a `reactions` entry:

```python
{
    "H2O": {
        "atoms": ...,          # ase.Atoms.todict()
        "pointgroup": "C2v",
        "spin": 0,
        "rotation_num": 2,
        "islinear": "nonlinear",
        "energy": ...,          # potential energy (eV)
        "vibrations": ...,      # vibrational energies (eV)
    },
    # ... more compounds ...
    "reactions": {
        0: {
            "reaction_string": "1 H2O + 1 H2CO = 1 H2 + 1 CH2O2",
            "reactants": {"H2O": 1, "H2CO": 1},
            "products": {"H2": 1, "CH2O2": 1},
        },
        # ... more reactions ...
    },
}
```

Reaction lists can be generated with the companion tool
[reactit](https://github.com/badw/reactit). The per-compound quantum data is
produced by the helper classes described below.

## From energies to weights

### Gibbs free energy and equilibrium constant

[`ReactionGibbsandEquilibrium`][arcs.generate.ReactionGibbsandEquilibrium] does the
thermodynamics. For each species it calls ASE's `IdealGasThermo` to obtain a Gibbs
free energy at the requested temperature and pressure, then combines them into a
reaction free energy:

$$
\Delta G_\mathrm{r} = \sum_\text{products} \nu_i\,G_i - \sum_\text{reactants} \nu_j\,G_j
$$

The equilibrium constant follows from $\Delta G_\mathrm{r} = -k_\mathrm{B} T \ln K$:

$$
K = \exp\!\left(\frac{-\Delta G_\mathrm{r}}{(k_\mathrm{B}/e)\,T}\right)
$$

With `log_K=True`, the exponential is skipped and the exponent itself is stored.

### The cost function

Edge weights come from [`cost_function`][arcs.generate.GraphGenerator.cost_function],
adapted from [Nature Communications 12, 2021](https://www.nature.com/articles/s41467-021-23339-x):

$$
\text{cost} = \ln\!\left(1 + \frac{273}{T}\,\exp\!\left(\frac{\Delta G_\mathrm{r}}{N_\text{atoms}}\right)\right)
$$

where $N_\text{atoms}$ is the number of reactant atoms (normalisation can be turned
off). Forward and reverse edges get the cost of their respective directions, so a
downhill reaction is cheap forwards and expensive backwards.

## Generating quantum data

If you are starting from raw electronic-structure calculations rather than the
bundled data, three helper classes turn calculator output into the per-compound
dictionary above:

- [`GetEnergyandVibrationsVASP`][arcs.generate.GetEnergyandVibrationsVASP] — reads a
  VASP relaxation (`POSCAR`/`OUTCAR`) and a DFPT vibrations directory.
- [`GetEnergyandVibrationsPsi4`][arcs.generate.GetEnergyandVibrationsPsi4] — takes an
  `ase.Atoms`, an energy (Hartree), and frequencies (cm⁻¹) from Psi4.
- [`GetEnergyandVibrationsAseCalc`][arcs.generate.GetEnergyandVibrationsAseCalc] —
  wraps a generic ASE atoms + vibrations calculation.

Each exposes an `as_dict()` method returning exactly the structure ARCS expects,
including point group, spin, and rotational symmetry number derived automatically.

See the full [`arcs.generate` API reference](../api/generate.md) for every method
and argument.
-->