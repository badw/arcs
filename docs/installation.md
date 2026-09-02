# Installation

## Requirements

ARCS requires **Python 3.10 or newer** and builds on a number of scientific
packages (ASE, pymatgen, chempy, NetworkX, SciPy, NumPy, pandas, and others).
These are installed automatically as dependencies.

## From source

The recommended way to install the development version is directly from the
repository:

```bash
git clone https://github.com/badw/arcs.git
cd arcs
pip install .
```

To install in editable mode while developing:

```bash
pip install -e .
```

## Backend data

ARCS ships with coupled-cluster reference data in
`src/arcs/data/quantum_data.json.gz`. Computed using [Psi4](https://psicode.org) 

This file contains the per-molecule energies and vibrational frequencies used to build reaction graphs. 

The list of included compounds and the calculation settings are documented in
`src/arcs/data/README.md`.
