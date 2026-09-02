# `arcs.generate`

The core module: it converts first-principles data into a weighted reaction graph
and provides the initial-concentration helpers.

## Graph construction

::: arcs.generate.GraphGenerator

## Thermodynamics

::: arcs.generate.ReactionGibbsandEquilibrium

## Initial concentrations

::: arcs.generate.GenerateInitialConcentrations

## Quantum-data helpers

These classes convert electronic-structure output into the per-compound
dictionaries ARCS consumes.

::: arcs.generate.GetEnergyandVibrationsVASP

::: arcs.generate.GetEnergyandVibrationsPsi4

::: arcs.generate.GetEnergyandVibrationsAseCalc

## Module functions

::: arcs.generate.parse_molecule

::: arcs.generate.get_compound_directory
