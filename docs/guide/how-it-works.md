# How ARCS works

<!--
ARCS models the impurity chemistry of stored CO<sub>2</sub> as a search over a
**reaction network**. Rather than integrating a full kinetic model, it treats the
problem probabilistically: reactions that are thermodynamically favourable and
involve species that are actually present are more likely to be visited.

## The three stages

### 1. Thermodynamics → graph

Each molecule has first-principles data: a total energy, vibrational frequencies,
point group, spin, and geometry. From these, ASE's `IdealGasThermo` gives the
Gibbs free energy of each species at a given temperature and pressure, and hence
the **Gibbs free energy of reaction** $\Delta G_\mathrm{r}$ and the
**equilibrium constant** $K$ for every candidate reaction.

These reactions become a directed graph. There are two kinds of node:

- **Compound nodes** (strings such as `"H2O"`)
- **Reaction nodes** (integer indices)

Edges connect compounds to the reactions they participate in and back again, in
both the forward and reverse directions. Every edge carries a **weight** derived
from a thermodynamic cost function, so "cheap" edges correspond to favourable
chemistry.

### 2. Random walks over the graph

Starting from an initial gas mixture, ARCS performs many independent
**random walks**. At each step of a walk it:

1. Picks a handful of compounds, weighted by their current concentration.
2. Finds the shortest paths between pairs of those compounds through the graph —
   these paths *are* candidate reactions.
3. Ranks the candidates by edge weight and reaction size, and chooses one
   probabilistically.
4. Solves that reaction to equilibrium with `chempy`, updating the concentrations.

Because the choice of compounds and reactions is random but concentration- and
thermodynamics-weighted, common pathways are visited often and rare ones
occasionally — exactly the statistics we want to collect.

### 3. Aggregate statistics

A single walk is noisy. Running thousands of them and pooling the results yields
robust statistics: how frequently each reaction occurs, and how each species'
concentration changes on average. That ensemble is the actual scientific output.

## Why a graph and Monte Carlo?

The full combinatorial space of reactions among dozens of species is enormous.
Representing it as a weighted graph lets ARCS use fast shortest-path algorithms to
propose only *chemically connected, atom-balanced* reactions, while the
Monte&nbsp;Carlo sampling explores the space in proportion to how likely each
pathway is — without ever enumerating all of it.

## Reading the rest of the guide

| Stage | Guide page | Key class |
| --- | --- | --- |
| Build the graph | [Building the reaction graph](building-the-graph.md) | [`GraphGenerator`][arcs.generate.GraphGenerator] |
| Seed concentrations | [Initial concentrations](concentrations.md) | [`GenerateInitialConcentrations`][arcs.generate.GenerateInitialConcentrations] |
| Sample | [Sampling the network](sampling.md) | [`Traversal`][arcs.traversal.Traversal] |
| Analyse | [Analysing results](analysis.md) | [`AnalyseSampling`][arcs.analysis.AnalyseSampling] |
-->