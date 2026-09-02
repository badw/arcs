# ARCS



**Automated Reactions for CO<sub>2</sub> Storage**

ARCS is a Python toolkit for exploring the network of chemical reactions that can
occur among impurities in stored CO<sub>2</sub>. It turns first-principles
thermochemistry into a weighted reaction graph, samples that graph with a
Monte&nbsp;Carlo random walk, and reports which reactions dominate and how species
concentrations evolve.

!!! note "Under development"
    These DOCS are under development 

## What it does

Given a set of coupled-cluster (or DFT) energies and vibrational frequencies for
each molecule, plus a list of candidate reactions, ARCS:

1. Computes the **Gibbs free energy of reaction** and **equilibrium constant** for
   every reaction at a chosen temperature and pressure.
2. Builds a weighted, directed **reaction graph** (`networkx.MultiDiGraph`) where
   compounds and reactions are nodes, and edge weights come from a thermodynamic
   cost function.
3. **Samples** the graph many times over, starting from an initial gas mixture,
   using equilibrium chemistry at each step to update concentrations.
4. **Analyses** the ensemble of runs to rank reactions and summarise how each
   compound's concentration changes.

## The pipeline at a glance

``` mermaid
graph LR
    A[Quantum data<br/>energies + vibrations] --> B[GraphGenerator]
    B --> C[Reaction graph<br/>MultiDiGraph]
    C --> D[Traversal.sample]
    E[GenerateInitialConcentrations] --> D
    D --> F[AnalyseSampling]
    F --> G[Reaction statistics<br/>+ average concentrations]
```

## Where to go next

<!-- This text is hidden in the rendered output 
<div class="grid cards" markdown>

-   :material-download: **[Installation](installation.md)**

    Install ARCS and its scientific dependencies.

-   :material-rocket-launch: **[Quickstart](quickstart.md)**

    Run an end-to-end simulation in a few lines.

-   :material-book-open-variant: **[User Guide](guide/how-it-works.md)**

    Understand the model, its parameters, and how to tune them.

-   :material-code-braces: **[API Reference](api/generate.md)**

    Full documentation of every module, class, and function.

</div>
-->
## Citation and credits

ARCS was funded by and carried out in collaboration with
[Equinor](https://www.equinor.com/no),
[TotalEnergies](https://ts.totalenergies.com), and
[Shell](https://www.shell.com).

The author is Benjamin A. D. Williamson (`benjamin.williamson@ntnu.no`). ARCS is
released under the MIT licence.
