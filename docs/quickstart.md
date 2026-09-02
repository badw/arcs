# Quickstart

This page walks through a complete ARCS run: building a reaction graph from the
bundled quantum data, seeding an initial gas mixture, sampling the network, and
analysing the outcome.

## 1. Build the reaction graph

[`GraphGenerator`][arcs.generate.GraphGenerator] reads the quantum data, evaluates
the thermodynamics of every reaction at your chosen conditions, and returns a
weighted `networkx.MultiDiGraph`.

```python
from arcs.generate import GraphGenerator

gg = GraphGenerator()
graph = gg.from_file(
    filename="src/arcs/data/quantum_data.json.gz",
    temperature=248,          # Kelvin
    pressure=20,              # bar
    max_reaction_length=4,    # max reactants + products per reaction
)
```


## 2. Set the initial concentrations

[`GenerateInitialConcentrations`][arcs.generate.GenerateInitialConcentrations]
produces a concentration dictionary keyed by every compound in the graph.
`update_ic` starts from an all-zero baseline and overrides just the species you
name (values are in ppm-like units).

```python
from arcs.generate import GenerateInitialConcentrations

concentrations = GenerateInitialConcentrations(graph=graph).update_ic(
    {"H2O": 30, "O2": 10, "SO2": 10, "H2S": 10, "NO2": 10}
)
```

## 3. Sample the network

[`Traversal`][arcs.traversal.Traversal] performs many independent random walks
over the graph. Each walk repeatedly picks likely compounds, ranks the reactions
that connect them, solves the resulting equilibrium with `chempy`, and updates
the concentrations. Sampling is parallelised across CPUs.

```python
from arcs.traversal import Traversal

t = Traversal(graph=graph)

results = t.sample(
    initial_concentrations=concentrations,
    nsamples=1000,   # number of independent random walks
    ncpus=4,         # parallel workers
)
```

`results` is a list with one entry per sample. Each entry holds the
`concentrations` at every step and the `reaction_statistics` for the reactions
that fired.

## 4. Analyse the results

[`AnalyseSampling`][arcs.analysis.AnalyseSampling] aggregates the ensemble.

### Reaction frequencies

Which reactions occurred most often across all samples:

```python
import pandas as pd
from arcs.analysis import AnalyseSampling

analysis = AnalyseSampling()
stats = pd.Series(analysis.reaction_statistics(results)).sort_values(ascending=False)
stats.head(10)
```

```text
1 H2 + 1 SO2 = 1 O2 + 1 H2S              369
1 H2O + 1 SO2 = 1 H2SO3                  270
2 H2 + 1 O2 = 2 H2O                      227
3 O2 + 2 H2S = 2 H2O + 2 SO2             163
...
```

### Average concentration changes

How each species shifted, on average, from its initial value:

```python
average_data = pd.DataFrame(analysis.average_sampling(results))
average_data = average_data.loc[~(average_data == 0).all(axis=1)]
average_data.sort_values(by="diff", inplace=True)
average_data.round(2)
```

```text
compound  initial  mean  diff  sem   std   var
H2S        10.0    4.88 -5.12  0.10  4.75  22.53
NO2        10.0    6.19 -3.81  0.10  4.85  23.48
O2         10.0    6.24 -3.76  0.12  5.76  33.18
...
```

### Visualise the reaction network

`result_to_pyvis` builds an interactive [pyvis](https://pyvis.readthedocs.io)
network of the most frequent reactions, coloured by frequency:

```python
net = analysis.result_to_pyvis(results, head=10, cmap="Reds")
net.show("reactions.html")
```

This produces a standalone HTML file you can open in a browser. Here is a live
example of the kind of network it generates — drag the nodes, scroll to zoom, and
hover a reaction to see its equation:

<iframe src="../examples/example_pyvis_graph.html"
        style="width: 100%; height: 600px; float:none; border: 1px solid var(--md-default-fg-color--lightest); border-radius: 0.2rem;"
        title="Interactive ARCS reaction network"
        loading="lazy"></iframe>

        
## Next steps

- [How ARCS works](guide/how-it-works.md) explains the model behind each step.
- [Sampling the network](guide/sampling.md) covers the `Traversal` parameters that
  control how walks explore the graph.
- The [API Reference](api/generate.md) documents every option in detail.
