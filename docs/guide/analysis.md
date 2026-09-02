# Analysing results
<!--
[`AnalyseSampling`][arcs.analysis.AnalyseSampling] turns the raw list of walks
returned by [`Traversal.sample`][arcs.traversal.Traversal.sample] into readable
statistics and visualisations.

```python
from arcs.analysis import AnalyseSampling

analysis = AnalyseSampling()
```

## Formatting options

The constructor accepts two display flags:

| Flag | Effect |
| --- | --- |
| `use_markdown` | Render chemical formulae with `<sub>` subscripts (for notebooks / HTML). |
| `use_latex` | Render formulae with LaTeX subscripts instead. |

If both are set, `use_markdown` wins (with a warning).

## Reaction statistics

[`reaction_statistics`][arcs.analysis.AnalyseSampling.reaction_statistics] counts
how often each reaction appears across all samples and returns a dictionary of
reaction → frequency.

```python
import pandas as pd

stats = pd.Series(analysis.reaction_statistics(results)).sort_values(ascending=False)
stats.head(10)
```

By default (`flip_reaction=True`) reactions are oriented so that the thermodynamically
favoured direction is reported, which merges a reaction and its reverse into a
single consistent entry.

## Average concentration changes

[`average_sampling`][arcs.analysis.AnalyseSampling.average_sampling] summarises how
every species moved from its initial value, returning the initial value, mean,
difference, standard error, standard deviation, and variance:

```python
average_data = pd.DataFrame(analysis.average_sampling(results))
average_data = average_data.loc[~(average_data == 0).all(axis=1)]
average_data.sort_values(by="diff").round(2)
```

A negative `diff` means the species was consumed on average; a positive `diff`
means it was produced.

## Path-length diagnostics

Two helpers describe how far the walks progressed:

- [`count_path_length`][arcs.analysis.AnalyseSampling.count_path_length] — a
  histogram of how many steps each walk took.
- [`reduce_data_by_minimum_path_length`][arcs.analysis.AnalyseSampling.reduce_data_by_minimum_path_length]
  — filter the dataset to walks that ran for at least a given number of steps,
  which is handy for excluding walks that stalled immediately.

```python
lengths = analysis.count_path_length(results)
long_runs = analysis.reduce_data_by_minimum_path_length(results, minimum_path_length=3)
```

## Interactive network visualisation

[`result_to_pyvis`][arcs.analysis.AnalyseSampling.result_to_pyvis] builds a
[pyvis](https://pyvis.readthedocs.io) network of the most frequent reactions.
Compound nodes are boxes, reaction nodes are circles coloured by frequency using a
Matplotlib colormap.

```python
net = analysis.result_to_pyvis(results, head=10, cmap="Reds")
net.show("reactions.html")
```

`head` sets how many of the top reactions to include, and any extra keyword
arguments are forwarded to the pyvis `Network` constructor (height, width,
`notebook`, and so on).

See the [`arcs.analysis` API reference](../api/analysis.md) for exact signatures.
-->