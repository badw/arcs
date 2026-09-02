# Sampling the network
<!--
[`Traversal`][arcs.traversal.Traversal] explores the reaction graph with repeated
random walks and collects statistics. It is constructed from a graph and then run
with [`sample`][arcs.traversal.Traversal.sample].

```python
from arcs.traversal import Traversal

t = Traversal(graph=graph)
results = t.sample(
    initial_concentrations=concentrations,
    nsamples=1000,
    ncpus=4,
)
```

## What one walk does

A single walk is performed by
[`random_walk`][arcs.traversal.Traversal.random_walk]. For up to `max_steps`
steps it repeats:

1. **Pick compounds** — [`get_weighted_random_compounds`][arcs.traversal.Traversal.get_weighted_random_compounds]
   samples a few species with probability proportional to their current
   concentration, above a discovery threshold.
2. **Rank reactions** — [`get_weighted_reaction_rankings`][arcs.traversal.Traversal.get_weighted_reaction_rankings]
   finds shortest paths between pairs of those compounds (each path is a reaction),
   checks that the atoms balance, and scores them by edge weight and reaction size.
3. **Choose one** — [`choose_reaction`][arcs.traversal.Traversal.choose_reaction]
   selects a reaction probabilistically, favouring better-ranked ones.
4. **Equilibrate** — a `chempy` equilibrium system is built and solved, and the
   resulting concentrations feed into the next step.

If no valid reaction can be found, the walk stops early. The result records the
concentrations at every step and which reactions fired.

## Tuning parameters

These are set as attributes on the `Traversal` instance (via the constructor's
`**kws`, or by assignment). Defaults are shown.

| Parameter | Default | Effect |
| --- | --- | --- |
| `max_steps` | `5` | Maximum reactions per walk. |
| `max_compounds` | `5` | Upper bound on compounds considered each step. |
| `discovery_threshold` | `5` | Minimum probability (%) for a species to be eligible. |
| `maximum_reaction_number` | `10` | Cap on ranked candidate reactions kept per step. |
| `exclude_co2` | `True` | Ignore CO<sub>2</sub> when weighting compounds (it is background). |
| `ceiling` | `2000` | Percentage above the median at which a concentration is deemed abnormally large. |
| `scale_largest` | `10` | Factor by which such large concentrations are scaled down. |
| `rank_small_reactions_higher` | `True` | Prefer reactions with fewer coefficients / atoms. |
| `rank_by_number_of_atoms` | `True` | When ranking by size, use atom count rather than coefficient count. |
| `shortest_path_method` | `"Djikstra"` | Algorithm passed to NetworkX for shortest paths. |

You can override any of these at construction time:

```python
t = Traversal(
    graph=graph,
    max_steps=10,
    discovery_threshold=3,
    exclude_co2=True,
)
```

## Running many walks

[`sample`][arcs.traversal.Traversal.sample] runs `nsamples` independent walks in
parallel across `ncpus` workers (via `tqdm_pathos`). Pass progress-bar options
through `tqdm_kws`:

```python
results = t.sample(
    initial_concentrations=concentrations,
    nsamples=5000,
    ncpus=8,
    tqdm_kws={"desc": "sampling"},
)
```

The return value is a list with one dictionary per walk, each containing
`concentrations` and `reaction_statistics`. Feed this straight into
[`AnalyseSampling`](analysis.md).

!!! tip "How many samples?"
    Reaction frequencies and average concentrations stabilise as `nsamples` grows.
    Start around 1000 for exploration and increase until the ranked statistics stop
    shifting.

Full signatures are in the [`arcs.traversal` API reference](../api/traversal.md).
-->