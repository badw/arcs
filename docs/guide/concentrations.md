# Initial concentrations
<!--
Every sampling run starts from a set of species concentrations. These are held in
a plain dictionary mapping each compound to a number (in ppm-like units), and
ARCS provides [`GenerateInitialConcentrations`][arcs.generate.GenerateInitialConcentrations]
to build one that is consistent with the graph you are sampling.

```python
from arcs.generate import GenerateInitialConcentrations

gic = GenerateInitialConcentrations(graph=graph)
```

The class inspects the graph so that the dictionary it returns covers exactly the
compound nodes present, which keeps the downstream equilibrium solver happy.

## Ways to build a starting mixture

### From an explicit mixture (most common)

[`update_ic`][arcs.generate.GenerateInitialConcentrations.update_ic] starts from an
all-zero baseline and sets only the species you specify:

```python
concentrations = gic.update_ic(
    {"H2O": 30, "O2": 10, "SO2": 10, "H2S": 10, "NO2": 10}
)
```

Everything not named stays at zero. This is the recommended way to define a
realistic impurity mixture.

### All zero

[`all_zero`][arcs.generate.GenerateInitialConcentrations.all_zero] returns a
dictionary with every compound set to `0` — a blank slate you can populate
yourself.

```python
concentrations = gic.all_zero()
```

### Fully random

[`all_random`][arcs.generate.GenerateInitialConcentrations.all_random] assigns a
small random concentration to every compound.

```python
concentrations = gic.all_random()
```

### Random over selected species

[`specific_random`][arcs.generate.GenerateInitialConcentrations.specific_random]
gives random concentrations to a chosen subset and zero to the rest:

```python
concentrations = gic.specific_random(compounds=["CO2", "H2O", "O2"])
```

## The `include_co2` flag

All four methods accept `include_co2` (default `True`). CO<sub>2</sub> is the bulk
background of the stored fluid, so during sampling it is treated as effectively
inexhaustible. Setting `include_co2=False` pins CO<sub>2</sub> to a reference value
of `1` instead of leaving it among the trace species.

See the [`arcs.generate` API reference](../api/generate.md#arcs.generate.GenerateInitialConcentrations)
for full signatures.
-->