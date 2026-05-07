from arcs.traversal import Traversal
from arcs.generate import GraphGenerator
import pytest
import numpy as np
import random


@pytest.fixture(autouse=True)
def set_random_seed():
    random.seed(42)
    np.random.seed(42)


def test_length_multiplier():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    t.length_multiplier(0)
    assert t.length_multiplier(0) == 6


def test_get_weighted_random_compounds():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    assert t.get_weighted_random_compounds(concentrations=ic) == [
        np.str_('H2O'), np.str_('O2')]


def test_check_reactant_atoms():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    weighted_random_compounds = t.get_weighted_random_compounds(
        concentrations=ic)
    assert t.check_reactant_atoms(
        0, weighted_random_compounds=weighted_random_compounds)


def test_scale_large_concentrations():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1000, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    assert t.scale_large_concentrations(
        ic) == {'H2O': 100.0, 'O2': 0.5, 'H2': 0.5, 'H2O2': 0.5}


def test_get_weighted_reaction_rankings():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    weighted_random_compounds = t.get_weighted_random_compounds(
        concentrations=ic)

    weighted_reaction_rankings = t.get_weighted_reaction_rankings(
        weighted_random_compounds)

    assert weighted_reaction_rankings == {0: pytest.approx(
        2.1714934121520555), 3: pytest.approx(6.251319894856691)}


def test_choose_reaction():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    weighted_random_compounds = t.get_weighted_random_compounds(
        concentrations=ic)

    weighted_reaction_rankings = t.get_weighted_reaction_rankings(
        weighted_random_compounds)

    assert t.choose_reaction(weighted_reaction_rankings) == 0


def test_generate_chempy_eqsystem():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    assert t.generate_chempy_eqsystem(
        0).string() == '2 H2 + O2 = 2 H2O; 8.03e+112\n'


def test_chempy_equilibrium_concentrations():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    chempy_eqsystem = t.generate_chempy_eqsystem(1)

    concs = t.chempy_equilibrium_concentrations(
        concentrations=ic, equilibrium_reaction=chempy_eqsystem, chempy_sane=True)

    assert concs == {'H2O': 1,
                     'O2': pytest.approx(3.105623900164973e-24),
                     'H2': pytest.approx(1.0381458774827786e-27),
                     'H2O2': pytest.approx(1.0)}

    chempy_eqsystem = t.generate_chempy_eqsystem(0)

    concs = t.chempy_equilibrium_concentrations(
        concentrations=ic, equilibrium_reaction=chempy_eqsystem, chempy_sane=True)

    assert concs is None

    concs = t.chempy_equilibrium_concentrations(
        concentrations=ic, equilibrium_reaction=chempy_eqsystem, chempy_sane=False)

    assert concs == {'H2O': pytest.approx(1.600000000010232),
                     'O2': pytest.approx(2.3699276358254308e-43),
                     'H2': pytest.approx(1.1599899546332952e-35),
                     'H2O2': 0.5}


def test_random_walk():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    assert t.random_walk(ic, chempy_sane=True) == {'concentrations': {0: {'H2O': 1, 'O2': 0.5, 'H2': 0.5, 'H2O2': 0.5},
                                                                      1: {'H2O': 1,
                                                                          'O2': pytest.approx(3.105623900164973e-24),
                                                                          'H2': pytest.approx(1.0381458774827786e-27),
                                                                          'H2O2': pytest.approx(1.0)},
                                                                      2: {'H2O': pytest.approx(1.0),
                                                                          'O2': pytest.approx(3.105623900164973e-24),
                                                                          'H2': pytest.approx(3.863637144171021e-63),
                                                                          'H2O2': pytest.approx(1.0)},
                                                                      3: {'H2O': pytest.approx(1.9999984518648553),
                                                                          'O2': pytest.approx(0.4999992259324276),
                                                                          'H2': pytest.approx(3.863637144171021e-63),
                                                                          'H2O2': pytest.approx(1.5481351446973076e-06)}},
                                                   'reaction_statistics': {0: None,
                                                                           1: {'reaction': {'reaction_string': '1 H2 + 1 O2 = 1 H2O2',
                                                                                            'reactants': {'H2': 1, 'O2': 1},
                                                                               'products': {'H2O2': 1}},
                                                                               'equilibrium_constant': pytest.approx(3.101649763848951e+50)},
                                                                           2: {'reaction': {'reaction_string': '1 H2 + 1 H2O2 = 2 H2O',
                                                                                            'reactants': {'H2': 1, 'H2O2': 1},
                                                                               'products': {'H2O': 2}},
                                                                               'equilibrium_constant': pytest.approx(2.5882347712405565e+62)},
                                                                           3: {'reaction': {'reaction_string': '1 O2 + 2 H2O = 2 H2O2',
                                                                                            'reactants': {'O2': 1, 'H2O': 2},
                                                                               'products': {'H2O2': 2}},
                                                                               'equilibrium_constant': pytest.approx(1.1983649235815957e-12)}}}
    # perhaps add chempy_sane=False?


def test_sampling_function():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    t.initial_concentrations = ic
    assert t.sampling_function(0) == {'concentrations': {0: {'H2O': 1, 'O2': 0.5, 'H2': 0.5, 'H2O2': 0.5},
                                                         1: {'H2O': 1,
                                                             'O2': pytest.approx(3.105623900164973e-24),
                                                             'H2': pytest.approx(1.0381458774827786e-27),
                                                             'H2O2': pytest.approx(1.0)},
                                                         2: {'H2O': pytest.approx(1.0),
                                                             'O2': pytest.approx(3.105623900164973e-24),
                                                             'H2': pytest.approx(3.863637144171021e-63),
                                                             'H2O2': pytest.approx(1.0)},
                                                         3: {'H2O': pytest.approx(1.9999984518648553),
                                                             'O2': pytest.approx(0.4999992259324276),
                                                             'H2': pytest.approx(3.863637144171021e-63),
                                                             'H2O2': pytest.approx(1.5481351446973076e-06)}},
                                      'reaction_statistics': {0: None,
                                                              1: {'reaction': {'reaction_string': '1 H2 + 1 O2 = 1 H2O2',
                                                                               'reactants': {'H2': 1, 'O2': 1},
                                                                  'products': {'H2O2': 1}},
                                                                  'equilibrium_constant': pytest.approx(3.101649763848951e+50)},
                                                              2: {'reaction': {'reaction_string': '1 H2 + 1 H2O2 = 2 H2O',
                                                                               'reactants': {'H2': 1, 'H2O2': 1},
                                                                  'products': {'H2O': 2}},
                                                                  'equilibrium_constant': pytest.approx(2.5882347712405565e+62)},
                                                              3: {'reaction': {'reaction_string': '1 O2 + 2 H2O = 2 H2O2',
                                                                               'reactants': {'O2': 1, 'H2O': 2},
                                                                  'products': {'H2O2': 2}},
                                                                  'equilibrium_constant': pytest.approx(1.1983649235815957e-12)}}}


def test_sample():
    gg = GraphGenerator()
    graph = gg.from_file(
        filename="./test_graph.json",
        temperature=248,
        pressure=20,
    )
    t = Traversal(graph=graph)
    ic = {"H2O": 1, "O2": 0.5, "H2": 0.5, "H2O2": 0.5}
    results = t.sample(initial_concentrations=ic,
                       ncpus=1,
                       nsamples=2,
                       tqdm_kws={"disable": True})

    assert results == [{'concentrations': {0: {'H2O': 1, 'O2': 0.5, 'H2': 0.5, 'H2O2': 0.5},
                                           1: {'H2O': 1,
                                               'O2': pytest.approx(3.105623900164973e-24),
                                               'H2': pytest.approx(1.0381458774827786e-27),
                                               'H2O2': pytest.approx(1.0)},
                                           2: {'H2O': pytest.approx(1.0),
                                               'O2': pytest.approx(3.105623900164973e-24),
                                               'H2': pytest.approx(3.863637144171021e-63),
                                               'H2O2': pytest.approx(1.0)},
                                           3: {'H2O': pytest.approx(1.9999984518648553),
                                               'O2': pytest.approx(0.4999992259324276),
                                               'H2': pytest.approx(3.863637144171021e-63),
                                               'H2O2': pytest.approx(1.5481351446973076e-06)}},
                        'reaction_statistics': {0: None,
                                                1: {'reaction': {'reaction_string': '1 H2 + 1 O2 = 1 H2O2',
                                                                 'reactants': {'H2': 1, 'O2': 1},
                                                                 'products': {'H2O2': 1}},
                                                    'equilibrium_constant': pytest.approx(3.101649763848951e+50)},
                                                2: {'reaction': {'reaction_string': '1 H2 + 1 H2O2 = 2 H2O',
                                                                 'reactants': {'H2': 1, 'H2O2': 1},
                                                    'products': {'H2O': 2}},
                                                    'equilibrium_constant': pytest.approx(2.5882347712405565e+62)},
                                                3: {'reaction': {'reaction_string': '1 O2 + 2 H2O = 2 H2O2',
                                                                 'reactants': {'O2': 1, 'H2O': 2},
                                                    'products': {'H2O2': 2}},
                                                    'equilibrium_constant': pytest.approx(1.1983649235815957e-12)}}},
                       {'concentrations': {0: {'H2O': 1, 'O2': 0.5, 'H2': 0.5, 'H2O2': 0.5},
                                           1: {'H2O': 1,
                                               'O2': pytest.approx(3.105623900164973e-24),
                                               'H2': pytest.approx(1.0381458774827786e-27),
                                               'H2O2': pytest.approx(1.0)},
                                           2: {'H2O': pytest.approx(1.0),
                                               'O2': pytest.approx(3.105623900164973e-24),
                                               'H2': pytest.approx(3.863637144171021e-63),
                                               'H2O2': pytest.approx(1.0)},
                                           3: {'H2O': pytest.approx(1.0),
                                               'O2': pytest.approx(3.105623900164973e-24),
                                               'H2': pytest.approx(3.863637144171021e-63),
                                               'H2O2': pytest.approx(1.0)},
                                           4: {'H2O': pytest.approx(1.0),
                                               'O2': pytest.approx(3.105623900164973e-24),
                                               'H2': pytest.approx(3.863637144171021e-63),
                                               'H2O2': pytest.approx(1.0)},
                                           5: {'H2O': pytest.approx(1.9999984518648553),
                                               'O2': pytest.approx(0.4999992259324276),
                                               'H2': pytest.approx(3.863637144171021e-63),
                                               'H2O2': pytest.approx(1.5481351446973076e-06)}},
                        'reaction_statistics': {0: None,
                                                1: {'reaction': {'reaction_string': '1 H2 + 1 O2 = 1 H2O2',
                                                                 'reactants': {'H2': 1, 'O2': 1},
                                                                 'products': {'H2O2': 1}},
                                                    'equilibrium_constant': pytest.approx(3.101649763848951e+50)},
                                                2: {'reaction': {'reaction_string': '1 H2 + 1 H2O2 = 2 H2O',
                                                                 'reactants': {'H2': 1, 'H2O2': 1},
                                                    'products': {'H2O': 2}},
                                                    'equilibrium_constant': pytest.approx(2.5882347712405565e+62)},
                                                3: {'reaction': {'reaction_string': '1 H2 + 1 H2O2 = 2 H2O',
                                                                 'reactants': {'H2': 1, 'H2O2': 1},
                                                    'products': {'H2O': 2}},
                                                    'equilibrium_constant': pytest.approx(2.5882347712405565e+62)},
                                                4: {'reaction': {'reaction_string': '1 H2 + 1 H2O2 = 2 H2O',
                                                                 'reactants': {'H2': 1, 'H2O2': 1},
                                                    'products': {'H2O': 2}},
                                                    'equilibrium_constant': pytest.approx(2.5882347712405565e+62)},
                                                5: {'reaction': {'reaction_string': '1 O2 + 2 H2O = 2 H2O2',
                                                                 'reactants': {'O2': 1, 'H2O': 2},
                                                    'products': {'H2O2': 2}},
                                                    'equilibrium_constant': pytest.approx(1.1983649235815957e-12)}}}]
