from csv import Error
import numpy as np
import json
from monty.serialization import loadfn
from arcs.generate import GetEnergyandVibrationsVASP
from arcs.generate import ReactionGibbsandEquilibrium
from arcs.generate import GraphGenerator
from arcs.generate import GenerateInitialConcentrations
from arcs.generate import parse_molecule
from ase.io import read
from ase.atoms import Atoms
import pytest


# GetEnergyandVibrationsVASP
def test_parse_molecule():
    assert parse_molecule("H2O") == {'H': 2, 'O': 1}
    assert parse_molecule("(CH3_2CHCH2NH2") == {'C': 4, 'H': 11, 'N': 1}


def test_get_initial_magnetic_moments():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    aseatoms = read('./vaspdata/relax/POSCAR')
    magmoms = gevv.get_initial_magnetic_moments(aseatoms)
    assert magmoms == [0, 0, 0]


def test_get_atoms():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    aseatoms = gevv.get_atoms()

    testaseatoms = Atoms.fromdict({'numbers': np.array([6, 8, 8]),
                                   'positions': np.array([[5.00000001, 5.00000002, 5.00000001],
                                                          [5.00000014, 5.00000002,
                                                           6.16346205],
                                                          [4.99999985, 4.99999996, 3.83653794]]),
                                   'initial_magmoms': np.array([0., 0., 0.]),
                                   'cell': np.array([[10.,  0.,  0.],
                                                     [0., 10.,  0.],
                                                     [0.,  0., 10.]]),
                                   'pbc': np.array([True,  True,  True])})

    assert list(aseatoms.todict()['numbers']) == list(
        testaseatoms.todict()['numbers'])
    # can add more here


def test_get_energy():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    energy = gevv.get_energy()

    assert energy == -28.38822683


def test_get_spin():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    spin = gevv.get_spin()

    assert spin == 0


def test_get_pointgroup():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    pointgroup = gevv.get_pointgroup()

    assert pointgroup == 'D*h'


def test_islinear():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    linear = gevv.islinear()

    assert linear == 'linear'


def test_get_rotation_num():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    rotation_num = gevv.get_rotation_num()

    assert rotation_num == 2


def test_get_vibrations():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    vibrations = gevv.get_vibrations()

    assert vibrations == [0.304843445,
                          0.172210317,
                          0.08320696,
                          0.082896424,
                          0.001993745,
                          0.000983726,
                          -7.8669e-05,
                          -0.000301184,
                          -0.001327542]


def test_as_dict():
    gevv = GetEnergyandVibrationsVASP(
        relax_directory='./vaspdata/relax/', vibrations_directory='./vaspdata/vibrations/')
    _dict = gevv.as_dict()

    # add an assert

# ReactionGibbsandEquilibrium


def test_json_data_loader():
    json.load(open('test_quantum_data.json'))


def test_Gibbs():
    quantum_data = loadfn('test_quantum_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=quantum_data)
    gibbs_free_energy = rge.Gibbs(
        quantum_dict=quantum_data['CO2'], temperature=100, pressure=1)
    assert gibbs_free_energy == pytest.approx(-5119.272188209697)


def test_reaction_gibbs():
    quantum_data = loadfn('test_quantum_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=quantum_data)
    reaction = quantum_data['reactions'][str(307)]
    reaction_gibbs = rge.reaction_gibbs(
        reaction=reaction,
        pressure=1,
        temperature=100
    )
    assert reaction_gibbs == pytest.approx(0.41032333892781025)


def test_equilibrium_constant():
    quantum_data = loadfn('test_quantum_data.json')
    rge = ReactionGibbsandEquilibrium(reaction_input=quantum_data)
    equilibrium_constant = rge.equilibrium_constant(gibbs_free_energy=-0.6184489112909475,
                                                    temperature=100)
    assert equilibrium_constant == 1.4738501117040624e+31


def test_get_reaction_gibbs_and_equilibrium():
    quantum_data = loadfn('test_quantum_data.json')
    reaction = quantum_data['reactions'][str(307)]
    rge = ReactionGibbsandEquilibrium(reaction_input=quantum_data)
    _dict = rge.get_reaction_gibbs_and_equilibrium(
        reaction=reaction, temperature=100, pressure=1)

    assert _dict == {'g': pytest.approx(
        0.41032333892781025), 'k': pytest.approx(2.0922514397133124e-21)}


#  test GraphGenerator

def test_cost_function():
    applied_reactions = [{'g': 0.6184489112909475,
                          'k': 6.784950464493314e-32,
                          'r': {'reaction_string': '1 CO2 + 1 H2 = 1 H2O + 1 CO',
                                'reactants': {'CO2': 1, 'H2': 1},
                                'products': {'H2O': 1, 'CO': 1}}}]
    gg = GraphGenerator()

    cost = gg.cost_function(
        gibbs_free_energy=applied_reactions[0]['g'],
        temperature=100,
        reactants=applied_reactions[0]['r']['reactants']
    )

    assert cost == pytest.approx(1.4084092098796381)


def test_generate_multidigraph():
    applied_reactions = [{'g': 4.747780176970167,
                          'k': 5.277282624067824e-240,
                          'r': {'reaction_string': '2 H2O = 1 O2 + 2 H2',
                                'reactants': {'H2O': 2},
                                'products': {'O2': 1, 'H2': 2}}},
                         {'g': 0.22399682291949574,
                          'k': 5.141111023488729e-12,
                          'r': {'reaction_string': '1 CO2 + 1 CH4 = 1 CH3COOH',
                                'reactants': {'CO2': 1, 'CH4': 1},
                                'products': {'CH3COOH': 1}}},
                         {'g': 0.7892106026904457,
                          'k': 1.6808780628806743e-40,
                          'r': {'reaction_string': '1 CO + 1 CH3OH = 2 H2CO',
                                'reactants': {'CO': 1, 'CH3OH': 1},
                                'products': {'H2CO': 2}}},
                         {'g': -1.2923694047965597,
                          'k': 1.356910917936724e+65,
                          'r': {'reaction_string': '1 O2 + 2 HNO2 = 2 HNO3',
                                'reactants': {'O2': 1, 'HNO2': 2},
                                'products': {'HNO3': 2}}},
                         {'g': -3.112355307239554,
                          'k': 7.174701288889934e+156,
                          'r': {'reaction_string': '1 O2 + 1 CH3COOH = 2 CH2O2',
                                'reactants': {'O2': 1, 'CH3COOH': 1},
                                'products': {'CH2O2': 2}}}]

    gg = GraphGenerator()
    graph = gg.generate_multidigraph(
        temperature=100, applied_reactions=applied_reactions)

    adjacency_map = {'H2O': {0: {0: {'weight': np.float64(1.9492014814446184)}}},
                     0: {'H2O': {0: {'weight': np.float64(0.8053103344762139)}},
                         'O2': {0: {'weight': np.float64(1.9492014814446184)}},
                         'H2': {0: {'weight': np.float64(1.9492014814446184)}}},
                     'O2': {0: {0: {'weight': np.float64(0.8053103344762139)}},
                            3: {0: {'weight': np.float64(1.2234902236516356)}},
                            4: {0: {'weight': np.float64(1.0985582213409335)}}},
                     'H2': {0: {0: {'weight': np.float64(0.8053103344762139)}}},
                     'CO2': {1: {0: {'weight': np.float64(1.3369778231103757)}}},
                     1: {'CO2': {0: {'weight': np.float64(1.2959924751414802)}},
                         'CH4': {0: {'weight': np.float64(1.2959924751414802)}},
                         'CH3COOH': {0: {'weight': np.float64(1.3369778231103757)}}},
                     'CH4': {1: {0: {'weight': np.float64(1.3369778231103757)}}},
                     'CH3COOH': {1: {0: {'weight': np.float64(1.2959924751414802)}},
                                 4: {0: {'weight': np.float64(1.0985582213409335)}}},
                     'CO': {2: {0: {'weight': np.float64(1.389551610603596)}}},
                     2: {'CO': {0: {'weight': np.float64(1.2451742189951533)}},
                         'CH3OH': {0: {'weight': np.float64(1.2451742189951533)}},
                         'H2CO': {0: {'weight': np.float64(1.389551610603596)}}},
                     'CH3OH': {2: {0: {'weight': np.float64(1.389551610603596)}}},
                     'H2CO': {2: {0: {'weight': np.float64(1.2451742189951533)}}},
                     3: {'O2': {0: {'weight': np.float64(1.4126027501246856)}},
                         'HNO2': {0: {'weight': np.float64(1.4126027501246856)}},
                         'HNO3': {0: {'weight': np.float64(1.2234902236516356)}}},
                     'HNO2': {3: {0: {'weight': np.float64(1.2234902236516356)}}},
                     'HNO3': {3: {0: {'weight': np.float64(1.4126027501246856)}}},
                     4: {'O2': {0: {'weight': np.float64(1.5532383356136172)}},
                         'CH3COOH': {0: {'weight': np.float64(1.5532383356136172)}},
                         'CH2O2': {0: {'weight': np.float64(1.0985582213409335)}}},
                     'CH2O2': {4: {0: {'weight': np.float64(1.5532383356136172)}}}}

    nodes = ['H2O',
             0,
             'O2',
             'H2',
             'CO2',
             1,
             'CH4',
             'CH3COOH',
             'CO',
             2,
             'CH3OH',
             'H2CO',
             3,
             'HNO2',
             'HNO3',
             4,
             'CH2O2']

    assert dict(graph.adjacency()) == adjacency_map

    assert list(graph.nodes()) == nodes

    def test_from_file():

        return (None)

    def test_from_dict():

        return (None)


# GenerateInitialConcentrations

def test_all_random():
    applied_reactions = [{'g': 4.747780176970167,
                          'k': 5.277282624067824e-240,
                          'r': {'reaction_string': '2 H2O = 1 O2 + 2 H2',
                                'reactants': {'H2O': 2},
                                'products': {'O2': 1, 'H2': 2}}},
                         {'g': 0.22399682291949574,
                          'k': 5.141111023488729e-12,
                          'r': {'reaction_string': '1 CO2 + 1 CH4 = 1 CH3COOH',
                                'reactants': {'CO2': 1, 'CH4': 1},
                                'products': {'CH3COOH': 1}}},
                         {'g': 0.7892106026904457,
                          'k': 1.6808780628806743e-40,
                          'r': {'reaction_string': '1 CO + 1 CH3OH = 2 H2CO',
                                'reactants': {'CO': 1, 'CH3OH': 1},
                                'products': {'H2CO': 2}}},
                         {'g': -1.2923694047965597,
                          'k': 1.356910917936724e+65,
                          'r': {'reaction_string': '1 O2 + 2 HNO2 = 2 HNO3',
                                'reactants': {'O2': 1, 'HNO2': 2},
                                'products': {'HNO3': 2}}},
                         {'g': -3.112355307239554,
                          'k': 7.174701288889934e+156,
                          'r': {'reaction_string': '1 O2 + 1 CH3COOH = 2 CH2O2',
                                'reactants': {'O2': 1, 'CH3COOH': 1},
                                'products': {'CH2O2': 2}}}]

    gg = GraphGenerator()
    graph = gg.generate_multidigraph(
        temperature=100, applied_reactions=applied_reactions)
    gic = GenerateInitialConcentrations(graph=graph)

    try:
        gic.all_random()
    except Exception:
        raise (Error)


def test_all_zero():
    applied_reactions = [{'g': 4.747780176970167,
                          'k': 5.277282624067824e-240,
                          'r': {'reaction_string': '2 H2O = 1 O2 + 2 H2',
                                'reactants': {'H2O': 2},
                                'products': {'O2': 1, 'H2': 2}}},
                         {'g': 0.22399682291949574,
                          'k': 5.141111023488729e-12,
                          'r': {'reaction_string': '1 CO2 + 1 CH4 = 1 CH3COOH',
                                'reactants': {'CO2': 1, 'CH4': 1},
                                'products': {'CH3COOH': 1}}},
                         {'g': 0.7892106026904457,
                          'k': 1.6808780628806743e-40,
                          'r': {'reaction_string': '1 CO + 1 CH3OH = 2 H2CO',
                                'reactants': {'CO': 1, 'CH3OH': 1},
                                'products': {'H2CO': 2}}},
                         {'g': -1.2923694047965597,
                          'k': 1.356910917936724e+65,
                          'r': {'reaction_string': '1 O2 + 2 HNO2 = 2 HNO3',
                                'reactants': {'O2': 1, 'HNO2': 2},
                                'products': {'HNO3': 2}}},
                         {'g': -3.112355307239554,
                          'k': 7.174701288889934e+156,
                          'r': {'reaction_string': '1 O2 + 1 CH3COOH = 2 CH2O2',
                                'reactants': {'O2': 1, 'CH3COOH': 1},
                                'products': {'CH2O2': 2}}}]

    gg = GraphGenerator()
    graph = gg.generate_multidigraph(
        temperature=100, applied_reactions=applied_reactions)
    gic = GenerateInitialConcentrations(graph=graph)
    all_zero = gic.all_zero()

    assert all_zero == {'H2O': 0,
                        'O2': 0,
                        'H2': 0,
                        'CO2': 0,
                        'CH4': 0,
                        'CH3COOH': 0,
                        'CO': 0,
                        'CH3OH': 0,
                        'H2CO': 0,
                        'HNO2': 0,
                        'HNO3': 0,
                        'CH2O2': 0}


def test_specific_random():
    applied_reactions = [{'g': 4.747780176970167,
                          'k': 5.277282624067824e-240,
                          'r': {'reaction_string': '2 H2O = 1 O2 + 2 H2',
                                'reactants': {'H2O': 2},
                                'products': {'O2': 1, 'H2': 2}}},
                         {'g': 0.22399682291949574,
                          'k': 5.141111023488729e-12,
                          'r': {'reaction_string': '1 CO2 + 1 CH4 = 1 CH3COOH',
                                'reactants': {'CO2': 1, 'CH4': 1},
                                'products': {'CH3COOH': 1}}},
                         {'g': 0.7892106026904457,
                          'k': 1.6808780628806743e-40,
                          'r': {'reaction_string': '1 CO + 1 CH3OH = 2 H2CO',
                                'reactants': {'CO': 1, 'CH3OH': 1},
                                'products': {'H2CO': 2}}},
                         {'g': -1.2923694047965597,
                          'k': 1.356910917936724e+65,
                          'r': {'reaction_string': '1 O2 + 2 HNO2 = 2 HNO3',
                                'reactants': {'O2': 1, 'HNO2': 2},
                                'products': {'HNO3': 2}}},
                         {'g': -3.112355307239554,
                          'k': 7.174701288889934e+156,
                          'r': {'reaction_string': '1 O2 + 1 CH3COOH = 2 CH2O2',
                                'reactants': {'O2': 1, 'CH3COOH': 1},
                                'products': {'CH2O2': 2}}}]

    gg = GraphGenerator()
    graph = gg.generate_multidigraph(
        temperature=100, applied_reactions=applied_reactions)
    gic = GenerateInitialConcentrations(graph=graph)
    specific_random = gic.specific_random(compounds=['H2O', 'CO'])

    assert specific_random['H2O'] > 0
    assert specific_random['CO'] > 0
    assert specific_random['O2'] == 0


def test_update_ic():
    applied_reactions = [{'g': 4.747780176970167,
                          'k': 5.277282624067824e-240,
                          'r': {'reaction_string': '2 H2O = 1 O2 + 2 H2',
                                'reactants': {'H2O': 2},
                                'products': {'O2': 1, 'H2': 2}}},
                         {'g': 0.22399682291949574,
                          'k': 5.141111023488729e-12,
                          'r': {'reaction_string': '1 CO2 + 1 CH4 = 1 CH3COOH',
                                'reactants': {'CO2': 1, 'CH4': 1},
                                'products': {'CH3COOH': 1}}},
                         {'g': 0.7892106026904457,
                          'k': 1.6808780628806743e-40,
                          'r': {'reaction_string': '1 CO + 1 CH3OH = 2 H2CO',
                                'reactants': {'CO': 1, 'CH3OH': 1},
                                'products': {'H2CO': 2}}},
                         {'g': -1.2923694047965597,
                          'k': 1.356910917936724e+65,
                          'r': {'reaction_string': '1 O2 + 2 HNO2 = 2 HNO3',
                                'reactants': {'O2': 1, 'HNO2': 2},
                                'products': {'HNO3': 2}}},
                         {'g': -3.112355307239554,
                          'k': 7.174701288889934e+156,
                          'r': {'reaction_string': '1 O2 + 1 CH3COOH = 2 CH2O2',
                                'reactants': {'O2': 1, 'CH3COOH': 1},
                                'products': {'CH2O2': 2}}}]

    gg = GraphGenerator()
    graph = gg.generate_multidigraph(
        temperature=100, applied_reactions=applied_reactions)
    gic = GenerateInitialConcentrations(graph=graph)
    specific_random = gic.update_ic(update_dict={'CO2': 1, 'H2O': 3})

    assert specific_random['H2O'] == 3
    assert specific_random['CO2'] == 1
    assert specific_random['O2'] == 0
