import random
import pytest
import sys
from os import path
# Add the alts directory to the Python path
sys.path.append(path.join(path.dirname(__file__), '..', '..'))
import numpy as np
from random import randint
from simulator import Simulation
from simulator import Individual
from models import blind_altruism
from models import single_gene_green_beard
from models import green_beard


@pytest.fixture(scope='function')
def simulation(base_general_configuration, ba_dom_model_configuration):
    sim = Simulation(*base_general_configuration, ba_dom_model_configuration)
    sim.populate_groups()
    return sim


@pytest.fixture(scope='function')
def simulation_sggb(base_general_configuration, sggb_dom_model_configuration):
    sim = Simulation(*base_general_configuration, sggb_dom_model_configuration)
    sim.populate_groups()
    return sim


@pytest.fixture(scope='function')
def simulation_gb(base_general_configuration, gb_dom_model_configuration):
    sim = Simulation(*base_general_configuration, gb_dom_model_configuration)
    sim.populate_groups()
    return sim


@pytest.fixture()
def possible_recipients(simulation):
    return blind_altruism.get_possible_recipients(simulation.groups[0][0], simulation.groups[0])


@pytest.fixture()
def altruist(simulation):
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'altruistic']]
    return individual


@pytest.fixture()
def recipient(simulation):
    individual = Individual(simulation)
    return individual


def test_get_possible_recipients_ba(simulation):
    original_group = simulation.groups[0]
    original_group_size = len(original_group)
    # A random individual will be selected, to test the function, this individual does not need to be altruistic
    random_ind = original_group[randint(0, len(original_group))]
    possible_recipients = blind_altruism.get_possible_recipients(random_ind, original_group)
    group_size = len(possible_recipients)
    in_group = random_ind in possible_recipients
    assert group_size == original_group_size - 1
    assert not in_group


def test_get_possible_recipients_sggb(simulation_sggb):
    original_group = simulation_sggb.groups[0]
    altruist_count = 0
    altruist_i = None
    for ind in original_group:
        print(ind.id, ind.phenotype)
        if ind.phenotype[0] == 'altruistic':
            altruist_count += 1
            altruist_i = original_group.index(ind)
    if not altruist_i:
        altruist = original_group[0]
    else:
        altruist_count -= 1
        altruist = original_group[altruist_i]
    possible_recipients = single_gene_green_beard.get_possible_recipients(altruist, original_group)
    group_size = len(possible_recipients)
    in_group = altruist in possible_recipients
    assert group_size == altruist_count
    assert not in_group


def test_get_possible_recipients_gb(simulation_gb):
    original_group = simulation_gb.groups[0]
    beard_count = 0
    altruist_i = None
    for ind in original_group:
        if ind.phenotype[1] == 'beard':
            beard_count += 1
        if ind.phenotype[0] == 'altruistic':
            altruist_i = original_group.index(ind)
    if not altruist_i:
        altruist = original_group[0]
        if altruist.phenotype[1] == 'beard':
            beard_count -= 1
    else:
        altruist = original_group[altruist_i]
    if altruist.phenotype[1] == 'beard':
        beard_count -= 1
    print(altruist.id)
    possible_recipients = green_beard.get_possible_recipients(altruist, original_group)
    group_size = len(possible_recipients)
    in_group = altruist in possible_recipients
    assert group_size == beard_count
    assert not in_group


# The rest of the functions are common for every model


@pytest.mark.parametrize('altruist_ancestry, recipient_ancestry, penetrance,'
                         'benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost,'
                         'expected_altruist_survival_prob, expected_recipient_survival_prob', [
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           1,       # benefit_relatedness_ratio
                           1,       # cost_benefit_ratio
                           0.2,     # minimum_benefit
                           0.5,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.8,     # expected_altruist_survival_prob
                           1.2),    # expected_recipient_survival_prob
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           0,       # benefit_relatedness_ratio
                           0.5,     # cost_benefit_ratio
                           0.4,     # minimum_benefit
                           0.4,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.8,     # expected_altruist_survival_prob
                           1.4),    # expected_recipient_survival_prob
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           0.1,     # penetrance
                           0.5,     # benefit_relatedness_ratio
                           0.8,     # cost_benefit_ratio
                           1,       # minimum_benefit
                           0.8,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.92,    # expected_altruist_survival_prob
                           1.1),    # expected_recipient_survival_prob
                          ([[2], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[3], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           1,       # benefit_relatedness_ratio
                           1,       # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.4,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.75,    # expected_altruist_survival_prob
                           1.25),   # expected_recipient_survival_prob
                          ([[2], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[3], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           0.5,     # penetrance
                           1,       # benefit_relatedness_ratio
                           1,       # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.5,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.875,   # expected_altruist_survival_prob
                           1.125),  # expected_recipient_survival_prob
                          ([[2], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[3], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           0.5,     # benefit_relatedness_ratio
                           1,       # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.4,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.875,   # expected_altruist_survival_prob
                           1.125),  # expected_recipient_survival_prob
                          ([[2], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[3], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           1,       # benefit_relatedness_ratio
                           0.5,     # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.5,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.875,   # expected_altruist_survival_prob
                           1.25),   # expected_recipient_survival_prob
                          ([[3], [0, 0], [0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[4], [0, 0], [0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           1,       # benefit_relatedness_ratio
                           1,       # cost_benefit_ratio
                           0.1,     # minimum_benefit
                           0.5,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.9,     # expected_altruist_survival_prob
                           1.1),    # expected_recipient_survival_prob
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           1,       # benefit_relatedness_ratio
                           1,       # cost_benefit_ratio
                           0.5,     # minimum_benefit
                           0.1,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.9,     # expected_altruist_survival_prob
                           1.1)     # expected_recipient_survival_prob
                         ])
def test_altruistic_act(altruist_ancestry, recipient_ancestry, penetrance,
                           benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost,
                           expected_altruist_survival_prob, expected_recipient_survival_prob,
                           altruist, recipient):
    altruist.id = altruist_ancestry[0]
    altruist.ancestry = altruist_ancestry[1:]
    recipient.id = recipient_ancestry[0]
    recipient.ancestry = recipient_ancestry[1:]
    blind_altruism.altruistic_act(altruist, [recipient], penetrance, 0,
                                  benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost)
    assert altruist.survival_probability == pytest.approx(expected_altruist_survival_prob)
    assert recipient.survival_probability == pytest.approx(expected_recipient_survival_prob)


@pytest.mark.parametrize('benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost, expected_sp', [
    (0.5, 1, 0.1, 0.4, [0.45, 0.1, 0.5, 0.725, 0.45, 0.55, 0.5, 0.1, 0.7, 0.925]),
    (0.75, 0.5, 0.2, 0.3, [0.775, 0.2, 0.5, 0.9, 0.775, 0.525, 0.5, 0.2, 0.7, 1.125])
])
def test_selection_dom(simulation, benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost,
                          expected_sp, one_group_general_configuration, ba_dom_model_configuration):
    random.seed(190223)
    simulation = Simulation(*one_group_general_configuration, ba_dom_model_configuration)
    simulation.populate_groups()
    for i in range(3):
        simulation.selection_event()
        simulation.reproduce()
        simulation.generation += 1
    blind_altruism.selection(simulation.groups,
                             benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost)
    survival_probabilities = []
    for ind in simulation.groups[0]:
        survival_probabilities.append(ind.survival_probability)
    print(survival_probabilities)
    assert survival_probabilities == pytest.approx(expected_sp)


@pytest.mark.parametrize('benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost, expected_sp', [
    (0.5, 1, 0.1, 0.2, [0.56875, 0.45, 0.6625, 0.4, 0.4875, 0.55, 0.54375, 0.4, 0.4375, 0.5]),
    (0.75, 0.5, 0.1, 0.3, [0.35, 0.509375, 0.846875, 0.715625, 0.4625, 0.846875, 0.490625, 0.490625, 0.8375, 0.5]),
])
def test_selection_cod(simulation, benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost,
                          expected_sp, one_group_general_configuration, ba_cod_model_configuration):
    random.seed(950603)
    simulation = Simulation(*one_group_general_configuration, ba_cod_model_configuration)
    simulation.populate_groups()
    for i in range(3):
        simulation.selection_event()
        simulation.reproduce()
        simulation.generation += 1
    blind_altruism.selection(simulation.groups,
                             benefit_relatedness_ratio, cost_benefit_ratio, minimum_benefit, maximum_cost)
    survival_probabilities = []
    for ind in simulation.groups[0]:
        survival_probabilities.append(ind.survival_probability)
    assert survival_probabilities == pytest.approx(expected_sp)


