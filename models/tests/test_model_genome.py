import random

import pytest
from random import randint
import sys
from os import path
# Add the alts directory to the Python path
sys.path.append(path.join(path.dirname(__file__), '..'))
sys.path.append(path.join(path.dirname(__file__), '..', '..'))
from simulator import Simulation
from simulator import Individual
from models import blind_altruism_genomes


@pytest.fixture(scope='function')
def simulation(base_general_configuration, ba_dom_model_configuration):
    sim = Simulation(*base_general_configuration, ba_dom_model_configuration, True)
    sim.populate_groups()
    return sim


@pytest.fixture()
def possible_recipients(simulation):
    return blind_altruism_genomes.get_possible_recipients(simulation.groups[0][0], simulation.groups[0])


@pytest.fixture()
def altruist(simulation):
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'altruistic'], ['neutral', 'neutral']]
    individual.generate_genome()
    return individual


@pytest.fixture()
def recipient(simulation):
    individual = Individual(simulation)
    return individual


@pytest.fixture()
def dom_group(one_group_general_configuration, ba_dom_model_configuration):
    simulation = Simulation(*one_group_general_configuration, ba_dom_model_configuration, True)
    group = []
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'altruistic'], ['neutral', 'neutral']]
    individual.id = 1
    individual.ancestry = [[0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'selfish'], ['neutral', 'neutral']]
    individual.id = 6
    individual.ancestry = [[0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 11
    individual.ancestry = [[7, 6], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'altruistic'], ['neutral', 'neutral']]
    individual.id = 15
    individual.ancestry = [[9, 7], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'selfish'], ['neutral', 'neutral']]
    individual.id = 19
    individual.ancestry = [[1, 6], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'selfish'], ['neutral', 'neutral']]
    individual.id = 20
    individual.ancestry = [[1, 6], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 21
    individual.ancestry = [[11, 6], [7, 6, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'altruistic'], ['neutral', 'neutral']]
    individual.id = 22
    individual.ancestry = [[1, 15], [0, 0, 9, 7], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 23
    individual.ancestry = [[11, 15], [7, 6, 9, 7], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 24
    individual.ancestry = [[11, 6], [7, 6, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    return [group]


@pytest.fixture()
def cod_group(one_group_general_configuration, ba_cod_model_configuration):
    simulation = Simulation(*one_group_general_configuration, ba_cod_model_configuration, True)
    group = []
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'selfish'], ['neutral', 'neutral']]
    individual.id = 1
    individual.ancestry = [[0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 6
    individual.ancestry = [[0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 11
    individual.ancestry = [[7, 6], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 15
    individual.ancestry = [[9, 7], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 19
    individual.ancestry = [[1, 6], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 20
    individual.ancestry = [[1, 6], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 21
    individual.ancestry = [[11, 6], [7, 6, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['altruistic', 'selfish'], ['neutral', 'neutral']]
    individual.id = 22
    individual.ancestry = [[1, 15], [0, 0, 9, 7], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 23
    individual.ancestry = [[11, 15], [7, 6, 9, 7], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    individual = Individual(simulation)
    individual.genotype = [['selfish', 'selfish'], ['neutral', 'neutral']]
    individual.id = 24
    individual.ancestry = [[11, 6], [7, 6, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
    group.append(individual)
    return [group]


def test_get_possible_recipients(simulation):
    original_group = simulation.groups[0]
    original_group_size = len(original_group)
    # A random individual will be selected, to test the function, this individual does not need to be altruistic
    random_ind = original_group[randint(0, len(original_group) - 1)]
    possible_recipients = blind_altruism_genomes.get_possible_recipients(random_ind, original_group)
    group_size = len(possible_recipients)
    in_group = random_ind in possible_recipients
    assert group_size == original_group_size - 1
    assert not in_group


# The rest of the functions are common for every model


@pytest.mark.parametrize('altruist_ancestry, recipient_ancestry, penetrance,'
                         'exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,'
                         'expected_altruist_survival_prob, expected_recipient_survival_prob', [
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           1,       # exp_factor
                           1,       # cost_benefit_ratio
                           0.2,     # minimum_benefit
                           0.5,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.8,     # expected_altruist_survival_prob
                           1.2),    # expected_recipient_survival_prob
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           2,       # exp_factor
                           0.5,     # cost_benefit_ratio
                           0.4,     # minimum_benefit
                           0.4,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.8,     # expected_altruist_survival_prob
                           1.4),    # expected_recipient_survival_prob
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           0.1,     # penetrance
                           4,       # exp_factor
                           0.8,     # cost_benefit_ratio
                           1,       # minimum_benefit
                           0.8,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.92,    # expected_altruist_survival_prob
                           1.1),    # expected_recipient_survival_prob
                          ([[2], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[3], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           2,       # exp_factor
                           1,       # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.4,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.898463676,   # expected_altruist_survival_prob
                           1.101536324),  # expected_recipient_survival_prob
                          ([[2], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[3], [1, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           0.5,     # penetrance
                           4,       # exp_factor
                           1,       # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.5,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.983970698,   # expected_altruist_survival_prob
                           1.016029302),  # expected_recipient_survival_prob
                          ([[3], [1, 2], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[4], [1, 2], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           2,       # exp_factor
                           1,       # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.4,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.731058579,   # expected_altruist_survival_prob
                           1.268941421),  # expected_recipient_survival_prob
                          ([[3], [1, 2], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[4], [1, 2], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,       # penetrance
                           2,       # exp_factor
                           1,     # cost_benefit_ratio
                           0.01,    # minimum_benefit
                           0.25,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.75,   # expected_altruist_survival_prob
                           1.25),   # expected_recipient_survival_prob
                          ([[2], [0, 0], [0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[3], [0, 0], [0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           1,        # penetrance
                           4,        # exp_factor
                           0.6,      # cost_benefit_ratio
                           0.05,     # minimum_benefit
                           0.5,      # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.97,     # expected_altruist_survival_prob
                           1.05),    # expected_recipient_survival_prob
                          ([[1], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # altruist_ancestry
                           [[2], [0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],  # recipient_ancestry
                           0.5,     # penetrance
                           10,      # exp_factor
                           1,       # cost_benefit_ratio
                           0.5,     # minimum_benefit
                           0.5,     # maximum_cost
                           # Base survival probability from base_general_configuration is 1
                           0.75,    # expected_altruist_survival_prob
                           1.25)    # expected_recipient_survival_prob
                         ])
def test_altruistic_act(altruist_ancestry, recipient_ancestry, penetrance,
                        exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,
                        expected_altruist_survival_prob, expected_recipient_survival_prob,
                        altruist, recipient):
    altruist.id = altruist_ancestry[0]
    altruist.ancestry = altruist_ancestry[1:]
    recipient.id = recipient_ancestry[0]
    recipient.ancestry = recipient_ancestry[1:]
    blind_altruism_genomes.altruistic_act(altruist, [recipient], penetrance, 0,
                                          exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost)
    assert altruist.survival_probability == pytest.approx(expected_altruist_survival_prob)
    assert recipient.survival_probability == pytest.approx(expected_recipient_survival_prob)


@pytest.mark.parametrize('exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost, expected_sp', [
    (2, 1, 0.1, 0.4, [1, 1.1, 1, 1, 1.268941, 1, 1, 0.631059, 1, 1]),
    # (4, 0.5, 0.1, 0.1, [1.040398539, 1.1, 1, 1, 1.119202922, 1, 1, 0.9, 1, 1])
])
def test_selection_dom(simulation, exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,
                       expected_sp, dom_group):
    random.seed(190913)
    blind_altruism_genomes.selection(dom_group, exp_factor, cost_benefit_ratio,
                                     minimum_benefit, maximum_cost)
    survival_probabilities = []
    for ind in dom_group[0]:
        survival_probabilities.append(ind.survival_probability)
    print(survival_probabilities)
    print(expected_sp)
    assert survival_probabilities == pytest.approx(expected_sp)


@pytest.mark.parametrize('exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost, expected_sp', [
    (2, 1, 0.1, 0.4, [1, 1.05, 1, 1, 1.134470711, 1, 1, 0.815529289, 1, 1]),
    (4, 0.5, 0.1, 0.1, [1.020199269, 1.05, 1, 1, 1.059601461, 1, 1, 0.95, 1, 1])
])
def test_selection_cod(simulation, exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,
                       expected_sp, cod_group):
    random.seed(190913)
    blind_altruism_genomes.selection(cod_group,
                                     exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost)
    survival_probabilities = []
    for ind in cod_group[0]:
        survival_probabilities.append(ind.survival_probability)
    assert survival_probabilities == pytest.approx(expected_sp)
