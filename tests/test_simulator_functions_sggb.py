import pytest
import sys
from os import path
# Add the alts directory to the Python path
sys.path.append(path.join(path.dirname(__file__), '..'))
from simulator import Simulation
from simulator import Individual
import random
import numpy as np


# [0        #0:  Generations in the simulation (First generation is number 0)
#  1,       #1:  Number of groups in the population
#  20,      #2:  Mean number of individuals inside each group
#  0,       #3:  Standard deviation of number of individual per group
#  0,       #4:  Proportion of migration per groups in each generation
#  0,       #5:  Proportion of emigrants of the whole population
#  0,       #6:  Proportion of immigrants of the whole population
#  None,    #7:  Default phenotype for immigrants
#  1,       #8:  Mean life expectancy
#  0,       #9:  Standard deviation of the life expectancy
#  1,       #10: Survival probability assigned to each individual at the beginning of the generation
#  3]       #11: Generations taking into account in genealogy


# Randomizes individual's genotype based on initial gene frequencies
@pytest.mark.parametrize('general_configuration, model_configuration, expected_result', [
    ([0, 1, 20, 0, 0, 0, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '1, 0',
                    'mutation_rate': '0'}},
     {"[['selfish', 'selfish']]": 1,
      "[['selfish', 'altruistic']]": 0,
      "[['altruistic', 'selfish']]": 0,
      "[['altruistic', 'altruistic']]": 0}
     ),
    ([0, 1, 20, 0, 0, 0, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'altruistic > selfish',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}},
     {"[['selfish', 'selfish']]": 0.25,
      "[['selfish', 'altruistic']]": 0.25,
      "[['altruistic', 'selfish']]": 0.25,
      "[['altruistic', 'altruistic']]": 0.25}
     ),
    ([0, 1, 20, 0, 0, 0, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish = altruistic',
                    'initial_frequencies': '0, 1',
                    'mutation_rate': '0'}},
     {"[['selfish', 'selfish']]": 0,
      "[['selfish', 'altruistic']]": 0,
      "[['altruistic', 'selfish']]": 0,
      "[['altruistic', 'altruistic']]": 1}
     )
])
def test_generate_individual_genotype_sggb(general_configuration, model_configuration, expected_result):
    # random.seed(829504)
    simulation = Simulation(*general_configuration, model_configuration)
    genotypes_dictionary = {"[['selfish', 'selfish']]": 0,
                            "[['selfish', 'altruistic']]": 0,
                            "[['altruistic', 'selfish']]": 0,
                            "[['altruistic', 'altruistic']]": 0}
    replicas = 5000
    for i in range(replicas):
        genotype = str(simulation.generate_individual_genotype())
        genotypes_dictionary[genotype] += 1
    # Transform counts to proportions
    genotypes_dictionary = {k: v / replicas for total in (sum(genotypes_dictionary.values()),)
                            for k, v in genotypes_dictionary.items()}
    assert genotypes_dictionary == pytest.approx(expected_result, rel=0.1)


@pytest.mark.parametrize('general_configuration, expected_mean, expected_sd', [
    ([0, 1000, 20, 1, 0, 1, 0, None, 1, 0, 1, 3],
     20, 1),
    ([0, 20, 5, 0, 0, 1, 0, None, 1, 0, 1, 3],
     5, 0),
])
def test_populate_groups_sggb(general_configuration, expected_mean, expected_sd, sggb_dom_model_configuration):
    simulation = Simulation(*general_configuration, sggb_dom_model_configuration)
    simulation.populate_groups()
    sizes = []
    for group in simulation.groups:
        sizes.append(len(group))
    mean = np.mean(sizes)
    sd = np.std(sizes)
    assert mean == pytest.approx(expected_mean, rel=0.15)
    assert sd == pytest.approx(expected_sd, rel=0.15)


@pytest.mark.parametrize('general_configuration, expected_deaths', [
    ([0, 100, 20, 0, 0, 1, 0, None, 1, 0, 1, 3], 0),
    ([0, 100, 20, 0, 0, 1, 0, None, 1, 0, 0.5, 3], 10),
    ([0, 100, 20, 0, 0, 1, 0, None, 1, 0, 0, 3], 20)
])
def test_selection_event_sggb(general_configuration, expected_deaths, sggb_dom_model_configuration):
    simulation = Simulation(*general_configuration, sggb_dom_model_configuration)
    simulation.populate_groups()
    initial_inds_per_group = []
    for group in simulation.groups:
        initial_inds_per_group.append(len(group))
    simulation.selection_event()
    survivors_per_group = []
    for group_i in range(len(simulation.groups)):
        survivors = initial_inds_per_group[group_i] - len(simulation.groups[group_i])
        survivors_per_group.append(survivors)
    result_deaths = np.round(np.mean(survivors_per_group))
    assert result_deaths == expected_deaths


# Generates an Individual object with the configured immigrant genotype
@pytest.mark.parametrize('general_configuration, model_configuration, expected_result', [
    ([0, 1, 20, 0, 0, 0, 1, 'altruistic', 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}},
     {"[['selfish', 'selfish']]": 0,
      "[['selfish', 'altruistic']]": 0,
      "[['altruistic', 'selfish']]": 0,
      "[['altruistic', 'altruistic']]": 1}
     ),
    ([0, 1, 20, 0, 0, 0, 1, 'selfish', 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}},
     {"[['selfish', 'selfish']]": 1,
      "[['selfish', 'altruistic']]": 0,
      "[['altruistic', 'selfish']]": 0,
      "[['altruistic', 'altruistic']]": 0}
     )
])
def test_generate_immigrant_sggb(general_configuration, model_configuration, expected_result):
    # random.seed(496296)
    simulation = Simulation(*general_configuration, model_configuration)
    replicas = 100
    genotypes_dictionary = {"[['selfish', 'selfish']]": 0,
                            "[['selfish', 'altruistic']]": 0,
                            "[['altruistic', 'selfish']]": 0,
                            "[['altruistic', 'altruistic']]": 0}
    for i in range(replicas):
        genotype = str(simulation.generate_immigrant().genotype)
        genotypes_dictionary[genotype] += 1
    # Transform counts to proportions
    genotypes_dictionary = {k: v / replicas for total in (sum(genotypes_dictionary.values()),)
                            for k, v in genotypes_dictionary.items()}
    assert genotypes_dictionary == expected_result


# The proportion of the population configured will emigrate out of the population
@pytest.mark.parametrize('general_configuration, model_configuration', [
    ([0, 1, 20, 0, 0, 1, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}}),
    ([0, 1, 20, 0, 0, 0.5, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}}),
    ([0, 1, 20, 0, 0, 0.1, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}}),
    ([0, 1, 20, 0, 0, 0, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}})
])
def test_emigration_sggb(general_configuration, model_configuration):
    # random.seed(588786)
    simulation = Simulation(*general_configuration, model_configuration)
    simulation.populate_groups()
    initial_population_size = 0
    for group in simulation.groups:
        initial_population_size += len(group)

    simulation.migration()

    result_population_size = 0
    for group in simulation.groups:
        result_population_size += len(group)

    expected_population_size = initial_population_size - initial_population_size * general_configuration[5]
    assert result_population_size == expected_population_size


# The proportion of the population configured will be immigrants with the correct phenotype
@pytest.mark.parametrize('general_configuration, model_configuration', [
    ([0, 1, 20, 0, 0, 0, 1, 'altruistic', 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}}),
    ([0, 1, 20, 0, 0, 0, 0.5, 'selfish', 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}}),
    ([0, 1, 20, 0, 0, 0, 0.1, 'altruistic', 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}}),
    ([0, 1, 20, 0, 0, 0, 0, 'selfish', 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}})
])
def test_immigration_sggb(general_configuration, model_configuration):
    # random.seed(313951)
    simulation = Simulation(*general_configuration, model_configuration)
    simulation.populate_groups()
    initial_phenotypes = {"['selfish']": 0,
                          "['altruistic']": 0}
    for group_i in range(len(simulation.groups)):
        simulation.groups[group_i] = simulation.groups[group_i][0:int(len(simulation.groups[group_i]) / 2)]
        for ind in simulation.groups[group_i]:
            initial_phenotypes[str(ind.phenotype)] += 1
    current_population_size = 0
    for count in initial_phenotypes.values():
        current_population_size += count
    full_population_size = general_configuration[2]
    immigrants = min(round(full_population_size * general_configuration[6]),
                     full_population_size - current_population_size)

    expected_phenotypes = initial_phenotypes
    expected_phenotypes[f"['{general_configuration[7]}']"] += immigrants

    simulation.migration()
    result_phenotypes = {"['selfish']": 0,
                         "['altruistic']": 0}
    for group in simulation.groups:
        for ind in group:
            result_phenotypes[str(ind.phenotype)] += 1
    assert result_phenotypes == expected_phenotypes


@pytest.mark.parametrize('gene_index, reproducer_allele, expected_allele', [
    (0, 'selfish', 'altruistic'),
    (0, 'altruistic', 'selfish')
])
def test_mutate_genotype_sggb(gene_index, reproducer_allele, expected_allele,
                              base_general_configuration, sggb_mutate_genotype_model_config):
    simulation = Simulation(*base_general_configuration, sggb_mutate_genotype_model_config)
    result_allele = simulation.mutate_genotype(reproducer_allele, gene_index)
    assert result_allele == expected_allele


@pytest.mark.parametrize('sire_genotype, dam_genotype, expected_genotype', [
    ([['altruistic', 'altruistic']], [['altruistic', 'altruistic']], [['altruistic', 'altruistic']]),
    ([['altruistic', 'altruistic']], [['selfish', 'selfish']], [['altruistic', 'selfish']]),
    ([['altruistic', 'selfish']], [['altruistic', 'selfish']], [['altruistic', 'selfish']]),
    ([['selfish', 'selfish']], [['selfish', 'selfish']], [['selfish', 'selfish']]),
])
def test_generate_offspring_genotype_sggb(sire_genotype, dam_genotype, expected_genotype,
                                          base_general_configuration, sggb_dom_model_configuration):
    np.random.seed(221317)
    random.seed(221317)
    simulation = Simulation(*base_general_configuration, sggb_dom_model_configuration)
    sire = Individual(simulation)
    sire.genotype = sire_genotype
    dam = Individual(simulation)
    dam.genotype = dam_genotype

    descendant = Individual(simulation)

    simulation.generate_offspring_genotype([sire, dam], descendant)
    assert descendant.genotype == expected_genotype


@pytest.mark.parametrize('sire_id, sire_ancestry, dam_id, dam_ancestry, expected_ancestry', [
    (1, [[0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],
     2, [[0, 0], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],
     [[1, 2], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]),
    (1, [[3, 4], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],
     2, [[3, 4], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],
     [[1, 2], [3, 4, 3, 4], [0, 0, 0, 0, 0, 0, 0, 0]]),
    (1, [[3, 0], [5, 6, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],
     2, [[4, 7], [0, 0, 8, 9], [0, 0, 0, 0, 0, 0, 0, 0]],
     [[1, 2], [3, 0, 4, 7], [5, 6, 0, 0, 0, 0, 8, 9]]),
])
def test_generate_offspring_ancestry_sggb(sire_id, sire_ancestry, dam_id, dam_ancestry, expected_ancestry,
                                          base_general_configuration, sggb_dom_model_configuration):
    simulation = Simulation(*base_general_configuration, sggb_dom_model_configuration)
    sire = Individual(simulation)
    sire.ancestry = sire_ancestry
    sire.id = sire_id
    dam = Individual(simulation)
    dam.ancestry = dam_ancestry
    dam.id = dam_id

    descendant = Individual(simulation)
    simulation.generate_offspring_ancestry([sire, dam], descendant)
    assert descendant.ancestry == expected_ancestry


def test_discard_old_individuals_sggb(sggb_dom_model_configuration):
    general_configuration = [0, 1, 5000, 0, 0, 0, 0, None, 20, 2, 1, 3]
    simulation = Simulation(*general_configuration, sggb_dom_model_configuration)
    simulation.populate_groups()
    gen_of_death = []
    for i in range(80):
        start_inds = len(simulation.groups[0])
        simulation.groups = simulation.discard_old_individuals()
        simulation.group_exchange()
        final_inds = len(simulation.groups[0])
        simulation.generation += 1
        for j in range(start_inds - final_inds):
            gen_of_death.append(simulation.current_generation)
    mean = np.mean(gen_of_death)
    sd = np.std(gen_of_death)
    assert mean == pytest.approx(20, abs=0.1)
    assert sd == pytest.approx(2, abs=0.1)


def test_reproduce_sggb(sggb_dom_model_configuration):
    # np.random.seed(310475)
    random.seed(310475)
    general_configuration = [0, 1, 10, 0, 0, 0, 0, None, 10, 0, 0.5, 3]
    simulation = Simulation(*general_configuration, sggb_dom_model_configuration)
    simulation.populate_groups()

    simulation.selection_event()
    # Last ind will die of old age but will leave offspring
    simulation.groups[0][-1].life_expectancy = 1
    simulation.reproduce()

    reproducers_ids = []
    for newborn in simulation.groups[0][7:]:
        reproducers_ids.append(newborn.ancestry[0])

    individuals = []
    for ind in simulation.groups[0]:
        individuals.append(ind.id)

    assert len(simulation.groups[0]) == 10
    assert reproducers_ids == [[9, 10], [3, 5], [10, 2]]
    assert individuals == [1, 2, 3, 5, 6, 8, 9, 11, 12, 13]


# The proportion of the population configured will emigrate out of the population
@pytest.mark.parametrize('general_configuration, model_configuration', [
    ([0, 5, 20, 5, 0.25, 0, 0, None, 1, 0, 1, 3],
     {'module': {'name': 'single_gene_green_beard'},
      'behaviour': {'alleles': 'selfish > altruistic',
                    'initial_frequencies': '0.5, 0.5',
                    'mutation_rate': '0'}})
])
def test_group_exchange_sggb(general_configuration, model_configuration):
    np.random.seed(162497)
    random.seed(162497)
    expected_groups = [[28, 44, 81, 79, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 17, 18, 19, 20, 21, 22],
                       [11, 23, 60, 91, 83, 24, 25, 27, 29, 30, 31, 34, 35, 36, 38, 39, 40, 41, 42, 43, 45, 46, 47, 48, 49],
                       [16, 32, 33, 100, 93, 98, 51, 53, 54, 55, 57, 58, 59, 61, 62, 63, 64, 66, 67, 69, 70],
                       [13, 10, 26, 56, 106, 71, 72, 73, 74, 75, 76, 78, 80, 82, 84, 85, 86, 87, 88, 89, 90],
                       [15, 50, 37, 52, 65, 68, 77, 92, 94, 95, 96, 97, 99, 101, 102, 103, 104, 105, 107, 108]]

    simulation = Simulation(*general_configuration, model_configuration)
    simulation.populate_groups()

    simulation.group_exchange()

    result_groups = []
    for group in simulation.groups:
        group_ids = []
        for ind in group:
            group_ids.append(ind.id)
        result_groups.append(group_ids)
    assert result_groups == expected_groups


def test_generations_sggb(sggb_dom_model_configuration):
    general_configuration = [9, 1, 50, 0, 0, 0, 0, None, 1, 0, 1, 3]
    simulation = Simulation(*general_configuration, sggb_dom_model_configuration)
    simulation.populate_groups()
    for i in range(simulation.generations):
        simulation.pass_generation()
    assert simulation.current_generation == 10
