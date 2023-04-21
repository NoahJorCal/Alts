import pytest
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Alts')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Alts', 'models')))


@pytest.fixture()
def base_general_configuration():
    return [0,     # 0:  Generations in the simulation (First generation is number 0)
            5,     # 1:  Number of groups in the population
            20,    # 2:  Mean number of individuals inside each group
            0,     # 3:  Standard deviation of number of individuals per group
            0,     # 4:  Proportion of migration per groups in each generation
            0,     # 5:  Proportion of emigrants of the whole population
            0,     # 6:  Proportion of immigrants of the whole population
            None,  # 7:  Default phenotype for immigrants
            1,     # 8:  Mean life expectancy
            0,     # 9:  Standard deviation of the life expectancy
            1,     # 10: Survival probability assigned to each individual at the beginning of the generation
            3]


@pytest.fixture()
def one_group_general_configuration():
    return [0,     # 0:  Generations in the simulation (First generation is number 0)
            1,     # 1:  Number of groups in the population
            10,    # 2:  Mean number of individuals inside each group
            0,     # 3:  Standard deviation of number of individuals per group
            0,     # 4:  Proportion of migration per groups in each generation
            0,     # 5:  Proportion of emigrants of the whole population
            0,     # 6:  Proportion of immigrants of the whole population
            None,  # 7:  Default phenotype for immigrants
            100,   # 8:  Mean life expectancy
            0,     # 9:  Standard deviation of the life expectancy
            0.5,   # 10: Survival probability assigned to each individual at the beginning of the generation
            3]


# ConfigParser object of configuration files is read in the same way as a dictionary
@pytest.fixture()
def ba_dom_model_configuration():
    return {'module': {'name': 'blind_altruism'},             # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'}               # Mutation rate
            }


@pytest.fixture()
def ba_cod_model_configuration():
    return {'module': {'name': 'blind_altruism'},             # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'}               # Mutation rate
            }


@pytest.fixture()
def sggb_dom_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'}               # Mutation rate
            }


@pytest.fixture()
def gb_dom_model_configuration():
    return {'module': {'name': 'green_beard'},
            'behaviour': {'alleles': 'selfish > altruistic',
                          'initial_frequencies': '0.5, 0.5',
                          'mutation_rate': '1'},
            'beard': {'alleles': 'nonbeard > beard',
                      'initial_frequencies': '0.5, 0.5',
                      'mutation_rate': '0'}}


@pytest.fixture()
def ba_mutate_genotype_model_config():
    return {'module': {'name': 'blind_altruism'},
            'behaviour': {'alleles': 'selfish > altruistic',
                          'initial_frequencies': '0.5, 0.5',
                          'mutation_rate': '1'}}


@pytest.fixture()
def sggb_mutate_genotype_model_config():
    return {'module': {'name': 'green_beard'},
            'behaviour': {'alleles': 'selfish > altruistic',
                          'initial_frequencies': '0.5, 0.5',
                          'mutation_rate': '1'},
            'beard': {'alleles': 'nonbeard > beard',
                      'initial_frequencies': '0.5, 0.5',
                      'mutation_rate': '0'}}


@pytest.fixture()
def gb_mutate_genotype_model_config():
    return {'module': {'name': 'green_beard'},
            'behaviour': {'alleles': 'selfish > altruistic',
                          'initial_frequencies': '0.5, 0.5',
                          'mutation_rate': '1'},
            'beard': {'alleles': 'nonbeard > beard',
                      'initial_frequencies': '0.5, 0.5',
                      'mutation_rate': '0'}}

# [0, 1, 20, 0, 0, 0, 0, None, 1, 0, 1]
# {'module': {'name': 'blind_altruism'}, 'behaviour': {'alleles': 'selfish ? altruistic', 'initial_frequencies': '1, 0', 'mutation_rate': '0'}}