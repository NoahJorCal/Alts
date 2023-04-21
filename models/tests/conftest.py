import pytest
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'models')))


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
            3,     # 11: Generations taking into account in individuals' ancestry
            10]    # 12: Size of the genome


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
            3,     # 11: Generations taking into account in individuals' ancestry
            10]    # 12: Size of the genome


# ConfigParser object of configuration files is read in the same way as a dictionary
@pytest.fixture()
def ba_dom_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},     # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Locus's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '0'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral > neutral',       # Locus's alleles and inheritance pattern
                        'initial_frequencies': '0.5, 0.5',    # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }


@pytest.fixture()
def ba_cod_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},     # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Locus's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '0'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral > neutral',       # Locus's alleles and inheritance pattern
                        'initial_frequencies': '0.5, 0.5',    # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }


@pytest.fixture()
def ba_mutate_genotype_model_config():
    return {'module': {'name': 'blind_altruism_genomes'},     # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Locus's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '1'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral > neutral',       # Locus's alleles and inheritance pattern
                        'initial_frequencies': '0.5, 0.5',    # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }
