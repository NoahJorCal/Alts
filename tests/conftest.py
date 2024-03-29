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
            100,   # 3:  Limit size of groups
            1000,  # 4:  Limit size of the whole population
            1,     # 5:  Number of descendants per survivor
            0,     # 6:  Proportion of migration between groups in each generation
            0,     # 7:  Proportion of emigrants of the whole population
            0,     # 8:  Proportion of immigrants of the whole population
            None,  # 9:  Default phenotype for immigrants
            1,     # 10:  Mean life expectancy
            0,     # 11:  Standard deviation of the life expectancy
            1,     # 12: Mean survival probability assigned to each individual at the beginning of the generation
            0,     # 13: Standard deviation of the survival probability
            3,     # 14: Generations taken into account in individuals' ancestry
            None]  # 15: Output file name


@pytest.fixture()
def one_group_general_configuration():
    return [0,     # 0:  Generations in the simulation (First generation is number 0)
            1,     # 1:  Number of groups in the population
            10,    # 2:  Mean number of individuals inside each group
            10,    # 3:  Limit size of groups
            10,    # 4:  Limit size of the whole population
            1,     # 5:  Number of descendants per survivor
            0,     # 6:  Proportion of migration between groups in each generation
            0,     # 7:  Proportion of emigrants of the whole population
            0,     # 8:  Proportion of immigrants of the whole population
            None,  # 9:  Default phenotype for immigrants
            100,   # 10:  Mean life expectancy
            0,     # 11:  Standard deviation of the life expectancy
            0.5,   # 12: Mean survival probability assigned to each individual at the beginning of the generation
            0,     # 13: Standard deviation of the survival probability
            3,     # 14: Generations taken into account in individuals' ancestry
            None]  # 15: Output file name


# ConfigParser object of configuration files is read in the same way as a dictionary
@pytest.fixture()
def ba_dom_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},     # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Locus' alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '0'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral',                 # Locus' alleles and inheritance pattern
                        'initial_frequencies': '1',           # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }


@pytest.fixture()
def ba_cod_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},     # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Locus' alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '0'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral',                 # Locus' alleles and inheritance pattern
                        'initial_frequencies': '1',           # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }


@pytest.fixture()
def ba_haplotypes_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},     # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Locus' alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '10',                 # Size of the loci for mutation and recombination
                          'mutation_rate': '0.1',             # Mutation probability on a position of the genome
                          'recombination_rate': '0.1'},       # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral',                 # Locus' alleles and inheritance pattern
                        'initial_frequencies': '1',           # Initial allele's frequencies in order of appearance
                        'locus_size': '10',                   # Size of the loci for mutation and recombination
                        'mutation_rate': '0.1',               # Mutation probability on a position of the genome
                        'recombination_rate': '0.1'}          # Recombination probability on a position of the genome
            }


@pytest.fixture()
def cod_mutate_genotype_model_config():
    return {'module': {'name': 'blind_altruism_genomes'},     # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Locus' alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '1'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral',                 # Locus' alleles and inheritance pattern
                        'initial_frequencies': '1',           # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }
