import pytest
import sys
from os import path
# Add the alts directory to the Python path
sys.path.append(path.join(path.dirname(__file__), '..'))
from simulator import Simulation
from simulator import Individual

### Individual object of simulator.py with blind_altruism model (ba) ###
''' Letter code of test fixtures' names
ba: blind_altruism
s: selfish
a: altruistic
n: nonbeard
b: beard
d: dominant
r: recessive
c: codominant
'''


# ConfigParser object of configuration files is read in the same way as a dictionary
@pytest.fixture(scope='function')
def sd_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},             # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '0'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral',       # Locus's alleles and inheritance pattern
                        'initial_frequencies': '1',    # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }


@pytest.fixture(scope='function')
def ad_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},             # Model name
            'behaviour': {'alleles': 'altruistic > selfish',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '0'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral',       # Locus's alleles and inheritance pattern
                        'initial_frequencies': '1',    # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }


@pytest.fixture(scope='function')
def c_model_configuration():
    return {'module': {'name': 'blind_altruism_genomes'},             # Model name
            'behaviour': {'alleles': 'altruistic = selfish',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'locus_size': '0',                  # Size of the loci for mutation and recombination
                          'mutation_rate': '0',               # Mutation probability on a position of the genome
                          'recombination_rate': '0'},         # Recombination probability on a position of the genome
            'neutral': {'alleles': 'neutral',       # Locus's alleles and inheritance pattern
                        'initial_frequencies': '1',    # Initial allele's frequencies in order of appearance
                        'locus_size': '0',                    # Size of the loci for mutation and recombination
                        'mutation_rate': '0',                 # Mutation probability on a position of the genome
                        'recombination_rate': '0'}            # Recombination probability on a position of the genome
            }


@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['neutral', 'neutral']],
     ['selfish', 'neutral']
     ),
    ([['selfish', 'altruistic'], ['neutral', 'neutral']],
     ['selfish', 'neutral']
     ),
    ([['altruistic', 'selfish'], ['neutral', 'neutral']],
     ['selfish', 'neutral']
     ),
    ([['altruistic', 'altruistic'], ['neutral', 'neutral']],
     ['altruistic', 'neutral']
     )
])
def test_generate_phenotype_sd_ba(base_general_configuration, sd_model_configuration, genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sd_model_configuration, True)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['neutral', 'neutral']],
     ['selfish', 'neutral']
     ),
    ([['selfish', 'altruistic'], ['neutral', 'neutral']],
     ['altruistic', 'neutral']
     ),
    ([['altruistic', 'selfish'], ['neutral', 'neutral']],
     ['altruistic', 'neutral']
     ),
    ([['altruistic', 'altruistic'], ['neutral', 'neutral']],
     ['altruistic', 'neutral']
     )
])
def test_generate_phenotype_ad_ba(base_general_configuration, ad_model_configuration, genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, ad_model_configuration, True)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['neutral', 'neutral']],
     ['selfish', 'neutral']
     ),
    ([['selfish', 'altruistic'], ['neutral', 'neutral']],
     ['altruistic_selfish', 'neutral']
     ),
    ([['altruistic', 'selfish'], ['neutral', 'neutral']],
     ['altruistic_selfish', 'neutral']
     ),
    ([['altruistic', 'altruistic'], ['neutral', 'neutral']],
     ['altruistic', 'neutral']
     )
])
def test_generate_phenotype_c_ba(base_general_configuration, c_model_configuration, genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, c_model_configuration, True)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


@pytest.mark.parametrize('genotype', [
    ([['altruistic', 'altruistic'], ['neutral', 'neutral']])
])
def test_generate_genome(genotype, base_general_configuration, ba_dom_model_configuration):
    simulation = Simulation(*base_general_configuration, ba_dom_model_configuration, True)
    individual = Individual(simulation)
    individual.genotype = genotype
    individual.generate_genome()
    result_genotype = []
    empty_snvs = True
    for locus in individual.genome.loci:
        region = []
        for chromosome in locus.chromosomes:
            region.append(chromosome.allele)
            if chromosome.snvs.size > 0:
                empty_snvs = False
        result_genotype.append(region)
    assert result_genotype == genotype
    assert empty_snvs


def test_age_individual_ba(sd_model_configuration):
    general_configuration = [9, 1, 50, 50, 50, 0, 0, 0, 0, None, 8, 0, 1, 0, 3, None]
    simulation = Simulation(*general_configuration, sd_model_configuration, True)
    individual = Individual(simulation)
    ages = []
    for i in range(8):
        ages.append(individual.age)
        dead = individual.age_individual()
    assert ages == list(range(8))
    assert dead
