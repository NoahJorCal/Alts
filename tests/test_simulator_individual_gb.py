import pytest
import sys
from os import path
# Add the alts directory to the Python path
sys.path.append(path.join(path.dirname(__file__), '..'))
from simulator import Simulation
from simulator import Individual

### Individual object of simulator.py with green_beard model (gb) ###
''' Letter code of test fixtures' names
gb: green_beard
s: selfish
a: altruistic
n: nonbeard
b: beard
d: dominant
r: recessive
c: codominant
'''


# ConfigParser object of configuration files is read in the same way as a dictionary
# Behaviour is fixed on selfish > altruistic
@pytest.fixture(scope='function')
def sd_nd_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'nonbeard > beard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


@pytest.fixture(scope='function')
def sd_bd_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'beard > nonbeard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


@pytest.fixture(scope='function')
def sd_nbc_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'nonbeard = beard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


# Behaviour is fixed on altruistic > selfish
@pytest.fixture(scope='function')
def ad_nd_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'altruistic > selfish',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'nonbeard > beard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


@pytest.fixture(scope='function')
def ad_bd_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'altruistic > selfish',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'beard > nonbeard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


@pytest.fixture(scope='function')
def ad_nbc_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'altruistic > selfish',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'nonbeard = beard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


# Behaviour is fixed on selfish = altruistic
@pytest.fixture(scope='function')
def sac_nd_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'nonbeard > beard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


@pytest.fixture(scope='function')
def sac_bd_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'beard > nonbeard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


@pytest.fixture(scope='function')
def sac_nbc_model_configuration():
    return {'module': {'name': 'single_gene_green_beard'},    # Model name
            'behaviour': {'alleles': 'selfish = altruistic',  # Behaviour gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'},              # Mutation rate
            'beard': {'alleles': 'nonbeard = beard',          # Beard gene's alleles and inheritance pattern
                      'initial_frequencies': '0.5, 0.5',      # Initial allele's frequencies in order of appearance
                      'mutation_rate': '0.0'}                 # Mutation rate
            }


### selfish dominant ###
# nonbeard dominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['selfish', 'nonbeard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_sd_nd_gb(base_general_configuration, sd_nd_model_configuration,
                                     genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sd_nd_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### selfish dominant ###
# beard dominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['selfish', 'beard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['selfish', 'beard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_sd_bd_gb(base_general_configuration, sd_bd_model_configuration,
                                     genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sd_bd_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### selfish dominant ###
# beard codominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['selfish', 'nonbeard_beard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['selfish', 'nonbeard_beard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_sd_nbc_gb(base_general_configuration, sd_nbc_model_configuration,
                                      genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sd_nbc_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### altruist dominant ###
# nonbeard dominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['altruistic', 'nonbeard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['altruistic', 'nonbeard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_ad_nd_gb(base_general_configuration, ad_nd_model_configuration,
                                     genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, ad_nd_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### altruist dominant ###
# beard dominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['altruistic', 'beard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['altruistic', 'beard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_ad_bd_gb(base_general_configuration, ad_bd_model_configuration,
                                     genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, ad_bd_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### altruist dominant ###
# beard codominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['altruistic', 'nonbeard_beard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['altruistic', 'nonbeard_beard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_ad_nbc_gb(base_general_configuration, ad_nbc_model_configuration,
                                      genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, ad_nbc_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### behaviour codominant ###
# nonbeard dominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['selfish_altruistic', 'nonbeard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['selfish_altruistic', 'nonbeard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_sac_nd_gb(base_general_configuration, sac_nd_model_configuration,
                                      genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sac_nd_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### behaviour codominant ###
# beard dominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['selfish_altruistic', 'beard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['selfish_altruistic', 'beard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_sac_bd_gb(base_general_configuration, sac_bd_model_configuration,
                                      genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sac_bd_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


### behaviour codominant ###
# beard codominant
@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish'], ['nonbeard', 'nonbeard']],
     ['selfish', 'nonbeard']
     ),
    ([['selfish', 'altruistic'], ['nonbeard', 'beard']],
     ['selfish_altruistic', 'nonbeard_beard']
     ),
    ([['altruistic', 'selfish'], ['beard', 'nonbeard']],
     ['selfish_altruistic', 'nonbeard_beard']
     ),
    ([['altruistic', 'altruistic'], ['beard', 'beard']],
     ['altruistic', 'beard']
     )
])
def test_generate_phenotype_sac_nbc_gb(base_general_configuration, sac_nbc_model_configuration,
                                       genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sac_nbc_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


def test_age_individual_gb(sd_nd_model_configuration):
    general_configuration = [9, 1, 50, 0, 0, 0, 0, None, 8, 0, 1, 3]
    simulation = Simulation(*general_configuration, sd_nd_model_configuration)
    individual = Individual(simulation)
    ages = []
    for i in range(8):
        ages.append(individual.age)
        dead = individual.age_individual()
    assert ages == list(range(8))
    assert dead
