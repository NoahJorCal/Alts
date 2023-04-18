import pytest
from Alts.simulator import Simulation
from Alts.simulator import Individual

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
    return {'module': {'name': 'blind_altruism'},             # Model name
            'behaviour': {'alleles': 'selfish > altruistic',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'}               # Mutation rate
            }


@pytest.fixture(scope='function')
def ad_model_configuration():
    return {'module': {'name': 'blind_altruism'},             # Model name
            'behaviour': {'alleles': 'altruistic > selfish',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'}               # Mutation rate
            }


@pytest.fixture(scope='function')
def c_model_configuration():
    return {'module': {'name': 'blind_altruism'},             # Model name
            'behaviour': {'alleles': 'altruistic = selfish',  # Gene's alleles and inheritance pattern
                          'initial_frequencies': '0.5, 0.5',  # Initial allele's frequencies in order of appearance
                          'mutation_rate': '0'}               # Mutation rate
            }


@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish']],
     ['selfish']
     ),
    ([['selfish', 'altruistic']],
     ['selfish']
     ),
    ([['altruistic', 'selfish']],
     ['selfish']
     ),
    ([['altruistic', 'altruistic']],
     ['altruistic']
     )
])
def test_generate_phenotype_sd_ba(base_general_configuration, sd_model_configuration, genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, sd_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish']],
     ['selfish']
     ),
    ([['selfish', 'altruistic']],
     ['altruistic']
     ),
    ([['altruistic', 'selfish']],
     ['altruistic']
     ),
    ([['altruistic', 'altruistic']],
     ['altruistic']
     )
])
def test_generate_phenotype_ad_ba(base_general_configuration, ad_model_configuration, genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, ad_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


@pytest.mark.parametrize('genotype, expected_phenotype', [
    ([['selfish', 'selfish']],
     ['selfish']
     ),
    ([['selfish', 'altruistic']],
     ['altruistic_selfish']
     ),
    ([['altruistic', 'selfish']],
     ['altruistic_selfish']
     ),
    ([['altruistic', 'altruistic']],
     ['altruistic']
     )
])
def test_generate_phenotype_c_ba(base_general_configuration, c_model_configuration, genotype, expected_phenotype):
    simulation = Simulation(*base_general_configuration, c_model_configuration)
    individual = Individual(simulation)
    # The genotype setter already calls the generate_phenotype() method, but I will call it again in case this changes
    individual.genotype = genotype
    individual.generate_phenotype()
    assert individual.phenotype == expected_phenotype


def test_age_individual_ba(sd_model_configuration):
    general_configuration = [9, 1, 50, 0, 0, 0, 0, None, 8, 0, 1, 3]
    simulation = Simulation(*general_configuration, sd_model_configuration)
    individual = Individual(simulation)
    ages = []
    for i in range(8):
        ages.append(individual.age)
        dead = individual.age_individual()
    assert ages == list(range(8))
    assert dead

