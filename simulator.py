#!/usr/bin/python3
from configparser import ConfigParser
from os import path, makedirs, get_terminal_size
import argparse
import h5py

import time
import warnings

import numpy as np
import bisect
import re
import itertools
import random
from scipy.stats import poisson, truncnorm


# Import general configuration
general_config = ConfigParser()
config_path = path.join(path.dirname(__file__), 'config.ini')
general_config.read(config_path)

# Import model configuration
model_config = ConfigParser()
model_config.read(path.join(path.dirname(__file__), 'models', general_config['simulation']['model'] + '.ini'))
py_module = 'models.' + model_config['module']['name']
model_module = __import__(py_module, fromlist=[''])

# Import variable parameters configuration
var_params_config = ConfigParser()
var_params_path = path.join(path.dirname(__file__), 'variable_config.ini')
var_params_config.read(var_params_path)


# Normal truncated distribution
def truncated_normal(arguments):
    mean = arguments[0]
    low = arguments[1]
    upp = arguments[2]
    sd = arguments[3]
    if sd == 0:
        return mean
    return truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs()


# Print HDF5 file contents
def print_name_type(name, obj):
    print(name, type(obj))


# Increments number of file name if already exists
def uniquify(path_name):
    """
    Uniquify the name of a file by adding a 0 to the end or incrementing the final number.
    :param str path_name: Path to the file which name is going to be uniquified.
    :return: Path for the new file name.
    """
    filename, extension = path.splitext(path_name)
    counter = 0

    while path.exists(path_name):
        path_name = filename + "_" + str(counter) + extension
        counter += 1

    return path_name


# Delete random elements from a list, used for emigration
def delete_random_elems(input_list, n):
    """
    Delete random elements from a list.
    :param list input_list: List from which random elements will be deleted.
    :param int n: Number of elements to delete.
    :return: New list with random elements deleted.
    """
    for i in range(n):
        group_index = random.randrange(len(input_list))
        ind_index = random.randrange(len(input_list[group_index]))
        input_list[group_index].pop(ind_index)
    return input_list


# Length of a nested list, used for calculating the size of a groups structured population
def nested_len(nested_list):
    """
    Calculate length of a nested list.
    :param list[list] nested_list: It is mandatory that the list contains only lists with the objects to be counted.
    :return: The number of elements in the lists of the input list.
    """
    length = 0
    for group in nested_list:
        length += len(group)
    return length


class Chromosome:
    """
    One copy of the two from a locus, contains a list of SNVs and an allele.
    :param str allele: Name of the allele of this Chromosome.
    :ivar np.array snvs: List of floats representing mutations in an infinite sited model genome.
    """
    def __init__(self, allele):
        self.__snvs = np.array([])
        self.__allele = allele

    @property
    def snvs(self):
        return self.__snvs

    @snvs.setter
    def snvs(self, value):
        self.__snvs = value

    @property
    def allele(self):
        return self.__allele

    def mutate(self, locus_size, mutation_rate):
        """
        Adds mutations as floating point numbers to the SNVs list.
        :param int locus_size:
        :param float mutation_rate:
        """
        # Number of mutations calculated based on a Poisson distribution
        mean_mutations = mutation_rate * locus_size
        n_mutations = poisson.rvs(mean_mutations)
        mutations = np.sort([random.random() for _ in range(n_mutations)])
        mutations_indexes = np.searchsorted(self.__snvs, mutations)
        self.__snvs = np.insert(self.__snvs, mutations_indexes, mutations)
        if not list(np.sort(self.__snvs)) == list(self.__snvs):
            print(np.sort(self.__snvs), self.__snvs)
            raise Exception('SNVs list is not ordered')

    # Discontinued: SNVs are now floats instead of objects for performance
    def snvs_to_positions(self):
        return np.array([snv.position for snv in self.__snvs])

    def snvs_to_sequence(self, locus_size):
        """
        Save SNVs lists to haplotype format with 0 representing ancestral variant and 1 representing mutated
        :param int locus_size: Size of the locus used for calculating the position of the floats mutations
        :return: The chromosome's SNVs as a haplotype
        """
        # Position of the mutations
        mutations = [round(snv * locus_size) for snv in self.__snvs]
        sequence = np.zeros(locus_size)
        for mutation in mutations:
            if mutation != locus_size:
                # If the position is already mutated, to avoid recurrent mutations, the next position is mutated
                if sequence[mutation] == 1 and mutation + 1 != locus_size:
                    sequence[mutation + 1] = 1
                else:
                    sequence[mutation] = 1
        return sequence


class Locus:
    """
    Locus object of the genome containing two chromosome, each one with an allele.
    :param list[Chromosome] chromosomes: Lists of the two Chromosome objects of the locus.
    :param str name: Name of the locus.
    :param int locus_size: Size of the locus for calculating mutation and recombination based on rate.
    :param float mutation_rate: Probability of mutation in every position in parts per unit.
    :param float recombination_rate: Probability of recombination in every position in parts per unit.
    """
    def __init__(self, chromosomes, name, locus_size, mutation_rate, recombination_rate):
        self.__chromosomes = chromosomes
        self.__name = name
        self.__locus_size = locus_size
        self.__mutation_rate = mutation_rate
        self.__recombination_rate = recombination_rate

    @property
    def chromosomes(self):
        return self.__chromosomes

    @chromosomes.setter
    def chromosomes(self, value):
        self.__chromosomes = value

    @property
    def name(self):
        return self.__name

    @property
    def locus_size(self):
        return self.__locus_size

    @property
    def mutation_rate(self):
        return self.__mutation_rate

    @property
    def recombination_rate(self):
        return self.__recombination_rate

    def recombine(self, crossovers=None):
        """
        Recombine the chromosomes of the locus based on its recombination rate and size.
        """
        mean_crossovers = self.__recombination_rate * self.__locus_size
        # Crossover can be determined only for testing, in the simulation they are calculated at random
        if not crossovers:
            number_crossovers = poisson.rvs(mean_crossovers)
            crossovers = [random.random() for _ in range(number_crossovers)]
        for crossover_point in crossovers:
            # Index of crossover in each of the chromosomes
            crossover_index_0 = bisect.bisect_left(self.__chromosomes[0].snvs, crossover_point)
            crossover_index_1 = bisect.bisect_left(self.__chromosomes[1].snvs, crossover_point)
            # Chromosomes 0 gets the first part of itself and the second part of chromosome 1
            chromosome_0 = np.concatenate((
                self.__chromosomes[0].snvs[:crossover_index_0],
                self.__chromosomes[1].snvs[crossover_index_1:]))
            # Chromosomes 1 gets the first part of itself and the second part of chromosome 0
            chromosome_1 = np.concatenate((
                self.__chromosomes[1].snvs[:crossover_index_1],
                self.__chromosomes[0].snvs[crossover_index_0:]))
            # Recombined chromosomes are set to the attribute to be recombined again
            self.__chromosomes[0].snvs = chromosome_0
            self.__chromosomes[1].snvs = chromosome_1


class Genome:
    """
    Genome object which has all the genetic information. It contains all the loci of the genetic model.
    :param list[Locus] loci: List of Locus objects.
    """
    def __init__(self, loci):
        self.__loci = loci

    @property
    def loci(self):
        return self.__loci

    @property
    def locus_names(self):
        locus_names = []
        for locus in self.__loci:
            locus_names.append(locus.name)
        return locus_names

    def haplotype(self):
        """
        Haplotypes of both the copies of each locus.
        :return: Tuple of haplotypes.
        """
        full_haplotype = []
        for locus in self.__loci:
            locus_haplotype = []
            for chromosome in locus.chromosomes:
                locus_haplotype.append(chromosome.snvs_to_sequence(locus.locus_size))
            full_haplotype.append(locus_haplotype)
        return tuple(full_haplotype)


class Individual:
    """
    Individual object.
    :param Simulation simulation: Object of the simulation which contains information needed from the Individual object.
    :ivar int age: Age of the individual, each generation increases in 1.
    :ivar int life_expectancy: Age at which the individual will die calculated from a normal distribution.
    :ivar Genome genome: Genome object initialized with empty Locus objects.
    :ivar list genotype: List of lists containing the alleles in each locus.
    :ivar list phenotype: Lists of the phenotype for each locus.
    :ivar float initial_survival_probability: The initial survival probability does not change across the simulation,
    is calculated from a normal distribution and is inherited.
    :ivar float survival_probability: Variable survival probability, can change with altruism.
    :ivar int id: ID of the individual used for calculating relatedness, automatic incremental number.
    :ivar list ancestry: List with the IDs of the ancestors of the individual, each sublist represents a generation
    """
    def __init__(self, simulation):
        self.__simulation = simulation
        self.__age = 0
        self.__life_expectancy = round(np.random.normal(simulation.life_expectancy, simulation.life_expectancy_sd))
        if self.__life_expectancy == 0:
            self.__life_expectancy = 1
        if simulation.simulate_genome:
            loci = []
            for locus in simulation.loci:
                locus_size = int(simulation.model_config_dict[locus]['locus_size'])
                mutation_rate = float(simulation.model_config_dict[locus]['mutation_rate'])
                recombination_rate = float(simulation.model_config_dict[locus]['recombination_rate'])
                loci.append(Locus([], locus, locus_size, mutation_rate, recombination_rate))
            self.__genome = Genome(loci)

        self.__genotype = []
        self.__phenotype = []
        self.__initial_survival_probability = np.random.normal(simulation.survival_probability_mean,
                                                               simulation.survival_probability_sd)
        self.__survival_probability = self.__initial_survival_probability
        self.__id = simulation.newest_ind_id + 1
        # The newest ID of the simulation is updated for the next initialization of an individual
        simulation.newest_ind_id = self.__id
        self.__ancestry = []
        for i in range(simulation.ancestry_generations):
            # List of lists with the ancestors in each generation
            self.__ancestry.append([0] * 2 ** (i + 1))

    @property
    def age(self):
        return self.__age

    @property
    def life_expectancy(self):
        return self.__life_expectancy

    @life_expectancy.setter
    def life_expectancy(self, value):
        self.__life_expectancy = value

    @property
    def genome(self):
        return self.__genome

    @genome.setter
    def genome(self, value):
        self.__genome = value

    @property
    def genotype(self):
        return self.__genotype

    @genotype.setter
    def genotype(self, value):
        """
        The genotype is a lists with a list for each locus with the alleles of that locus. With the genotype,
        the phenotype is calculated and the Genome object is completed
        :param list[list[str]] value: Genotype of the individual as a list of lists containing the alleles of each locus
        """
        self.__genotype = value
        self.generate_phenotype()
        if self.__simulation.simulate_genome:
            self.generate_genome()

    @property
    def phenotype(self):
        return self.__phenotype

    @phenotype.setter
    def phenotype(self, value):
        self.__phenotype = value

    @property
    def initial_survival_probability(self):
        return self.__initial_survival_probability

    @initial_survival_probability.setter
    def initial_survival_probability(self, value):
        self.__initial_survival_probability = value
        self.__survival_probability = value

    @property
    def survival_probability(self):
        return self.__survival_probability

    @survival_probability.setter
    def survival_probability(self, value):
        self.__survival_probability = value

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, value):
        self.__id = value

    @property
    def ancestry(self):
        return self.__ancestry

    @ancestry.setter
    def ancestry(self, value):
        self.__ancestry = value

    def generate_phenotype(self):
        """
        Calculates the phenotype considering the genotype and inheritance pattern of the locus.\n
        Inheritance symbols are:\n
        > : Left dominant\n
        = : codominant
        """
        phenotype = []
        for i in range(len(self.__simulation.loci)):
            locus = self.__simulation.loci[i]
            # If both alleles are the same the phenotype is that allele
            if self.__genotype[i][0] == self.__genotype[i][1]:
                phenotype.append(self.__genotype[i][0])
            else:
                # ?: Non capturing group
                # Matches >, =, any whitespace character, end of line or beginning of line
                inheritance_divider = r'(?:>|=|\s|$|^)'
                characters_between_alleles = re.search(
                    f'{inheritance_divider}{self.__genotype[i][0]}{inheritance_divider}(.*){inheritance_divider}'
                    f'{self.__genotype[i][1]}{inheritance_divider}',
                    self.__simulation.model_config_dict[locus]['alleles'])
                reverse = False
                if not characters_between_alleles:
                    characters_between_alleles = re.search(
                        f'{inheritance_divider}{self.__genotype[i][1]}{inheritance_divider}(.*){inheritance_divider}'
                        f'{self.__genotype[i][0]}{inheritance_divider}',
                        self.__simulation.model_config_dict[locus]['alleles'])
                    reverse = True
                # Dominance-recessiveness
                if '>' in characters_between_alleles.group(1):
                    if reverse:
                        phenotype.append(self.__genotype[i][1])
                    else:
                        phenotype.append(self.__genotype[i][0])
                # Codominance
                else:
                    chosen_phenotype = self.__genotype[i][0] + '_' + self.__genotype[i][1]
                    # [locus][0] is equal to locus's alleles
                    if chosen_phenotype in self.__simulation.loci_properties[locus][0]:
                        phenotype.append(chosen_phenotype)
                    else:
                        phenotype.append(self.__genotype[i][1] + '_' + self.__genotype[i][0])
        self.phenotype = phenotype

    def generate_genome(self):
        """
        Adds Chromosome objects to the genome with the genotype assigned to the individual.
        :return:
        """
        for locus, genotype in zip(self.__genome.loci, self.__genotype):
            chromosomes = []
            for allele in genotype:
                chromosomes.append(Chromosome(allele))
            locus.chromosomes = chromosomes

    def age_individual(self):
        """
        Adds one to the age of the individual.
        :return: True if individual has reached its life expectancy.
        """
        self.__age += 1
        return self.__age == self.__life_expectancy


class Simulation:
    """
    Simulation object that contains all the information that the simulations
    needs to run and all the individuals structured in groups.
    :param int generations: Number of generations in the simulation.
    :param int group_number: Number of groups in the simulation.
    :param int group_size: Initial size of all groups in generation 0.
    :param int group_size_limit: Size limit of a group.
    :param int population_size_limit: Limit of the size of the whole population.
    :param float descendants_per_survivor: Average number of descendants that each survivor will have.
    :param float group_migration: Proportion of migrating individuals in each group.
    :param float emigration: Proportion of emigrants of the whole populations.
    :param float immigration: Proportion of immigrants of the whole populations.
    :param str immigration_phenotype: List of phenotypes per locus of the immigrants.
    :param int life_expectancy: Mean of the life expectancy normal distribution.
    :param float life_expectancy_sd: Standard deviation of the life expectancy normal distribution.
    :param float survival_probability_mean: Mean of the survival probability normal distribution.
    :param float life_expectancy_sd: Standard deviation of the survival probability normal distribution.
    :param int ancestry_generations: Number of generations taken into account in the pedigree for relatedness.
    :param str output_file_name: Name of the HDF5 file where the summary of the simulation will be stored.
    :param dictionary model_config_dict: Dictionary with the configuration of the genetic model.
    :param list[float] altruism_config: List with altruism configuration: exp_factor, cost_benefit_ratio,
    minimum_benefit, maximum_cost, help_higher_sp_probability, help_lower_sp_probability
    :param list[float] selfishness_config: List with selfishness configuration: gained_lost_ratio,
    gained_per_competition, maximum_gained, compete_higher_sp_probability, compete_lower_sp_probability
    :param bool save_data: If True, saves the data of each simulation in the output file.
    :param bool simulate_genome: If True, the genomes with SNVs will be simulated.
    """
    def __init__(
            self,
            generations,
            group_number,
            group_size,
            group_size_limit,
            population_size_limit,
            descendants_per_survivor,
            group_migration,
            emigration,
            immigration,
            immigration_phenotype,
            life_expectancy,
            life_expectancy_sd,
            survival_probability_mean,
            survival_probability_sd,
            ancestry_generations,
            output_file_name,
            model_config_dict,
            altruism_config,
            selfishness_config,
            save_data,
            simulate_genome):
        self.__generations = generations + 1
        self.__group_number = group_number
        self.__group_size = group_size
        self.__groups = []
        self.__calc_avg_survival_prob = []
        self.__exp_avg_survival_prob = []
        self.__groups_indices = np.arange(self.__group_number)
        self.__total_groups = self.__group_number
        self.__group_size_limit = group_size_limit
        self.__population_size_limit = population_size_limit
        self.__descendants_per_survivor = descendants_per_survivor
        self.__group_migration = group_migration
        self.__emigration = emigration
        self.__immigration = immigration
        if immigration_phenotype:
            self.__immigration_phenotype = immigration_phenotype.replace(' ', '').split(',')
        self.__life_expectancy = life_expectancy
        self.__life_expectancy_sd = life_expectancy_sd
        self.__survival_probability_mean = survival_probability_mean
        self.__survival_probability_sd = survival_probability_sd
        self.__generation = 0
        self.__ancestry_generations = ancestry_generations
        self.__model_config_dict = model_config_dict
        self.__loci_properties = {}
        self.__alleles_combinations = []
        self.__altruism_config = altruism_config
        self.__selfishness_config = selfishness_config
        self.__newest_ind_id = 0
        self.__output_file_name = output_file_name
        self.__save_data = save_data
        self.__simulate_genome = simulate_genome
        self.__stop = False
        '''
        Example of model_config_dict:
        {'module': {'name': 'blind_altruism_genomes'},    # Model name
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
        } '''

        for locus in self.__model_config_dict.keys():
            # Ignore module name
            if locus != 'module':
                # Invalid allele's name
                incorrect_character = re.search(r'[^\w\s=>]', self.__model_config_dict[locus]['alleles'])
                if incorrect_character:
                    raise Exception(f'Character "{incorrect_character.group(0)}" in allele not valid')
                # Saving allele names for each locus
                allele_options = re.split(r'\s*>\s*|\s*=\s*', self.__model_config_dict[locus]['alleles'])
                noncodominant_allele_options = allele_options.copy()
                # Saving codominant phenotypes for each locus
                codominant_alleles_raw = re.findall(r'[\w=]*=[\w=]*',
                                                    self.__model_config_dict[locus]['alleles'].replace(' ', ''))
                for element in codominant_alleles_raw:
                    for combination in itertools.combinations(element.split('='), 2):
                        codominant_phenotype = combination[0] + '_' + combination[1]
                        allele_options.append(codominant_phenotype)
                for allele in allele_options:
                    if re.search(r'\s', allele):
                        raise Exception(f'Missing inheritance symbol in "{allele}", allele names can\'t have spaces')
                initial_frequencies = re.split(r'\s*,\s*', self.__model_config_dict[locus]['initial_frequencies'])
                if len(initial_frequencies) != len(noncodominant_allele_options) and \
                        len(initial_frequencies) != len(noncodominant_allele_options) - 1:
                    raise Exception(f'Incorrect number of allele frequencies given in locus "{locus}",'
                                    f'there must be same number as alleles or one less')
                # Initial frequencies are defined by a portion of the range from 0 to 1
                # First frequency range is set manually as it needn't be calculated
                initial_frequencies_ranges = [float(initial_frequencies[0])]
                for i in range(len(initial_frequencies))[1:]:
                    new_range = float(initial_frequencies[i]) + initial_frequencies_ranges[i - 1]
                    if new_range > 1:
                        raise Exception(f'Error in locus "{locus}", the sum of allele frequencies cannot exceed 1')
                    initial_frequencies_ranges.append(new_range)
                if len(initial_frequencies_ranges) != len(noncodominant_allele_options):
                    initial_frequencies_ranges.append(1)
                elif initial_frequencies_ranges[-1] < 1:
                    raise Exception(f'Error in locus "{locus}", the sum of allele frequencies must be 1')
                self.__loci_properties[locus] = (allele_options, initial_frequencies_ranges)

        # Saving allele combinations of the different loci and their order
        all_alleles = []
        for locus in self.loci:
            # [locus][0] are the locus's alleles
            all_alleles.append(self.__loci_properties[locus][0])
        all_alleles_combinations = list(itertools.product(*all_alleles))
        for i in range(len(all_alleles_combinations)):
            allele_combination_name = ''
            for element in all_alleles_combinations[i]:
                allele_combination_name += element + '&'
            if allele_combination_name[:-1] not in self.__alleles_combinations:
                self.__alleles_combinations.append(allele_combination_name[:-1])

        if self.__immigration != 0:
            if type(self.__immigration_phenotype) == str:
                self.__immigration_phenotype = [self.__immigration_phenotype]
            self.__immigrants_genotype = []
            if len(self.loci) != len(self.__immigration_phenotype):
                raise Exception('Incorrect number of phenotypes specified for immigrants')
            for allele, locus_i in zip(self.__immigration_phenotype, range(len(self.__immigration_phenotype))):
                if allele not in self.__loci_properties[self.loci[locus_i]][0]:
                    raise Exception(f'{allele} invalid immigrant phenotype in locus {self.loci[locus_i]}')
                self.__immigrants_genotype.append([allele, allele])

    @property
    def generations(self):
        return self.__generations

    @property
    def current_generation(self):
        return self.__generation

    @property
    def generation(self):
        return self.__generation

    @generation.setter
    def generation(self, value):
        self.__generation = value

    @property
    def ancestry_generations(self):
        return self.__ancestry_generations

    @property
    def groups(self):
        return self.__groups

    @groups.setter
    def groups(self, value):
        self.__groups = value

    @property
    def total_groups(self):
        return self.__total_groups

    @property
    def life_expectancy(self):
        return self.__life_expectancy

    @property
    def life_expectancy_sd(self):
        return self.__life_expectancy_sd

    @property
    def survival_probability_mean(self):
        return self.__survival_probability_mean

    @property
    def survival_probability_sd(self):
        return self.__survival_probability_sd

    @property
    def loci_properties(self):
        return self.__loci_properties

    @property
    def loci(self):
        return list(self.__loci_properties.keys())

    @property
    def dict_allele_options(self):
        alleles = {}
        for locus in self.loci:
            # [locus][0] is equal to locus's alleles
            alleles[locus] = self.__loci_properties[locus][0]
        return alleles

    @property
    def alleles_combinations(self):
        return self.__alleles_combinations

    @property
    def newest_ind_id(self):
        return self.__newest_ind_id

    @newest_ind_id.setter
    def newest_ind_id(self, value):
        self.__newest_ind_id = value

    @property
    def model_config_dict(self):
        return self.__model_config_dict

    @property
    def stop(self):
        return self.__stop

    @property
    def simulate_genome(self):
        return self.__simulate_genome

    def generate_individual_genotype(self):
        """
        Randomize a genotype based on initial locus frequencies.
        """
        allele_pairs = []
        for locus in self.loci:
            alleles = []
            for i in range(2):
                random_picker = random.random()
                # [locus][1] is equal to locus's alleles ranges
                for j in range(len(self.__loci_properties[locus][1])):
                    if random_picker < self.__loci_properties[locus][1][j]:
                        # [locus][0] is equal to locus's alleles
                        alleles.append(self.__loci_properties[locus][0][j])
                        break
            allele_pairs.append(alleles)
        return allele_pairs

    def populate_groups(self):
        """
        Initializes the groups with the specified number of individuals in each one and creates the output HDF5 file.
        """
        if self.__save_data:
            # Creates the output file if the class was not called from a test
            with h5py.File(self.__output_file_name, 'w') as _:
                pass
        for group_i in range(self.__group_number):
            group = []
            for ind_i in range(self.__group_size):
                individual = Individual(self)
                individual.genotype = self.generate_individual_genotype()
                group.append(individual)
            self.__groups.append(group)

    def reset_survival_prob(self):
        """
        Resets the survival probability of each individual to its initial state.
        The effect of altruism in one generation only lasts for that generation
        """
        for group in self.__groups:
            for ind in group:
                ind.survival_probability = ind.initial_survival_probability

    def selection_event(self):
        """
        Survival or death based on the survival probability of the individual. Represents foraging, predators, etc.
        """
        survivors = []
        # np.array for saving which individuals have survived and which haven't
        survived = np.zeros(0)
        survived_list_index = 0
        viable_groups = False
        for group in self.__groups:
            survivors_groups = []
            # The np.array is resized to add the rest of the individuals
            new_survived_size = survived_list_index + len(group)
            survived.resize((new_survived_size,), refcheck=False)
            for ind_i in range(len(group)):
                picker = random.random()
                # The individual will survive if a random generated number
                # between 0 and 1 is below its survival probability
                if picker < group[ind_i].survival_probability:
                    survivors_groups.append(group[ind_i])
                    # 1 indicates that the individual has survived
                    survived[survived_list_index + ind_i] = 1
                else:
                    # 0 indicates that the individual has died
                    survived[survived_list_index + ind_i] = 0
            survivors.append(survivors_groups)
            survived_list_index += len(group)
            # If there is at least one viable group the simulation will continue
            if len(survivors_groups) > 2:
                viable_groups = True
        if not viable_groups:
            self.__stop = 2
        self.__groups = survivors
        # If the method was not called from a test, the survivors data is stored in the output file
        if self.__save_data:
            with h5py.File(self.__output_file_name, 'a') as f:
                generation_group = f[f'/generation_{str(self.current_generation).zfill(len(str(self.__generations)))}']
                generation_group.create_dataset('survivors', data=survived)

    def save_avg_survival_prob(self):
        """
        The average survival probability of the individuals before and after altruism is stored. Groups where the
        difference between expected and gotten is higher will be rewarded because they have helped each other
        for the benefit of the group.
        """
        calc_avg_survival_prob = []
        exp_avg_survival_prob = []
        for group in self.__groups:
            calc_survival_prob = np.zeros(len(group))
            exp_survival_prob = np.zeros(len(group))
            for ind_i in range(len(group)):
                # survival_probability can be modified by altruism
                calc_survival_prob[ind_i] = group[ind_i].survival_probability
                # initial_survival_probability can not be modified by altruism and is inherited
                exp_survival_prob[ind_i] = group[ind_i].initial_survival_probability
            calc_avg_survival_prob.append(np.mean(calc_survival_prob))
            exp_avg_survival_prob.append(np.mean(exp_survival_prob))
        self.__calc_avg_survival_prob = calc_avg_survival_prob
        self.__exp_avg_survival_prob = exp_avg_survival_prob

    def generate_immigrant(self):
        """
        Generates a new individuals with the configured immigrants' genotype.
        :return: Immigrant individual with the configured genotype and no SNVs in genome
        """
        individual = Individual(self)
        # Immigrants are homozygous
        genotype = [[allele, allele] for allele in self.__immigration_phenotype]
        individual.genotype = genotype
        return individual

    def migration(self):
        """
        Emigration and immigration at the population level. The number of emigrants and immigrants are calculated
        based on the population size and the proportion configured. The emigrants come from random groups and
        immigrants enter random groups, the size of the group is not considered.
        """
        population_size = nested_len(self.__groups)

        if self.__emigration != 0:
            # Groups are selected at random
            self.__groups = delete_random_elems(self.__groups, round(population_size * self.__emigration))

        if self.__immigration != 0:
            immigrants = round(population_size * self.__immigration)
            for i in range(immigrants):
                # Groups are selected at random
                group_index = random.choice(range(len(self.__groups)))
                self.__groups[group_index].append(self.generate_immigrant())

    def delete_empty_groups(self):
        """
        Groups with zero individuals are deleted from the groups list. As the group of each individual is saved
        as output data, the original indices of the groups are stored.
        """
        for group_i in range(len(self.__groups) - 1, - 1, - 1):
            if len(self.__groups[group_i]) < 2:
                self.__groups.pop(group_i)
                np.delete(self.__groups_indices, group_i)

    def discard_old_individuals(self):
        """
        Deleted the individuals that will reach their life expectancy in the next generation. These individuals will
        leave offspring but will be deleted before the beginning of the next generation
        :return: Groups without the individuals that will reach their life expectancy in the next generation
        """
        # If all the individuals have a life expectancy of 1,
        # the generations are not overlapped and all of them will die
        if self.__life_expectancy == 1 and self.__life_expectancy_sd == 0:
            new_groups = [[] for _ in range(len(self.__groups))]
        else:
            new_groups = []
            for group in self.__groups:
                new_group = []
                for ind in group:
                    reached_life_expectancy = ind.age_individual()
                    if not reached_life_expectancy:
                        new_group.append(ind)
                new_groups.append(new_group)
        return new_groups

    def calc_new_groups_missing(self, new_groups):
        """
        Calculates the number of individuals that will be born in each group.
        :param list[list[Individual]] new_groups: Groups with the survivors of the generations
        and without the individuals that will
        die of old age in the next generation.
        :return: List of missing individuals per group.
        """
        # Descendants per survivor
        expect_new_inds = [len(group) * self.__descendants_per_survivor for group in self.__groups]
        # Modified survival probabilities / Initial survival probabilities ratio
        calc_exp_sp_ratios = [self.__calc_avg_survival_prob[group_i] / self.__exp_avg_survival_prob[group_i]
                              for group_i in range(len(self.__groups))]
        # Modified number of newborns based on ratio Modified survival probabilities / Initial survival probabilities
        calc_new_inds = [expect_new_inds[group_i] * calc_exp_sp_ratios[group_i]
                         for group_i in range(len(self.__groups))]
        # Scaled number of newborns to be the same as the expected calculated
        scaled_new_inds = [round(group_size * sum(expect_new_inds) / sum(calc_new_inds))
                           for group_size in calc_new_inds]
        missing_population = self.__population_size_limit - nested_len(new_groups)
        # If adding the newborns is going to make the population exceed the size limit, all the sizes are scaled down.
        if sum(scaled_new_inds) > missing_population:
            scaled_new_inds = [round(group_size * missing_population / sum(scaled_new_inds))
                               for group_size in scaled_new_inds]
        return scaled_new_inds

    def generate_offspring_genome(self, reproducers, descendant):
        """
        Generates the genome of the offspring based on the reproducers genomes.
        For each locus, there is a 50/50 chance of inherit one of the two chromosomes per parent.
        Before the chromosome is inherited, it is mutated and recombined.
        :param list[Individual] reproducers: List with the sire and the dam of the descendant, the order is irrelevant.
        :param Individual descendant: Offspring that the generated genome will be assigned to.
        """
        # Genotype in name of alleles format used in the altruistic step
        genotype = []
        # Actual Locus objects for generating the Genome object
        loci = []
        for locus_index in range(len(self.loci)):
            genotype.append([])
            chromosomes = []
            # For each locus, each reproducer contribute with one allele (Chromosome object)
            for reproducer in reproducers:
                # The chromosome that will be inherited is selected at random between the two
                chromosome_index = random.randint(0, 1)
                if self.__simulate_genome:
                    # The locus of the parent is cloned
                    locus = reproducer.genome.loci.copy()[locus_index]
                    if locus.recombination_rate != 0:
                        locus.recombine()
                    chromosome = locus.chromosomes[chromosome_index]
                    if locus.mutation_rate != 0:
                        if locus.mutation_rate * locus.locus_size == 0:
                            warnings.warn(f'Mutation rate times locus size is equal to 0 in locus {locus.name}, '
                                          f'there will be no mutations')
                        chromosome.mutate(locus.locus_size, locus.mutation_rate)
                    chromosomes.append(chromosome)
                genotype[locus_index].append(reproducer.genotype[locus_index][chromosome_index])

            if self.__simulate_genome:
                # The chromosomes list is ready to initialize the Locus object
                name = self.loci[locus_index]
                # The information of the locus is extracted from one of the parents, it is the same for every individual
                locus_size = reproducers[0].genome.loci[locus_index].locus_size
                mutation_rate = reproducers[0].genome.loci[locus_index].mutation_rate
                recombination_rate = reproducers[0].genome.loci[locus_index].recombination_rate
                # The Locus object is initialized and stored in the list that will be passed to the Genome object
                loci.append(Locus(chromosomes, name, locus_size, mutation_rate, recombination_rate))
        descendant.genotype = genotype
        if self.__simulate_genome:
            descendant.genome = Genome(loci)

    def generate_offspring_ancestry(self, reproducers, descendant):
        """
        Generates the ancestry list of the offspring based on the reproducers ancestries. The first sublist will be
        the IDs of the reproducers and the rest will be extracted from their own ancestry.
        :param list[Individual] reproducers: List with the sire and the dam of the descendant, the order is irrelevant.
        :param Individual descendant: Offspring that the generated ancestry will be assigned to.
        """
        for generation in range(self.__ancestry_generations):
            if generation == 0:
                # First generation of the ancestry are the parents
                descendant.ancestry[generation] = [reproducers[0].id, reproducers[1].id]
            else:
                # The rest are extracted from the parents, the first half is one parent and the second the other
                descendant.ancestry[generation] = reproducers[0].ancestry[generation - 1] + \
                                                  reproducers[1].ancestry[generation - 1]

    def divide_big_groups(self, groups):
        """
        Groups that have reach the group size limit are divided in two at random. Families are not kept together.
        :param list[Individuals] groups: Groups of individuals.
        """
        divided_groups = []
        # If any group is divided the indices for saving the summary of the simulation need to be updated
        new_groups_indices = []
        for group_i in range(len(groups)):
            if len(groups[group_i]) > self.__group_size_limit:
                # The index of the group to be divided is deleted, the two new groups will have new indices
                np.delete(self.__groups_indices, group_i)
                # Individuals are added to the group by age, the group is shuffled to avoid grouping
                # older individuals with older individuals and younger with younger
                random.shuffle(groups[group_i])
                # First half of the group is added as another group
                divided_groups.append(groups[group_i][:round(len(groups[group_i]) / 2)].copy())
                new_groups_indices.append(self.__groups_indices[group_i])
                # The original group index will be the second half
                divided_groups.append(groups[group_i][round(len(groups[group_i]) / 2):].copy())
                new_groups_indices.append(max(max(self.__groups_indices), max(new_groups_indices)) + 1)
                self.__total_groups += 1
            else:
                divided_groups.append(groups[group_i])
                new_groups_indices.append(self.__groups_indices[group_i])
        self.__groups_indices = np.array(new_groups_indices)
        self.__groups = divided_groups

    def group_exchange(self, survivors, newborns):
        """
        Migration at group level. The individuals born in the current generation, if any, will be the migrants.
        :param list[list[Individual]] survivors: List of groups of individuals that will pass to the next generation.
        :param list[list[Individual]] newborns: List of groups of individuals born in the current generation.
        :return: Lists of joined survivors and newborns with group exchange.
        """
        # If there is more than one group and there is group migration, there will be exchange between groups
        if len(survivors) > 1 and self.__group_migration > 0:
            for from_group_index in range(len(survivors)):
                group_size = len(survivors[from_group_index]) + len(newborns[from_group_index])
                # Number of emigrants of the group
                emigrants = round(group_size * self.__group_migration)
                target_indexes = list(range(len(survivors)))
                # Possible target groups
                target_indexes.remove(from_group_index)
                # If there have been newborn in the group, the migrants will be selected from them
                if len(newborns[from_group_index]) != 0:
                    # If there are enough newborns
                    if len(newborns[from_group_index]) > emigrants:
                        # The emigrants are selected at random from the group and sent to a random group
                        for i in range(emigrants):
                            emigrant = newborns[from_group_index].pop(random.randrange(len(newborns[from_group_index])))
                            target_group_index = random.choice(target_indexes)
                            survivors[target_group_index].append(emigrant)
                    # If there are not enough newborns, the rest will be selected from the individuals
                    # that were already in the population
                    else:
                        emigrants_left = emigrants - len(newborns[from_group_index])
                        for i in range(len(newborns[from_group_index])):
                            emigrant = newborns[from_group_index].pop(random.randrange(len(newborns[from_group_index])))
                            target_group_index = random.choice(target_indexes)
                            survivors[target_group_index].append(emigrant)
                        for i in range(emigrants_left):
                            emigrant = survivors[from_group_index].pop(
                                random.randrange(len(survivors[from_group_index])))
                            target_group_index = random.choice(target_indexes)
                            survivors[target_group_index].append(emigrant)
                # If there are no newborns, but there are individuals from the previous generation (no individuals
                # were born) the migrants will be selected from them.
                elif len(survivors[from_group_index]) != 0:
                    for i in range(emigrants):
                        emigrant = survivors[from_group_index].pop(random.randrange(len(survivors[from_group_index])))
                        target_group_index = random.choice(target_indexes)
                        survivors[target_group_index].append(emigrant)
            # The survivors groups are joined with the newborns after the migration
            for group_i in range(len(survivors)):
                survivors[group_i].extend(newborns[group_i])
            return survivors
        # There will be no migration
        else:
            for group_i in range(len(survivors)):
                survivors[group_i].extend(newborns[group_i])
            return survivors

    def reproduce(self):
        """
        The reproduction occurs inside groups and is independent of the other groups. Parents are selected at random
        from the individuals of the group.
        """
        self.delete_empty_groups()
        new_groups = self.discard_old_individuals()
        offspring = [[] for _ in range(len(self.__groups))]
        new_group_missing = self.calc_new_groups_missing(new_groups)
        for group_index in range(len(self.__groups)):
            group_size = new_group_missing[group_index]
            while len(offspring[group_index]) < group_size:
                reproducers = random.sample(self.__groups[group_index], 2)
                # Survival probability is a mean of the initial survival probabilities of the parents
                survival_probability = (reproducers[0].initial_survival_probability +
                                        reproducers[1].initial_survival_probability) / 2
                new_individual = Individual(self)
                new_individual.initial_survival_probability = survival_probability
                self.generate_offspring_genome(reproducers, new_individual)
                self.generate_offspring_ancestry(reproducers, new_individual)
                offspring[group_index].append(new_individual)
        new_groups = self.group_exchange(new_groups, offspring)
        self.divide_big_groups(new_groups)

    def save_generation_data(self):
        """
        From each individual of the generation it is saved the group it belongs to, the phenotype and the genotype. The
        information of survival is saved in the selection_event method.
        :return:
        """
        with h5py.File(self.__output_file_name, 'a') as f:
            # The HDF5 file will have a group for each generation
            generation_group = f.create_group(
                f'generation_{str(self.current_generation).zfill(len(str(self.__generations)))}')
            # Alleles of each locus to save the genotype as integers
            loci_alleles = []
            for locus_i in range(len(self.loci)):
                loci_alleles.append(self.__loci_properties[self.loci[locus_i]][0])
            groups = np.zeros(0)
            phenotypes = np.zeros(0)
            genotypes = []
            for _ in self.loci:
                genotypes.append(np.zeros(0))
            # For resizing the np arrays
            individuals_list_index = 0
            alleles_list_index = 0
            for group in self.__groups:
                new_phenotypes_size = individuals_list_index + len(group)
                phenotypes.resize((new_phenotypes_size,), refcheck=False)
                new_groups_size = individuals_list_index + len(group)
                groups.resize((new_groups_size,), refcheck=False)
                new_genotypes_size = alleles_list_index + (len(group) * 2)
                for i in range(len(self.loci)):
                    genotypes[i].resize((new_genotypes_size,), refcheck=False)
                genotypes_index = 0
                for ind_i in range(len(group)):
                    # Index of the individual's phenotype in the alleles_combinations list
                    phenotypes[individuals_list_index + ind_i] = self.__alleles_combinations.index(
                        '&'.join(group[ind_i].phenotype))
                    # Index of individual's group
                    groups[individuals_list_index + ind_i] = self.__groups_indices[self.__groups.index(group)]
                    for allele_i in range(2):
                        for locus_i in range(len(self.loci)):
                            # Index of the individual's aleles in the loci_alleles list
                            genotypes[locus_i][alleles_list_index + genotypes_index] = loci_alleles[locus_i].index(
                                group[ind_i].genotype[locus_i][allele_i])
                        genotypes_index += 1
                individuals_list_index += len(group)
                alleles_list_index += len(group) * 2
            else:
                # If the simulation does not stop, the data is saved in the HDF5 file
                generation_group.create_dataset('group', data=groups)
                generation_group.create_dataset('phenotype', data=phenotypes)
                generation_group['phenotype'].attrs['phenotype_names'] = self.__alleles_combinations
                for locus_i in range(len(self.loci)):
                    generation_group.create_dataset(f'locus_{self.loci[locus_i]}', data=genotypes[locus_i])

    def check_stop_simulation(self):
        altruists = 0
        selfish = 0
        for group in self.__groups:
            for ind in group:
                if 'altruistic' in ind.phenotype[0]:
                    altruists += 1
                elif 'selfish' in ind.phenotype[0]:
                    selfish += 1
        # If there are no altruists or no selfish the simulation will stop
        # print(altruists, selfish)
        if not altruists or not selfish:
            self.__stop = 1

    def pass_generation(self):
        """
        Simulation one generation:\n
        1. Resets the survival probabilities of the individuals\n
        2. Altruists help other individuals\n
        3. The average initial and modified survival probabilities is calculated\n
        4. Individuals die or survive based on their survival probability\n
        5. Migration at the population level\n
        6. Reproduction, which include migration between groups\n
        7. The data of the generation is saved
        """
        d = {'selfish': [0, 0],
             'selfish_altruistic': [0, 0],
             'altruistic': [0, 0]}
        print(self.current_generation)
        for g in self.__groups:
            for i in g:
                d[str(i.phenotype[0])][0] += 1
        suma = 0
        for i in d.values():
            suma += i[0]
        for key, value in d.items():
            d[key][1] = round(value[0] / suma, 4)
        print('summary', d)
        self.reset_survival_prob()
        model_module.selection(self.groups, self.__altruism_config, self.__selfishness_config)
        self.save_avg_survival_prob()
        self.selection_event()
        if not self.__stop:
            self.migration()
            self.reproduce()
            self.__generation += 1
            if self.__save_data:
                self.save_generation_data()
            self.check_stop_simulation()

    def save_haplotypes(self):
        haplotypes_set = [set() for _ in self.loci]
        for group in self.__groups:
            for individual in group:
                loci = individual.genome.loci
                for locus_index in range(len(self.loci)):
                    for chromosome in loci[locus_index].chromosomes:
                        chrom_haplotype = chromosome.snvs_to_sequence(loci[locus_index].locus_size)
                        haplotypes_set[locus_index].add(tuple(chrom_haplotype))
        haplotypes = []
        for haplotype_set in haplotypes_set:
            haplotypes.append(np.array(list(haplotype_set)))

        if self.__save_data:
            with h5py.File(self.__output_file_name, 'a') as f:
                f.create_group('haplotypes')
                for locus_i in range(len(self.loci)):
                    f.create_dataset(f'haplotypes/haplotypes_{self.loci[locus_i]}', data=haplotypes[locus_i])

    def save_snvs(self):
        """
        Saves the SNVs lists of the individuals in the current generation as a CSV file for each locus,
        each row being each individual.
        """
        snvs_lists = [[] for _ in range(len(self.loci))]
        for group in self.__groups:
            for individual in group:
                loci = individual.genome.loci
                for locus_index in range(len(self.loci)):
                    # For each chromosome of each locus of each individual of each group, the SNVs lists are stored
                    for chromosome in loci[locus_index].chromosomes:
                        snvs_lists[locus_index].append(chromosome.snvs)
        # The directory is created if it does not exist
        snvs_dir = path.join(path.dirname(__file__), 'snvs')
        makedirs(snvs_dir, exist_ok=True)
        # The haplotypes of each locus will be stored in different files
        for locus_i in range(len(self.loci)):
            file_name = f'{self.loci[locus_i]}_snvs.csv'
            result = path.join(path.dirname(__file__), 'snvs', file_name)
            result = uniquify(result)
            with open(result, 'a') as f:
                for haplotype in snvs_lists[locus_i]:
                    f.write(",".join(str(snv) for snv in haplotype))
                    f.write('\n')

    def altruist_perc(self):
        altruists = 0
        for group in self.__groups:
            for ind in group:
                if 'altruistic' in ind.phenotype[0]:
                    altruists += 1
        population_size = nested_len(self.__groups)
        altruist_percentage = altruists / population_size
        return altruist_percentage


def simulator_main(save_data, output_dir, output_file, simulate_genome, sim_seed=None, quiet=False):
    """
    Main function of the simulation that gets the configuration from the config.ini file, initializes the Simulation
    object, runs the number of generations configured and adds the metadata of the HDF5 file.
    :param bool save_data: If True, saves the data of each simulation in the output file.
    :param str output_dir: Output directory name where the HDF5 files will be stored.
    :param str output_file: Output HDF5 file name where the data will be stored.
    :param bool simulate_genome: If True, the genomes with SNVs will be simulated.
    :param int sim_seed: Seed used to fix simulation's result, defaults to None.
    :param bool quiet: Flag to avoid printing the program's feedback, defaults to False.
    :return: A tuple with whether the simulation ended because of lack of altruists or individuals and the name of the
    output file, as it could have changed because of uniquify.
    """
    # Generate random seed for parameters initialization
    random.seed()
    np.random.seed()

    if save_data:
        # The output directory is created if it did not exist
        makedirs(output_dir, exist_ok=True)
        # The extension is added to the HDF5 file if it was not present and the name is uniquified
        if '.h5' not in output_file and '.hdf5' not in output_file and \
           '.h5p' not in output_file and '.he5' not in output_file and \
           '.h5m' not in output_file and '.h5z' not in output_file:
            output_file += '.h5'
    output_file_name = path.join(path.dirname(__file__), output_dir, output_file)
    output_file_name = uniquify(output_file_name)

    # Initial parameters from the configuration file
    generations = int(general_config['simulation']['generations'])
    group_number = int(general_config['population']['group_number'])
    descendants_per_survivor = float(general_config['population']['descendants_per_survivor'])
    emigration = float(general_config['population']['emigration'])
    immigration = float(general_config['population']['immigration'])
    immigration_phenotype = general_config['population']['immigration_phenotype']
    life_expectancy = int(general_config['population']['life_expectancy'])
    life_expectancy_sd = float(general_config['population']['life_expectancy_sd'])
    ancestry_generations = int(general_config['population']['ancestry_generations'])
    # Randomly generating variable parameters
    group_size = int(truncated_normal([float(elem) for elem in var_params_config['parameters']['group_size']
                                      .replace(' ', '').split(',')]))
    group_size_limit = int(truncated_normal([float(elem) for elem in var_params_config['parameters']['group_size_limit']
                                            .replace(' ', '').split(',')]))
    population_size_limit = int(truncated_normal([float(elem) for elem in
                                                  var_params_config['parameters']['population_size_limit']
                                                 .replace(' ', '').split(',')]))
    group_migration = truncated_normal([float(elem) for elem in var_params_config['parameters']['group_migration']
                                       .replace(' ', '').split(',')])
    survival_probability_mean = truncated_normal([float(elem)for elem in
                                                  var_params_config['parameters']['survival_probability_mean']
                                                 .replace(' ', '').split(',')])
    survival_probability_sd = truncated_normal([float(elem) for elem in
                                                var_params_config['parameters']['survival_probability_sd'].
                                               replace(' ', '').split(',')])
    altruism_exp_factor = truncated_normal([float(elem) for elem in
                                            var_params_config['parameters']['benefit_relatedness_exp_factor']
                                           .replace(' ', '').split(',')])
    cost_benefit_ratio = truncated_normal([float(elem) for elem in var_params_config['parameters']['cost_benefit_ratio']
                                          .replace(' ', '').split(',')])
    minimum_benefit = truncated_normal([float(elem) for elem in var_params_config['parameters']['minimum_benefit']
                                       .replace(' ', '').split(',')])
    maximum_cost = truncated_normal([float(elem) for elem in var_params_config['parameters']['maximum_cost']
                                    .replace(' ', '').split(',')])
    help_higher_sp_probability = truncated_normal([float(elem) for elem in
                                                   var_params_config['parameters']['help_higher_sp_probability']
                                                  .replace(' ', '').split(',')])
    help_lower_sp_probability = truncated_normal([float(elem) for elem in
                                                  var_params_config['parameters']['help_lower_sp_probability']
                                                 .replace(' ', '').split(',')])
    gained_lost_ratio = truncated_normal([float(elem) for elem in var_params_config['parameters']['gained_lost_ratio']
                                         .replace(' ', '').split(',')])
    gained_per_competition = truncated_normal([float(elem) for elem in
                                               var_params_config['parameters']['gained_per_competition']
                                              .replace(' ', '').split(',')])
    maximum_gained = truncated_normal([float(elem) for elem in var_params_config['parameters']['maximum_gained']
                                      .replace(' ', '').split(',')])
    compete_higher_sp_probability = truncated_normal([float(elem) for elem in
                                                      var_params_config['parameters']['compete_higher_sp_probability']
                                                     .replace(' ', '').split(',')])
    compete_lower_sp_probability = truncated_normal([float(elem) for elem in
                                                     var_params_config['parameters']['compete_lower_sp_probability']
                                                    .replace(' ', '').split(',')])

    altruism_configuration = [altruism_exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,
                              help_higher_sp_probability, help_lower_sp_probability]
    selfishness_configuration = [gained_lost_ratio, gained_per_competition, maximum_gained,
                                 compete_higher_sp_probability, compete_lower_sp_probability]
    # The model configuration is passed to the Simulation object as a dictionary
    model_config_dict = {s: dict(model_config.items(s)) for s in model_config.sections()}
    altruism_initial_freq = truncated_normal([float(elem) for elem in
                                              var_params_config['parameters']['altruism_initial_freq']
                                             .replace(' ', '').split(',')])
    model_config_dict['behaviour']['initial_frequencies'] = str(1 - altruism_initial_freq)

    # Set seed for getting the same results everytime
    if sim_seed:
        random.seed(sim_seed)
        np.random.seed(sim_seed)

    parameters = [group_size, group_size_limit, population_size_limit, group_migration, survival_probability_mean,
                  survival_probability_sd, *altruism_configuration, *selfishness_configuration, altruism_initial_freq]
    # Initial parameters are saved in the output file
    with open('output.csv', 'a') as out_file:
        out_file.write(str(parameters)[1:-1].replace(' ', '') + ',')

    # The simulation is initialized and set up by populating the groups and
    # saving the state of the simulation at the beginning of generation 0
    simulation = Simulation(generations,
                            group_number,
                            group_size,
                            group_size_limit,
                            population_size_limit,
                            descendants_per_survivor,
                            group_migration,
                            emigration,
                            immigration,
                            immigration_phenotype,
                            life_expectancy,
                            life_expectancy_sd,
                            survival_probability_mean,
                            survival_probability_sd,
                            ancestry_generations,
                            output_file_name,
                            model_config_dict,
                            altruism_configuration,
                            selfishness_configuration,
                            save_data,
                            simulate_genome)
    simulation.populate_groups()
    if save_data:
        simulation.save_generation_data()

        # The duration of each generation is stored and saved in the output file
        generations_duration = np.zeros(generations)
        generation_start = time.perf_counter()

    # Generations loop
    for i in range(generations):
        simulation.pass_generation()
        # If there are no more altruists or individuals the simulation stops
        if simulation.stop:
            break
        if save_data:
            generation_end = time.perf_counter()
            generations_duration[i] = generation_end - generation_start
            generation_start = time.perf_counter()
        if not quiet:
            # Progress bar
            bar_msg = 'Simulation progress: '
            cols = get_terminal_size().columns - len(bar_msg)
            progress = cols*i/generations
            bar_char = '█'
            bar_end_chars = ' ▏▎▍▌▋▊▉'
            print('\r' + bar_msg + bar_char*int(progress) + bar_end_chars[int((progress - int(progress))*8)] +
                  ' ' * (cols - int(progress) - 1), end='')

    # At the end of the simulation, the SNVs are saved in the corresponding files
    if simulation.simulate_genome:
        simulation.save_snvs()

    if not simulation.stop or simulation.stop == 1:
        result = [simulation.altruist_perc(), simulation.current_generation]
        with open('output.csv', 'a') as out_file:
            out_file.write(str(result)[1:-1].replace(' ', '') + '\n')

        if save_data:
            # The metadata for plotting the data is stored
            phenotypes_names = []
            alleles_names = []
            for locus, alleles in simulation.loci_properties.items():
                phenotypes_names.append(alleles[0])
                alleles_names.append(alleles[0][:len(alleles[1])])

            with h5py.File(output_file_name, 'a') as f:
                f.attrs['generations'] = generations + 1
                f.attrs['n_loci'] = len(simulation.loci)
                f.attrs['groups'] = simulation.total_groups
                f.attrs['loci'] = simulation.loci
                f.attrs['phenotype_names'] = simulation.alleles_combinations
                f.attrs['phenotypes_names'] = str(phenotypes_names)
                f.attrs['alleles_names'] = str(alleles_names)
                f.create_dataset('duration', data=generations_duration)

        # Return whether the simulation was aborted because of lack of individuals
        return False, output_file_name
    else:
        return True, output_file_name


if __name__ == '__main__':
    # Argument parser
    parser = argparse.ArgumentParser(prog='Alts simulator',
                                     description='Main simulation of Alts program')
    parser.add_argument('-ns', '--no-save', dest='save_data', action='store_false',
                        help='The data will not be collected')
    parser.add_argument('-d', '--directory', default='output', help='Output directory where individual'
                                                                    'simulation results will be stored')
    parser.add_argument('-o', '--output', default='simulation.h5', help='Output file where the'
                                                                        'simulation results will be stored')
    parser.add_argument('-g', '--genome', dest='simulate_genome', action='store_true',
                        help='The genomes will be simulated')
    parser.add_argument('-s', '--seed', required=False, type=int, help='Set a seed for the simulation')

    args = parser.parse_args()
    global_output_file_name = args.output
    if '.h5' not in global_output_file_name and '.hdf5' not in global_output_file_name and \
       '.h5p' not in global_output_file_name and '.he5' not in global_output_file_name and \
       '.h5m' not in global_output_file_name and '.h5z' not in global_output_file_name:
        global_output_file_name += '.h5'
    seed = args.seed

    simulator_main(args.save_data, args.directory, args.output, args.simulate_genome, args.seed, args.quiet)
