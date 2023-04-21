#!/usr/bin/python3
from os import path, get_terminal_size
import sys
import re
import random
from typing import List, Any

import numpy
import itertools
from configparser import ConfigParser

# Add the models directory to the Python path
sys.path.append(path.join(path.dirname(__file__), 'models'))

# Import general configuration
general_config = ConfigParser()
config_path = path.join(path.dirname(__file__), 'config.ini')
general_config.read(config_path)

# Import model configuration
model_config = ConfigParser()
model_config.read(path.join(path.dirname(__file__), 'models', general_config['simulation']['model'] + '.ini'))
py_module = model_config['module']['name']
model_module = __import__(py_module, fromlist=[''])


def delete_random_elems(input_list, n):
    for i in range(n):
        group_index = random.randrange(len(input_list))
        ind_index = random.randrange(len(input_list[group_index]))
        input_list[group_index].pop(ind_index)
    return input_list


def nested_len(nested_list):
    length = 0
    for group in nested_list:
        length += len(group)
    return length


class Simulation:
    def __init__(
            self,
            generations,
            group_number,
            group_size,
            group_size_sd,
            group_migration,
            emigration,
            immigration,
            immigration_phenotype,
            life_expectancy,
            life_expectancy_sd,
            survival_probability,
            ancestry_generations,
            model_config_dict):
        self.__generations = generations + 1
        self.__group_number = group_number
        self.__group_size = group_size
        self.__group_size_sd = group_size_sd
        self.__groups = []
        self.__group_migration = group_migration
        self.__emigration = emigration
        self.__immigration = immigration
        self.__immigration_phenotype = immigration_phenotype
        self.__life_expectancy = life_expectancy
        self.__life_expectancy_sd = life_expectancy_sd
        self.__survival_probability = survival_probability
        self.__generation = 0
        self.__ancestry_generations = ancestry_generations
        self.__model_config_dict = model_config_dict
        self.__genes_properties = {}
        # self.__survival_rate_mean = 0
        self.__survivors_per_generation = []
        ''' First list of the summary is survivors per generation and second the number of individuals
        at the beginning of a generation '''
        self.__simulation_summary = [[], []]
        self.__alleles_combinations_indexes = []
        self.__newest_ind_id = 0
        for gene in self.__model_config_dict.keys():
            if gene != 'module':
                incorrect_character = re.search(r'[^\w\s=>]', self.__model_config_dict[gene]['alleles'])
                if incorrect_character:
                    raise Exception(f'Character "{incorrect_character.group(0)}" in allele not valid')
                # Saving allele names for each gene
                allele_options = re.split(r'\s*>\s*|\s*=\s*', self.__model_config_dict[gene]['alleles'])
                noncodominant_allele_options = allele_options.copy()
                # Saving codominant phenotypes for each gene
                codominant_alleles_raw = re.findall(r'[\w=]*=[\w=]*', self.__model_config_dict[gene]['alleles'].replace(' ', ''))
                for element in codominant_alleles_raw:
                    for combination in itertools.combinations(element.split('='), 2):
                        codominant_phenotype = combination[0] + '_' + combination[1]
                        allele_options.append(codominant_phenotype)
                for allele in allele_options:
                    if re.search(r'\s', allele):
                        raise Exception(f'Missing inheritance symbol in "{allele}", allele names can\'t have spaces')
                initial_frequencies = re.split(r'\s*,\s*', self.__model_config_dict[gene]['initial_frequencies'])
                if len(initial_frequencies) != len(noncodominant_allele_options) and \
                        len(initial_frequencies) != len(noncodominant_allele_options) - 1:
                    raise Exception(f'Incorrect number of allele frequencies given in gene "{gene}",'
                                    f'there must be same number as alleles or one less')
                # Initial frequencies are defined by a portion of the range from 0 to 1
                # First frequency range is set manually as it needn't be calculated
                initial_frequencies_ranges = [float(initial_frequencies[0])]
                for i in range(len(initial_frequencies))[1:]:
                    new_range = float(initial_frequencies[i]) + initial_frequencies_ranges[i - 1]
                    if new_range > 1:
                        raise Exception(f'Error in gene "{gene}", the sum of allele frequencies cannot exceed 1')
                    initial_frequencies_ranges.append(new_range)
                if len(initial_frequencies_ranges) != len(noncodominant_allele_options):
                    initial_frequencies_ranges.append(1)
                elif initial_frequencies_ranges[-1] < 1:
                    raise Exception(f'Error in gene "{gene}", the sum of allele frequencies must be 1')
                self.__genes_properties[gene] = (allele_options, initial_frequencies_ranges)

        # Saving allele combinations of the different genes and their order
        alleles_combinations_indexes = {}
        all_alleles = []
        for gene in self.genes:
            # [gene][0] are the gene's alleles
            all_alleles.append(self.__genes_properties[gene][0])
        all_alleles_combinations = list(itertools.product(*all_alleles))
        for i in range(len(all_alleles_combinations)):
            allele_combination_name = ''
            for element in all_alleles_combinations[i]:
                allele_combination_name += element + '&'
            alleles_combinations_indexes[allele_combination_name[:-1]] = i
        self.__alleles_combinations_indexes = alleles_combinations_indexes

        if self.__immigration != 0:
            if type(self.__immigration_phenotype) == str:
                self.__immigration_phenotype = [self.__immigration_phenotype]
            self.__immigrants_genotype = []
            if len(self.genes) != len(self.__immigration_phenotype):
                raise Exception('Incorrect number of phenotypes specified for immigrants')
            for allele in self.__immigration_phenotype:
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
    def life_expectancy(self):
        return self.__life_expectancy

    @property
    def life_expectancy_sd(self):
        return self.__life_expectancy_sd

    @property
    def survival_probability(self):
        return self.__survival_probability

    @property
    def genes_properties(self):
        return self.__genes_properties

    @property
    def genes(self):
        return list(self.__genes_properties.keys())

    @property
    def simulation_summary(self):
        return self.__simulation_summary

    @simulation_summary.setter
    def simulation_summary(self, value):
        self.__simulation_summary = value

    @property
    def dict_allele_options(self):
        alleles = {}
        for gene in self.genes:
            # [gene][0] is equal to gene's alleles
            alleles[gene] = self.__genes_properties[gene][0]
        return alleles

    @property
    def alleles_combinations_indexes(self):
        return self.__alleles_combinations_indexes

    @property
    def newest_ind_id(self):
        return self.__newest_ind_id

    @newest_ind_id.setter
    def newest_ind_id(self, value):
        self.__newest_ind_id = value

    @property
    def model_config_dict(self):
        return self.__model_config_dict

    # Randomize individual's genotype based on initial gene frequencies
    def generate_individual_genotype(self):
        allele_pairs = []
        for gene in self.genes:
            alleles = []
            for i in range(2):
                random_picker = random.random()
                # [gene][1] is equal to gene's alleles ranges
                for j in range(len(self.__genes_properties[gene][1])):
                    if random_picker < self.__genes_properties[gene][1][j]:
                        # [gene][0] is equal to gene's alleles
                        alleles.append(self.__genes_properties[gene][0][j])
                        break
            allele_pairs.append(alleles)
        return allele_pairs

    def populate_groups(self):
        for group_i in range(self.__group_number):
            group = []
            group_size_normal = round(numpy.random.normal(self.__group_size, self.__group_size_sd))
            for ind_i in range(group_size_normal):
                individual = Individual(self)
                individual.genotype = self.generate_individual_genotype()
                self.__newest_ind_id = individual.id
                group.append(individual)
            self.__groups.append(group)

    # Survival or death based on the survival probability of the individual. Represents foraging, predators, etc.
    def selection_event(self):
        self.__survivors_per_generation = []
        survival_probabilities = []
        for group in self.__groups:
            survivors_groups = []
            for individual in group:
                survival_probabilities.append(individual.survival_probability)
                picker = random.random()
                if picker < individual.survival_probability:
                    survivors_groups.append(individual)
            self.__survivors_per_generation.append(survivors_groups)
        self.__groups = self.__survivors_per_generation
        # self.__survival_rate_mean = sum(survival_probabilities)/len(survival_probabilities)

    # Generates a new individuals with the given immigrant's genotype
    def generate_immigrant(self):
        individual = Individual(self)
        self.__newest_ind_id = individual.id
        individual.genotype = self.__immigrants_genotype
        return individual

    def migration(self):
        full_population_size = self.__group_size * self.__group_number
        real_population_size = nested_len(self.__groups)

        if self.__emigration != 0:
            self.__groups = delete_random_elems(self.__groups, round(real_population_size * self.__emigration))

        missing_population = full_population_size - real_population_size
        if self.__immigration != 0:
            immigrants = min(round(full_population_size * self.__immigration), missing_population)
            for i in range(immigrants):
                group_index = random.choice(range(len(self.__groups)))
                self.__groups[group_index].append(self.generate_immigrant())

    def mutate_genotype(self, reproducer_allele, gene_index):
        mutation_rate = float(self.__model_config_dict[self.genes[gene_index]]['mutation_rate'])
        if random.random() < mutation_rate:
            allele_options = self.__genes_properties[self.genes[gene_index]][0].copy()
            allele_options.remove(reproducer_allele)
            return random.choice(allele_options)
        else:
            return reproducer_allele

    def generate_offspring_genotype(self, reproducers, descendant):
        genotype = []
        # For each gene there's a probability of mutation
        for gene_index in range(len(reproducers[0].genotype)):
            # Chose one allele of each gene of each reproducer
            gene_alleles = []
            for reproducer_index in range(2):
                reproducer_allele = random.choice(reproducers[reproducer_index].genotype[gene_index])
                if float(self.__model_config_dict[self.genes[gene_index]]['mutation_rate']) != 0:
                    gene_alleles.append(self.mutate_genotype(reproducer_allele, gene_index))
                else:
                    gene_alleles.append(reproducer_allele)
            genotype.append(gene_alleles)
        descendant.genotype = genotype

    def generate_offspring_ancestry(self, reproducers, descendant):
        for generation in range(self.__ancestry_generations):
            if generation == 0:
                descendant.ancestry[generation] = [reproducers[0].id, reproducers[1].id]
            else:
                descendant.ancestry[generation] = reproducers[0].ancestry[generation - 1] + \
                                                  reproducers[1].ancestry[generation - 1]

    def discard_old_individuals(self):
        if self.__life_expectancy == 1 and self.__life_expectancy_sd == 0:
            new_groups = [[] for i in range(self.__group_number)]
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

    def reproduce(self):
        new_groups = self.discard_old_individuals()
        for group, group_index in zip(self.__groups, range(len(self.__groups))):
            group_size_normal = numpy.random.normal(self.__group_size, self.__group_size_sd)
            while len(new_groups[group_index]) < group_size_normal:
                reproducers = random.sample(group, 2)
                new_individual = Individual(self)
                self.generate_offspring_genotype(reproducers, new_individual)
                self.generate_offspring_ancestry(reproducers, new_individual)
                self.__newest_ind_id = new_individual.id
                new_groups[group_index].append(new_individual)
        self.__groups = new_groups

    def group_exchange(self):
        new_groups = [[] for i in range(self.__group_number)]
        for from_group_index in range(len(self.__groups)):
            emigrants = round(len(self.__groups[from_group_index]) * self.__group_migration)
            target_indexes = list(range(len(self.__groups)))
            target_indexes.remove(from_group_index)
            for i in range(emigrants):
                emigrant = self.__groups[from_group_index].pop(random.randrange(len(self.__groups[from_group_index])))
                target_group_index = random.choice(target_indexes)
                new_groups[target_group_index].append(emigrant)
        for group_index in range(len(self.__groups)):
            new_groups[group_index].extend(self.__groups[group_index])
        self.__groups = new_groups

    def save_generation_data(self, simulation_state):
        """ Possible alleles are added in a list, if more than 1 gene is considered all allele combinations
        will also be added """
        combined_phenotypes_list = list(self.__alleles_combinations_indexes.keys())
        # Proportion of individuals per phenotype before selection event
        if simulation_state == 0:
            # On the first generation, the simulation summary object is initialized
            if self.current_generation == 0:
                summary = [0 for i in range(self.generations)]
                for i in range(len(combined_phenotypes_list)):
                    self.__simulation_summary[0].append(summary.copy())
                    self.__simulation_summary[1].append(summary.copy())
            # Find index of individual's phenotype
            population_size = nested_len(self.__groups)
            for group in self.__groups:
                for individual in group:
                    match_indexes = list(range(len(combined_phenotypes_list)))
                    for phenotype in individual.phenotype:
                        new_match_indexes = []
                        for index in match_indexes:
                            combination = combined_phenotypes_list[index]
                            if re.search(r'(?:&|^)('+phenotype+')(?:&|$)', combination):
                                new_match_indexes.append(index)
                        match_indexes = new_match_indexes
                    self.__simulation_summary[1][match_indexes[0]][self.current_generation] += 1/population_size

        # Survivors per generation (after selection event)
        if simulation_state == 1:
            for individual in self.__survivors_per_generation:
                match_indexes = list(range(len(combined_phenotypes_list)))
                for phenotype in individual.phenotype:
                    new_match_indexes = []
                    for index in match_indexes:
                        allele = combined_phenotypes_list[index]
                        if re.search(r'(?:&|^)('+phenotype+')(?:&|$)', allele):
                            new_match_indexes.append(index)
                    match_indexes = new_match_indexes
                self.__simulation_summary[0][match_indexes[0]][self.current_generation] += 1

    def pass_generation(self):
        # self.save_generation_data(0)
        model_module.selection(self.groups)
        self.selection_event()
        self.migration()
        self.reproduce()
        self.group_exchange()
        # self.save_generation_data(1)
        self.__generation += 1


class Individual:
    def __init__(self, simulation):
        self.__simulation = simulation
        self.__age = 0
        self.__life_expectancy = round(numpy.random.normal(simulation.life_expectancy, simulation.life_expectancy_sd))
        if self.__life_expectancy == 0:
            self.__life_expectancy = 1
        self.__genotype = []
        self.__genes_properties = {}
        self.__phenotype = []
        self.__survival_probability = simulation.survival_probability
        self.__id = simulation.newest_ind_id + 1
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
    def genotype(self):
        return self.__genotype

    @genotype.setter
    def genotype(self, value):
        self.__genotype = value
        self.generate_phenotype()

    @property
    def phenotype(self):
        return self.__phenotype

    @phenotype.setter
    def phenotype(self, value):
        self.__phenotype = value

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

    # Calculates the phenotype considering the genotype and inheritance pattern
    def generate_phenotype(self):
        phenotype = []
        for i in range(len(self.__simulation.genes)):
            gene = self.__simulation.genes[i]
            # If both alleles are the same the phenotype is that allele
            if self.__genotype[i][0] == self.__genotype[i][1]:
                phenotype.append(self.__genotype[i][0])
            else:
                inheritance_divider = r'(?:>|=|\s|$|^)'
                characters_between_alleles = re.search(
                    f'{inheritance_divider}{self.__genotype[i][0]}{inheritance_divider}(.*){inheritance_divider}'
                    f'{self.__genotype[i][1]}{inheritance_divider}',
                    self.__simulation.model_config_dict[gene]['alleles'])
                reverse = False
                if not characters_between_alleles:
                    characters_between_alleles = re.search(
                        f'{inheritance_divider}{self.__genotype[i][1]}{inheritance_divider}(.*){inheritance_divider}'
                        f'{self.__genotype[i][0]}{inheritance_divider}',
                        self.__simulation.model_config_dict[gene]['alleles'])
                    reverse = True
                if '>' in characters_between_alleles.group(1):
                    if reverse:
                        phenotype.append(self.__genotype[i][1])
                    else:
                        phenotype.append(self.__genotype[i][0])
                else:
                    chosen_phenotype = self.__genotype[i][0] + '_' + self.__genotype[i][1]
                    # [gene][0] is equal to gene's alleles
                    if chosen_phenotype in self.__simulation.genes_properties[gene][0]:
                        phenotype.append(chosen_phenotype)
                    else:
                        phenotype.append(self.__genotype[i][1] + '_' + self.__genotype[i][0])
        self.phenotype = phenotype

    def age_individual(self):
        self.__age += 1
        return self.__age == self.__life_expectancy


class Genome:
    def __init__(self, genes):
        self.__genes = genes


class Gene:
    def __init__(self, chromosomes):
        self.__chromosomes = chromosomes


class Allele:
    def __init__(self, alleles):
        self.__genes = alleles


class SNVs:
    def __init__(self, phenotype):
        self.__phenotype = phenotype
        self.__snvs_list = []

    @property
    def snvs(self):
        return self.__snvs_list


def simulator_main():
    generations = int(general_config['simulation']['generations'])
    group_number = int(general_config['population']['group_number'])
    group_size = int(general_config['population']['group_size'])
    group_size_sd = float(general_config['population']['group_size_sd'])
    group_migration = float(general_config['population']['group_migration'])
    emigration = float(general_config['population']['emigration'])
    immigration = float(general_config['population']['immigration'])
    immigration_phenotype = general_config['population']['immigration_phenotype'].replace(' ', '').split(',')
    life_expectancy = int(general_config['population']['life_expectancy'])
    life_expectancy_sd = float(general_config['population']['life_expectancy_sd'])
    survival_probability = float(general_config['population']['survival_probability'])
    ancestry_generations = int(general_config['population']['ancestry_generations'])
    model_config_dict = {s: dict(model_config.items(s)) for s in model_config.sections()}

    simulation = Simulation(generations,
                            group_number,
                            group_size,
                            group_size_sd,
                            group_migration,
                            emigration,
                            immigration,
                            immigration_phenotype,
                            life_expectancy,
                            life_expectancy_sd,
                            survival_probability,
                            ancestry_generations,
                            model_config_dict)
    simulation.populate_groups()

    # Progress bar
    bar_msg = 'Simulation progress: '
    cols = get_terminal_size().columns - len(bar_msg)
    bar_char = '█'
    bar_end_chars = ' ▏▎▍▌▋▊▉'
    for i in range(simulation.generations):
        simulation.pass_generation()
        progress = cols*i/generations
        print('\033[K\r' + bar_msg + bar_char*int(progress) + bar_end_chars[int((progress - int(progress))*8)] +
              ' ' * (cols - int(progress) - 1), end='')

    # Rounds the values of the proportions of individuals per phenotype
    for phenotype_index in range(len(simulation.simulation_summary[1])):
        for generation_index in range(len(simulation.simulation_summary[1][phenotype_index])):
            simulation.simulation_summary[1][phenotype_index][generation_index] = \
                round(simulation.simulation_summary[1][phenotype_index][generation_index], 5)

    return simulation.simulation_summary, simulation.alleles_combinations_indexes, simulation.dict_allele_options


if __name__ == '__main__':
    simulator_main()
