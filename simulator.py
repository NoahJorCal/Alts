#!/usr/bin/python3
import time
import warnings
from os import path, get_terminal_size
import re
import random
import bisect

import h5py
import numpy as np
import pandas as pd
from scipy.stats import poisson
import itertools
import argparse
from configparser import ConfigParser
import os

last_max_group = 5
# global_survival_prob = []
# global_surv = []

if __name__ == '__main__':
    # Argument parser
    parser = argparse.ArgumentParser(prog='Alts simulator',
                                     description='Main simulation of Alts program')
    parser.add_argument('-d', '--directory', default='output', help='Output directory where individual'
                                                                    'simulation results will be stored')
    parser.add_argument('-o', '--output', default='simulation.h5', help='Output file where the'
                                                                        'simulation results will be stored')
    parser.add_argument('-s', '--seed', required=False, type=int, help='Set a seed for the simulation')

    args = parser.parse_args()
    global_output_file_name = args.output
    if '.h5' not in global_output_file_name and '.hdf5' not in global_output_file_name and \
       '.h5p' not in global_output_file_name and '.he5' not in global_output_file_name and \
       '.h5m' not in global_output_file_name and '.h5z' not in global_output_file_name:
        global_output_file_name += '.h5'
    seed = args.seed

# Import general configuration
general_config = ConfigParser()
config_path = os.path.join(os.path.dirname(__file__), 'config.ini')
general_config.read(config_path)

# Import model configuration
model_config = ConfigParser()
model_config.read(path.join(os.path.dirname(__file__), 'models', general_config['simulation']['model'] + '.ini'))
py_module = 'models.' + model_config['module']['name']
model_module = __import__(py_module, fromlist=[''])


def print_name_type(name, obj):
    print(name, type(obj))


def uniquify(path_name):
    filename, extension = os.path.splitext(path_name)
    counter = 0

    while os.path.exists(path_name):
        path_name = filename + "_" + str(counter) + extension
        counter += 1

    return path_name


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


class Chromosome:
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
        mean_mutations = mutation_rate * locus_size
        n_mutations = poisson.rvs(mean_mutations)
        mutations = np.sort([random.random() for _ in range(n_mutations)])
        mutations_indexes = np.searchsorted(self.__snvs, mutations)
        self.__snvs = np.insert(self.__snvs, mutations_indexes, mutations)
        if not list(np.sort(self.__snvs)) == list(self.__snvs):
            print(np.sort(self.__snvs), self.__snvs)
            raise Exception

    def snvs_to_positions(self):
        return np.array([snv.position for snv in self.__snvs])

    def snvs_to_sequence(self, locus_size):
        mutations = [round(snv * locus_size) for snv in self.__snvs]
        sequence = np.zeros(locus_size)
        for mutation in mutations:
            if mutation != locus_size:
                if sequence[mutation] == 1 and mutation + 1 != locus_size:
                    sequence[mutation + 1] = 1
                else:
                    sequence[mutation] = 1
        return sequence


class Locus:
    def __init__(self, chromosomes, name: str, locus_size: int, mutation_rate: float, recombination_rate: float):
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
        mean_crossovers = self.__recombination_rate * self.__locus_size
        if not crossovers:
            number_crossovers = poisson.rvs(mean_crossovers)
            crossovers = [random.random() for _ in range(number_crossovers)]
        for crossover_point in crossovers:
            crossover_index_0 = bisect.bisect_left(self.__chromosomes[0].snvs, crossover_point)
            crossover_index_1 = bisect.bisect_left(self.__chromosomes[1].snvs, crossover_point)
            chromosome_0 = np.concatenate((
                self.__chromosomes[0].snvs[:crossover_index_0],
                self.__chromosomes[1].snvs[crossover_index_1:]))
            chromosome_1 = np.concatenate((
                self.__chromosomes[1].snvs[:crossover_index_1],
                self.__chromosomes[0].snvs[crossover_index_0:]))
            self.__chromosomes[0].snvs = chromosome_0
            self.__chromosomes[1].snvs = chromosome_1


class Genome:
    def __init__(self, loci):
        self.__loci = loci

    @property
    def loci(self):
        loci = []
        for locus in self.__loci:
            loci.append(locus)
        return loci

    @property
    def locus_names(self):
        locus_names = []
        for locus in self.__loci:
            locus_names.append(locus.name)
        return locus_names

    def haplotype(self):
        full_haplotype = []
        for locus in self.__loci:
            locus_haplotype = []
            for chromosome in locus.chromosomes:
                locus_haplotype.append(chromosome.snvs_to_sequence(locus.locus_size))
            full_haplotype.append(locus_haplotype)
        return tuple(full_haplotype)


class Individual:
    def __init__(self, simulation, initial_survival_probability):
        self.__simulation = simulation
        self.__age = 0
        self.__life_expectancy = round(np.random.normal(simulation.life_expectancy, simulation.life_expectancy_sd))
        if self.__life_expectancy == 0:
            self.__life_expectancy = 1

        loci = []
        for locus in simulation.loci:
            locus_size = int(simulation.model_config_dict[locus]['locus_size'])
            mutation_rate = float(simulation.model_config_dict[locus]['mutation_rate'])
            recombination_rate = float(simulation.model_config_dict[locus]['recombination_rate'])
            loci.append(Locus([], locus, locus_size, mutation_rate, recombination_rate))
        self.__genome = Genome(loci)

        self.__genotype = []
        self.__phenotype = []
        self.__initial_survival_probability = initial_survival_probability
        self.__survival_probability = self.__initial_survival_probability
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
        self.__genotype = value
        self.generate_phenotype()

    @property
    def phenotype(self):
        return self.__phenotype

    @phenotype.setter
    def phenotype(self, value):
        self.__phenotype = value

    @property
    def initial_survival_probability(self):
        return self.__initial_survival_probability

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
        for i in range(len(self.__simulation.loci)):
            locus = self.__simulation.loci[i]
            # If both alleles are the same the phenotype is that allele
            if self.__genotype[i][0] == self.__genotype[i][1]:
                phenotype.append(self.__genotype[i][0])
            else:
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
                if '>' in characters_between_alleles.group(1):
                    if reverse:
                        phenotype.append(self.__genotype[i][1])
                    else:
                        phenotype.append(self.__genotype[i][0])
                else:
                    chosen_phenotype = self.__genotype[i][0] + '_' + self.__genotype[i][1]
                    # [locus][0] is equal to locus's alleles
                    if chosen_phenotype in self.__simulation.loci_properties[locus][0]:
                        phenotype.append(chosen_phenotype)
                    else:
                        phenotype.append(self.__genotype[i][1] + '_' + self.__genotype[i][0])
        self.phenotype = phenotype

    def generate_genome(self):
        for locus, genotype in zip(self.__genome.loci, self.__genotype):
            chromosomes = []
            for allele in genotype:
                chromosomes.append(Chromosome(allele))
            locus.chromosomes = chromosomes

    def age_individual(self):
        self.__age += 1
        return self.__age == self.__life_expectancy


class Simulation:
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
            model_config_dict,
            output_file_name='simulation_output.h5'):
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
        self.__newest_ind_id = 0
        self.__output_file_name = output_file_name
        self.__stop = False

        self.last_max_group = 5

        for locus in self.__model_config_dict.keys():
            if locus != 'module':
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

    # Randomize individual's genotype based on initial locus frequencies
    def generate_individual_genotype(self):
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
        with h5py.File(self.__output_file_name, 'w') as _:
            pass
        for group_i in range(self.__group_number):
            group = []
            for ind_i in range(self.__group_size):
                survival_probability = np.random.normal(self.__survival_probability_mean,
                                                        self.__survival_probability_sd)
                individual = Individual(self, survival_probability)
                individual.genotype = self.generate_individual_genotype()
                individual.generate_genome()
                self.__newest_ind_id = individual.id
                group.append(individual)
            self.__groups.append(group)

    def reset_survival_prob(self):
        for group in self.__groups:
            for ind in group:
                ind.survival_probability = ind.initial_survival_probability

    # Survival or death based on the survival probability of the individual. Represents foraging, predators, etc.
    def selection_event(self):
        alt_perc = []
        for group in self.__groups:
            alts = 0
            for ind in group:
                if 'altruistic' in ind.phenotype[0]:
                    alts += 1
            alt_perc.append(round(alts/len(group)*100))
        # print('alt perc', alt_perc)
        with h5py.File(self.__output_file_name, 'a') as f:
            survivors = []
            survived = np.zeros(0)
            survived_list_index = 0
            s = []
            for group in self.__groups:
                survivors_groups = []
                new_survived_size = survived_list_index + len(group)
                survived.resize((new_survived_size,), refcheck=False)
                for ind_i in range(len(group)):
                    picker = random.random()
                    if picker < group[ind_i].survival_probability:
                        survivors_groups.append(group[ind_i])
                        survived[survived_list_index + ind_i] = 1
                    else:
                        survived[survived_list_index + ind_i] = 0
                survivors.append(survivors_groups)
                s.append(round(len(survivors_groups)/len(group)*100))
                # global_surv.append(round(len(survivors_groups)/len(group)*100))
                survived_list_index += len(group)
            # print(self.__groups_indices)
            # print('survival_prob', list(np.around(self.__calc_avg_survival_prob, decimals=2)))
            # print('survi', s)
            # print()
            self.__groups = survivors
            generation_group = f[f'/generation_{str(self.current_generation).zfill(len(str(self.__generations)))}']
            generation_group.create_dataset('survivors', data=survived)

    def save_avg_survival_prob(self):
        calc_avg_survival_prob = []
        exp_avg_survival_prob = []
        for group in self.__groups:
            calc_survival_prob = np.zeros(len(group))
            exp_survival_prob = np.zeros(len(group))
            for ind_i in range(len(group)):
                calc_survival_prob[ind_i] = group[ind_i].survival_probability
                exp_survival_prob[ind_i] = group[ind_i].initial_survival_probability
            calc_avg_survival_prob.append(np.mean(calc_survival_prob))
            exp_avg_survival_prob.append(np.mean(exp_survival_prob))
            # global_survival_prob.append(np.mean(survival_prob))
        self.__calc_avg_survival_prob = calc_avg_survival_prob
        self.__exp_avg_survival_prob = exp_avg_survival_prob

    # Generates a new individuals with the given immigrant's genotype
    def generate_immigrant(self):
        survival_probability = np.random.normal(self.__survival_probability_mean,
                                                self.__survival_probability_sd)
        individual = Individual(self, survival_probability)
        individual.genotype = self.__immigrants_genotype
        individual.generate_genome()
        self.__newest_ind_id = individual.id
        return individual

    def migration(self):
        full_population_size = self.__group_size * len(self.__groups)
        real_population_size = nested_len(self.__groups)

        if self.__emigration != 0:
            self.__groups = delete_random_elems(self.__groups, round(real_population_size * self.__emigration))

        missing_population = full_population_size - real_population_size
        if self.__immigration != 0:
            immigrants = min(round(full_population_size * self.__immigration), missing_population)
            for i in range(immigrants):
                group_index = random.choice(range(len(self.__groups)))
                self.__groups[group_index].append(self.generate_immigrant())

    def generate_offspring_genome(self, reproducers, descendant):
        # Names for altruistic step
        genotype = []
        # Actual Locus objects
        loci = []
        for locus_index in range(len(reproducers[0].genome.loci)):
            genotype.append([])
            chromosomes = []
            # For each locus, each reproducer contribute with one allele (Chromosome object)

            for reproducer in reproducers:
                locus = reproducer.genome.loci.copy()[locus_index]
                if locus.recombination_rate != 0:
                    locus.recombine()
                # The chromosome that contains the inherited locus is selected at random
                chromosome = random.choice(locus.chromosomes)
                if locus.mutation_rate != 0:
                    if locus.mutation_rate * locus.locus_size == 0:
                        warnings.warn(f'Mutation rate times locus size is equal to 0 in locus {locus.name}, '
                                      f'there will be no mutations')
                    chromosome.mutate(locus.locus_size, locus.mutation_rate)
                chromosomes.append(chromosome)
                genotype[locus_index].append(chromosome.allele)

            name = self.loci[locus_index]
            locus_size = reproducers[0].genome.loci[locus_index].locus_size
            mutation_rate = reproducers[0].genome.loci[locus_index].mutation_rate
            recombination_rate = reproducers[0].genome.loci[locus_index].recombination_rate

            loci.append(Locus(chromosomes, name, locus_size, mutation_rate, recombination_rate))
        descendant.genotype = genotype
        descendant.genome = Genome(loci)

    def generate_offspring_ancestry(self, reproducers, descendant):
        for generation in range(self.__ancestry_generations):
            if generation == 0:
                descendant.ancestry[generation] = [reproducers[0].id, reproducers[1].id]
            else:
                descendant.ancestry[generation] = reproducers[0].ancestry[generation - 1] + \
                                                  reproducers[1].ancestry[generation - 1]

    def discard_old_individuals(self):
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

    def delete_empty_groups(self):
        for group_i in range(len(self.__groups) - 1, - 1, - 1):
            if len(self.__groups[group_i]) < 2:
                self.__groups.pop(group_i)
                np.delete(self.__groups_indices, group_i)

    def divide_big_groups(self, groups):
        split_groups = []
        new_groups_indices = []
        for group_i in range(len(groups)):
            if len(groups[group_i]) > self.__group_size_limit:
                # print(f'group {self.__groups_indices[group_i]} has been split')
                np.delete(self.__groups_indices, group_i)
                # Individuals are added in the group by age, the group is shuffled to avoid grouping
                # older individuals with older individuals and younger with younger
                random.shuffle(groups[group_i])
                # First half of the group is added as another group
                split_groups.append(groups[group_i][:round(len(groups[group_i]) / 2)].copy())
                new_groups_indices.append(self.__groups_indices[group_i])
                # The original group will be the second half
                split_groups.append(groups[group_i][round(len(groups[group_i]) / 2):].copy())
                new_groups_indices.append(max(max(self.__groups_indices), max(new_groups_indices)) + 1)
                self.__total_groups += 1
            else:
                split_groups.append(groups[group_i])
                new_groups_indices.append(self.__groups_indices[group_i])
        # print(new_groups_indices)
        self.__groups_indices = np.array(new_groups_indices)
        self.__groups = split_groups

    def reproduce(self):
        self.delete_empty_groups()
        print([len(group) for group in self.__groups])
        print([len(group) + len(group) * self.__descendants_per_survivor
                               for group in self.__groups])
        new_groups = self.discard_old_individuals()
        offspring = [[] for _ in range(len(self.__groups))]
        # print('before', [len(group) + len(group) * self.__descendants_per_survivor for group in self.__groups])
        # Survivors + descendants per survivor
        exp_group_sizes = [len(group) + len(group) * self.__descendants_per_survivor
                           for group in self.__groups]
        calc_exp_sp_ratios = [self.__calc_avg_survival_prob[group_i] / self.__exp_avg_survival_prob[group_i]
                              for group_i in range(len(self.__groups))]
        print(calc_exp_sp_ratios)
        mean_ratios = sum(calc_exp_sp_ratios) / len(calc_exp_sp_ratios)
        print(mean_ratios)
        calc_exp_sp_ratios = [(calc_exp_sp_ratios[group_i] - mean_ratios) * 2 + 1
                              for group_i in range(len(self.__groups))]
        print(calc_exp_sp_ratios)
        new_group_sizes = [round(exp_group_sizes[group_i] * calc_exp_sp_ratios[group_i])
                           for group_i in range(len(self.__groups))]
        print(new_group_sizes)
        # print('after', new_group_sizes)
        # print()
        # new_group_sizes = [len(group) + len(group) * self.__descendants_per_survivor for group in self.__groups]
        if sum(new_group_sizes) > self.__population_size_limit:
            new_group_sizes = [round(group_size * self.__population_size_limit / sum(new_group_sizes))
                               for group_size in new_group_sizes]
        for group_index in range(len(self.__groups)):
            group_size = new_group_sizes[group_index] - len(new_groups[group_index])

            while len(offspring[group_index]) < group_size:
                reproducers = random.sample(self.__groups[group_index], 2)
                survival_probability = (reproducers[0].initial_survival_probability +
                                        reproducers[1].initial_survival_probability) / 2
                new_individual = Individual(self, survival_probability)
                self.generate_offspring_genome(reproducers, new_individual)
                self.generate_offspring_ancestry(reproducers, new_individual)
                self.__newest_ind_id = new_individual.id
                # noinspection PyTypeChecker
                offspring[group_index].append(new_individual)
        new_groups = self.group_exchange(new_groups, offspring, new_group_sizes)
        self.divide_big_groups(new_groups)

    def group_exchange(self, survivors, newborns, group_sizes):
        if len(survivors) > 1:
            # new_groups = [[] for _ in range(len(self.__groups))]
            for from_group_index in range(len(survivors)):
                # Number of emigrants of the group
                emigrants = round(group_sizes[from_group_index] * self.__group_migration)
                target_indexes = list(range(len(survivors)))
                # Possible target groups
                target_indexes.remove(from_group_index)
                if len(newborns[from_group_index]) != 0:
                    if len(newborns[from_group_index]) > emigrants:
                        for i in range(emigrants):
                            emigrant = newborns[from_group_index].pop(random.randrange(len(newborns[from_group_index])))
                            target_group_index = random.choice(target_indexes)
                            survivors[target_group_index].append(emigrant)
                    else:
                        emigrants_left = emigrants - len(newborns[from_group_index])
                        for i in range(len(newborns[from_group_index])):
                            emigrant = newborns[from_group_index].pop(random.randrange(len(newborns[from_group_index])))
                            target_group_index = random.choice(target_indexes)
                            survivors[target_group_index].append(emigrant)
                        for i in range(emigrants_left):
                            emigrant = survivors[from_group_index].pop(random.randrange(len(survivors[from_group_index])))
                            target_group_index = random.choice(target_indexes)
                            survivors[target_group_index].append(emigrant)
                elif len(survivors[from_group_index]) != 0:
                    for i in range(emigrants):
                        emigrant = survivors[from_group_index].pop(random.randrange(len(survivors[from_group_index])))
                        target_group_index = random.choice(target_indexes)
                        survivors[target_group_index].append(emigrant)
            for group_i in range(len(survivors)):
                survivors[group_i].extend(newborns[group_i])
            return survivors

    def save_generation_data(self):
        with h5py.File(self.__output_file_name, 'a') as f:
            # The h5 file will have a group for each generation and inside that group,
            # a group for each group of the population
            generation_group = f.create_group(
                f'generation_{str(self.current_generation).zfill(len(str(self.__generations)))}')
            loci_alleles = []
            for locus_i in range(len(self.loci)):
                loci_alleles.append(self.__loci_properties[self.loci[locus_i]][0])
            groups = np.zeros(0)
            phenotypes = np.zeros(0)
            genotypes = []
            for _ in self.loci:
                genotypes.append(np.zeros(0))
            individuals_list_index = 0
            alleles_list_index = 0
            # print(self.__groups_indices)
            # print(np.arange(len(self.__groups)))
            # print()
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
                    phenotypes[individuals_list_index + ind_i] = self.__alleles_combinations.index(
                        '&'.join(group[ind_i].phenotype))
                    groups[individuals_list_index + ind_i] = self.__groups_indices[self.__groups.index(group)]
                    for allele_i in range(2):
                        for locus_i in range(len(self.loci)):
                            genotypes[locus_i][alleles_list_index + genotypes_index] = loci_alleles[locus_i].index(
                                group[ind_i].genotype[locus_i][allele_i])
                        genotypes_index += 1
                individuals_list_index += len(group)
                alleles_list_index += len(group) * 2
            altruistic_indexes = [i for i, item in enumerate(self.__alleles_combinations) if re.search(
                r'altruistic', item)]
            altruists = 0
            for altruistic_index in altruistic_indexes:
                altruists += np.count_nonzero(phenotypes == altruistic_index)
            # if not altruists:
            #     self.__stop = True
            # if self.last_max_group > round(max(groups)) + 1:
            #     raise Exception(self.current_generation)
            self.last_max_group = round(max(groups)) + 1
            generation_group.create_dataset('group', data=groups)
            generation_group.create_dataset('phenotype', data=phenotypes)
            generation_group['phenotype'].attrs['phenotype_names'] = self.__alleles_combinations
            for locus_i in range(len(self.loci)):
                generation_group.create_dataset(f'locus_{self.loci[locus_i]}', data=genotypes[locus_i])

    def pass_generation(self):
        print(self.current_generation)
        self.reset_survival_prob()
        model_module.selection(self.groups)
        self.save_avg_survival_prob()
        self.selection_event()
        self.migration()
        self.reproduce()
        # self.group_exchange()
        self.__generation += 1
        self.save_generation_data()

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

        with h5py.File(self.__output_file_name, 'a') as f:
            f.create_group('haplotypes')
            for locus_i in range(len(self.loci)):
                f.create_dataset(f'haplotypes/haplotypes_{self.loci[locus_i]}', data=haplotypes[locus_i])

    def save_snvs(self):
        snvs_lists = [[] for _ in range(len(self.loci))]
        for group in self.__groups:
            for individual in group:
                loci = individual.genome.loci
                for locus_index in range(len(self.loci)):
                    for chromosome in loci[locus_index].chromosomes:
                        snvs_lists[locus_index].append(chromosome.snvs)
        snvs_dir = os.path.join(os.path.dirname(__file__), 'snvs')
        os.makedirs(snvs_dir, exist_ok=True)
        for locus_i in range(len(self.loci)):
            file_name = f'{self.loci[locus_i]}_snvs.csv'
            result = os.path.join(os.path.dirname(__file__), 'snvs', file_name)
            result = uniquify(result)
            locus_array = snvs_lists[locus_i]
            pd.DataFrame(locus_array).to_csv(result, header=False)


def simulator_main(output_dir, output_file, sim_seed=None, quiet=False):
    if sim_seed:
        random.seed(sim_seed)
        np.random.seed(sim_seed)
    os.makedirs(output_dir, exist_ok=True)
    if '.h5' not in output_file and '.hdf5' not in output_file and \
       '.h5p' not in output_file and '.he5' not in output_file and \
       '.h5m' not in output_file and '.h5z' not in output_file:
        output_file += '.h5'
    output_file_name = os.path.join(os.path.dirname(__file__), output_dir, output_file)
    output_file_name = uniquify(output_file_name)
    # if os.path.exists(output_file_name):
    #     raise Exception(f'File {output_file_name} already exists')

    generations = int(general_config['simulation']['generations'])
    group_number = int(general_config['population']['group_number'])
    group_size = int(general_config['population']['group_size'])
    group_size_limit = int(general_config['population']['group_size_limit'])
    population_size_limit = int(general_config['population']['population_size_limit'])
    descendants_per_survivor = float(general_config['population']['descendants_per_survivor'])
    group_migration = float(general_config['population']['group_migration'])
    emigration = float(general_config['population']['emigration'])
    immigration = float(general_config['population']['immigration'])
    immigration_phenotype = general_config['population']['immigration_phenotype']
    life_expectancy = int(general_config['population']['life_expectancy'])
    life_expectancy_sd = float(general_config['population']['life_expectancy_sd'])
    survival_probability_mean = float(general_config['population']['survival_probability_mean'])
    survival_probability_sd = float(general_config['population']['survival_probability_sd'])
    ancestry_generations = int(general_config['population']['ancestry_generations'])

    model_config_dict = {s: dict(model_config.items(s)) for s in model_config.sections()}

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
                            model_config_dict,
                            output_file_name)
    simulation.populate_groups()
    simulation.save_generation_data()

    generations_duration = np.zeros(generations)
    generation_start = time.perf_counter()

    # Progress bar
    bar_msg = 'Simulation progress: '
    cols = get_terminal_size().columns - len(bar_msg)
    bar_char = '█'
    bar_end_chars = ' ▏▎▍▌▋▊▉'
    for i in range(generations):
        simulation.pass_generation()
        if simulation.stop:
            return False
        progress = cols*i/generations
        if not quiet:
            print('\r' + bar_msg + bar_char*int(progress) + bar_end_chars[int((progress - int(progress))*8)] +
                  ' ' * (cols - int(progress) - 1), end='')
        generation_end = time.perf_counter()
        generations_duration[i] = generation_end - generation_start
        generation_start = time.perf_counter()

    # simulation.save_haplotypes()
    simulation.save_snvs()

    # with open('survivors.txt', 'w') as f:
    #     f.write(str(global_survival_prob))
    #     f.write('\n')
    #     f.write(str(global_surv))
    #     f.write('\n')

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
        print(simulation.alleles_combinations)
        f.attrs['phenotypes_names'] = str(phenotypes_names)
        print(str(phenotypes_names))
        f.attrs['alleles_names'] = str(alleles_names)
        str(alleles_names)
        f.create_dataset('duration', data=generations_duration)

    return output_file_name


if __name__ == '__main__':
    # Across-threads counter for output file
    simulator_main(args.directory, args.output, args.seed, args.quiet)
