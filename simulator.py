#!/usr/bin/python3

from os import path, get_terminal_size
import re
import random
import itertools

from configparser import ConfigParser

# Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

# Import model configuration
model_config = ConfigParser()
model_config.read(path.join('models', general_config['simulation']['model'] + '.ini'))
py_module = 'models.' + model_config['module']['name']
model = __import__(py_module, fromlist=[''])


def delete_random_elems(input_list, n):
    to_delete = set(random.sample(range(len(input_list)), n))
    return [x for i, x in enumerate(input_list) if i not in to_delete]


class Simulation:
    def __init__(
        self,
        population_size,
        generations,
        lifespan,
        emigration,
        immigration,
        immigration_phenotype,
        model,
        reproduction,
        selection_group_size,
        survival_range,
        altruism_cost_benefit,
        altruism_probability):
            self.__population_size = population_size
            self.__generations = generations
            self.__lifespan = lifespan
            self.__emigration = emigration
            self.__immigration = immigration
            self.__immigration_phenotype = immigration_phenotype
            self.__generation = 0
            self.__population = []
            self.__model = model
            self.__genes_properties = {}
            self.__reproduction = reproduction
            self.__selection_group_size = selection_group_size
            self.__groups = []
            self.__survival_range = survival_range
            self.__survival_rate_mean = 0
            self.__survivors_per_generation = []
            self.__altruism_cost_benefit = altruism_cost_benefit
            self.__altruism_probability = altruism_probability
            ''' First list of the summary is survivors per generation and second the number of individuals
            at the beginning of a generation '''
            self.__simulation_summary = [[],[]]
            self.__alleles_combinations_indexes = []
            # Number of digits for adding 0s at the beginning of the individuals ids
            self.__population_digits = len(str(self.__population_size))
            self.__generations_digits = len(str(self.__generations))
            self.__id_iter = itertools.count(start=1)

            for gene in model_config.sections():
                if gene != 'module':
                    incorrect_character = re.search('[^\w\s=>]', model_config[gene]['alleles'])
                    if incorrect_character:
                        raise Exception(f'Character "{incorrect_character.group(0)}" in allele not valid')
                    # Saving allele names for each gene
                    allele_options = re.split('\s*>\s*|\s*=\s*', model_config[gene]['alleles'])
                    noncodominant_allele_options = allele_options.copy()
                    # Saving codominant phenotypes for each gene
                    codominant_alleles_raw = re.findall('[\w=]*=[\w=]*', model_config[gene]['alleles'].replace(' ', ''))
                    for element in codominant_alleles_raw:
                        for combination in itertools.combinations(element.split('='), 2):
                            codominant_phenotype = combination[0]+'_'+combination[1]
                            allele_options.append(codominant_phenotype)
                    for allele in allele_options:
                        if re.search('\s', allele):
                            raise Exception(f'Missing inheritance symbol in "{allele}", allele names can\'t have spaces')
                    initial_frequencies = re.split('\s*,\s*', model_config[gene]['initial_frequencies'])
                    if len(initial_frequencies) != len(noncodominant_allele_options) and len(initial_frequencies) != len(noncodominant_allele_options)-1:
                        raise Exception(f'Incorrect number of allele frequencies given in gene "{gene}", there must be same number as alleles or one less')
                    # Initial frequencies are defined by a portion of the range from 0 to 1
                    # First frequency range is set manually as it needn't be calculated
                    initial_frequencies_ranges = [float(initial_frequencies[0])]
                    for i in range(len(initial_frequencies))[1:]:
                        new_range = float(initial_frequencies[i])+initial_frequencies_ranges[i-1]
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
                    allele_combination_name += element+'&'
                alleles_combinations_indexes[allele_combination_name[:-1]] = i
            self.__alleles_combinations_indexes = alleles_combinations_indexes

            if self.__immigration != 0:
                self.__immigrants_genotype = []
                if len(self.genes) != len(self.__immigration_phenotype):
                    raise Exception('Incorrect number of phenotypes specified for immigrants')
                for allele in self.__immigration_phenotype:
                    self.__immigrants_genotype.append([allele, allele])

    @property
    def population_size(self):
        return self.__population_size

    @property
    def generations(self):
        return self.__generations

    @property
    def current_generation(self):
        return self.__generation

    @property
    def population(self):
        return self.__population

    @population.setter
    def population(self,value):
        self.__population = value

    @property
    def genes(self):
        return list(self.__genes_properties.keys())

    @property
    def groups(self):
        return self.__groups

    @property
    def simulation_summary(self):
        return self.__simulation_summary

    @simulation_summary.setter
    def simulation_summary(self,value):
        self.__simulation_summary = value

    @property
    def genes_properties(self):
        return self.__genes_properties

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

    def populate(self):
        for i in range(self.__population_size):
            individual = Individual(self, next(self.__id_iter))
            individual.genotype = self.generate_individual_genotype()
            individual.generate_initial_id(self.current_generation, self.__population_digits, self.__generations_digits)
            self.__population.append(individual)

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
    
    # Generates a new individuals with the given immigrant's genotype
    def generate_immigrant(self):
        individual = Individual(self, next(self.__id_iter))
        individual.genotype = self.__immigrants_genotype
        # Although the individual has arrived in the current generation, it is considered it starts living in the next
        individual.generate_initial_id(self.current_generation + 1, self.__population_digits, self.__generations_digits)
        return individual
    
    def group_individuals(self):
        population_copy = self.__population.copy()
        groups = []
        # If the group size is higher that the population size there is just one group
        if len(population_copy) <= self.__selection_group_size:
            groups.append(population_copy)
            self.__groups = groups
        else:
            # Individuals assigned to a group will be removed from population_copy
            while len(population_copy) > self.__selection_group_size:
                group = []
                for i in range(self.__selection_group_size):
                    group.append(population_copy.pop(random.randint(0, len(population_copy) - 1)))
                groups.append(group)
            # The remaining individuals will be placed in random groups
            while len(population_copy) != 0:
                groups[random.randint(0, len(groups) - 1)].append(population_copy.pop(0))
            self.__groups = groups
    
    # Random survival probability based on the given range
    def assign_individuals_survival_probability(self):
        for individual in self.__population:
            individual.survival_probability = random.uniform(self.__survival_range[0], self.__survival_range[1])

    # Survival or death based on the survival probability of the individual. Represents foraging, predators, etc.
    def selection_event(self):
        self.__survivors_per_generation = []
        survival_probabilities = []

        for individual in self.__population:
            survival_probabilities.append(individual.survival_probability)
            picker = random.random()
            if picker < individual.survival_probability:
                self.__survivors_per_generation.append(individual)

        self.population = self.__survivors_per_generation
        self.__survival_rate_mean = sum(survival_probabilities)/len(survival_probabilities)

    def reproduce(self):
        new_population = []
        if self.__emigration != 0:
            self.__population = delete_random_elems(self.__population, round(len(self.__population)*self.__emigration))
        if self.__lifespan <= 1:
            # All previous individuals die
            missing_population = self.__population_size
        else:
            # Old individuals die
            survivors_population = [individual for individual in self.__population if individual.age < self.__lifespan]
            missing_population = self.__population_size - len(survivors_population)
        if self.__immigration != 0:
            immigrants = min(round(self.__population_size * self.__immigration), missing_population)
            for i in range(immigrants):
                new_population.append(self.generate_immigrant())

        # Panmixia
        if self.__reproduction == 0:
            while len(new_population) < missing_population:
                reproducers = random.sample(self.__population, 2)
                new_individual_genotype = []
                # For each gene there's a probability of mutation
                for i in range(len(reproducers[0].genotype)):
                    mutation_rate = float(model_config[self.genes[i]]['mutation_rate'])
                    gene_alleles = []
                    # Chose one allele of each gene of each reproducer
                    for reproducer_index in range(2):
                        reproducer_allele = random.choice(reproducers[reproducer_index].genotype[i])
                        # If a mutations occurs, chose from random another allele of the gene
                        if random.random() < mutation_rate:
                            allele_possibilities = self.__genes_properties[self.genes[i]][0].copy()
                            allele_possibilities.remove(reproducer_allele)
                            reproducer_allele = random.choice(allele_possibilities)
                        gene_alleles.append(reproducer_allele)
                    new_individual_genotype.append(gene_alleles)
                new_individual = Individual(self, next(self.__id_iter))
                new_individual.genotype = new_individual_genotype
                new_individual.generate_id(reproducers, self.__population_digits, self.__generations_digits)
                new_population.append(new_individual)

        else:
            if self.__reproduction > len(self.__population[0].phenotype):
                raise Exception('Number of identical alleles needed for reproduction cannot exceed number of genes')
            count = 0
            while len(new_population) < missing_population:
                possible_reproducers = self.__population.copy()
                first_reproducer = random.choice(possible_reproducers)
                possible_reproducers.remove(first_reproducer)
                same_phenotype_count = 0
                while same_phenotype_count < self.__reproduction:
                    if len(possible_reproducers) == 0:
                        break
                    same_phenotype_count = 0
                    second_reproducer = random.choice(possible_reproducers)
                    for phenotype_index in range(len(first_reproducer.phenotype)):
                        if first_reproducer.phenotype[phenotype_index] == second_reproducer.phenotype[phenotype_index]:
                            same_phenotype_count += 1
                    possible_reproducers.remove(second_reproducer)
                if same_phenotype_count >= self.__reproduction:
                    new_individual_genotype = []
                    reproducers = [first_reproducer, second_reproducer]
                    for i in range(len(reproducers[0].genotype)):
                        mutation_rate = float(model_config[self.genes[i]]['mutation_rate'])
                        gene_alleles = []
                        for reproducer_index in range(2):
                            reproducer_allele = random.choice(reproducers[reproducer_index].genotype[i])
                            if random.random() < mutation_rate:
                                # [gene][0] is equal to gene's alleles
                                allele_possibilities = self.__genes_properties[self.genes[i]][0].copy()
                                allele_possibilities.remove(reproducer_allele)
                                reproducer_allele = random.choice(allele_possibilities)
                            gene_alleles.append(reproducer_allele)
                        new_individual_genotype.append(gene_alleles)
                    new_individual = Individual(self, next(self.__id_iter))
                    new_individual.genotype = new_individual_genotype
                    new_population.append(new_individual)
                count += 1

        # Replace all the previous population if the lifespan is 1
        if self.__lifespan <= 1:
            self.__population = new_population
        else:
            # Increase age in individuals of existent population
            for individual in survivors_population:
                individual.age_individual()
            # Add the newborn individuals
            survivors_population.extend(new_population)
            self.__population = survivors_population

    def save_generation_data(self, simulation_state):
        ''' Possible alleles are added in a list, if more than 1 gene is considered all allele combinations
        will also be added '''
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
            for individual in self.__population:
                match_indexes = list(range(len(combined_phenotypes_list)))
                for phenotype in individual.phenotype:
                    new_match_indexes = []
                    for index in match_indexes:
                        combination = combined_phenotypes_list[index]
                        if re.search('(?:&|^)('+phenotype+')(?:&|$)', combination):
                            new_match_indexes.append(index)
                    match_indexes = new_match_indexes
                self.__simulation_summary[1][match_indexes[0]][self.current_generation] += 1/self.__population_size

        # Survivors per generation (after selection event)
        if simulation_state == 1:
            for individual in self.__survivors_per_generation:
                match_indexes = list(range(len(combined_phenotypes_list)))
                for phenotype in individual.phenotype:
                    new_match_indexes = []
                    for index in match_indexes:
                        allele = combined_phenotypes_list[index]
                        if re.search('(?:&|^)('+phenotype+')(?:&|$)', allele):
                            new_match_indexes.append(index)
                    match_indexes = new_match_indexes
                self.__simulation_summary[0][match_indexes[0]][self.current_generation] += 1

    def pass_generation(self):
        # Data at the start of the generation
        self.save_generation_data(0)
        self.assign_individuals_survival_probability()
        self.group_individuals()
        model.selection(self.groups)
        self.selection_event()
        ''' Iterative counter for individual's ids. The number is reset before generating the individuals
        of the next generation, immigrants or newborns '''
        self.__id_iter = itertools.count(start=1)
        self.reproduce()
        # Data after selection event
        self.save_generation_data(1)
        self.__generation += 1


class Individual:
    def __init__(self, simulation, id_num):
        self.__simulation = simulation
        self.__age = 0
        self.__genotype = []
        self.__genes_properties = {}
        self.__phenotype = []
        self.__survival_probability = 0
        self.__id = id_num

    @property
    def age(self):
        return self.__age

    @property
    def genotype(self):
        return self.__genotype

    @genotype.setter
    def genotype(self, value):
        self.__genotype = value
        self.choose_phenotype()

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

    # Calculates the phenotype considering the genotype and inheritance pattern
    def choose_phenotype(self):
        phenotype = []
        for i in range(len(self.__simulation.genes)):
            gene = self.__simulation.genes[i]
            # If both alleles are the same the phenotype is that allele
            if self.__genotype[i][0] == self.__genotype[i][1]:
                phenotype.append(self.__genotype[i][0])
            else:
                inheritance_divider = '(?:>|=|\s|$|^)'
                characters_between_alleles = re.search(
                    f'{inheritance_divider}{self.__genotype[i][0]}{inheritance_divider}(.*){inheritance_divider}{self.__genotype[i][1]}{inheritance_divider}',
                    model_config[gene]['alleles'])
                reverse = False
                if not characters_between_alleles or characters_between_alleles == None:
                    characters_between_alleles = re.search(
                        f'{inheritance_divider}{self.__genotype[i][1]}{inheritance_divider}(.*){inheritance_divider}{self.__genotype[i][0]}{inheritance_divider}',
                        model_config[gene]['alleles'])
                    reverse = True
                if '>' in characters_between_alleles.group(1):
                    if reverse:
                        phenotype.append(self.__genotype[i][1])
                    else:
                        phenotype.append(self.__genotype[i][0])
                else:
                    chosen_phenotype = self.__genotype[i][0]+'_'+self.__genotype[i][1]
                    # [gene][0] is equal to gene's alleles
                    if chosen_phenotype in self.__simulation.genes_properties[gene][0]:
                        phenotype.append(chosen_phenotype)
                    else:
                        phenotype.append(self.__genotype[i][1]+'_'+self.__genotype[i][0])
        self.phenotype = phenotype

    # Generates the id for individuals with unknown parents
    def generate_initial_id(self, generation, population_digits, generations_digits):
        # First individuals have unknown parents represented with an id with only 0s
        parents_id = '0' * (population_digits + generations_digits) * 2
        generation_id = str(generation).zfill(generations_digits)
        population_id = str(self.__id).zfill(population_digits)
        ind_id = parents_id + generation_id + population_id
        self.__id = ind_id

    # Generates the id for individuals with known parents
    def generate_id(self, parents, population_digits, generations_digits):
        # The generation and id number of each parent is saved
        parents_id = parents[0].id[(generations_digits + population_digits) * 2:] + \
                     parents[1].id[(generations_digits + population_digits) * 2:]
        # Although the individual is born in the current generation, it is considered it starts living in the next one
        generation_id = str(self.__simulation.current_generation + 1).zfill(generations_digits)
        population_id = str(self.__id).zfill(population_digits)
        ind_id = parents_id + generation_id + population_id
        self.__id = ind_id

    def age_individual(self):
        self.__age += 1


def simulator_main():
    generations = int(general_config['simulation']['generations'])+1
    population_size = int(general_config['population']['size'])
    lifespan = int(general_config['population']['lifespan'])
    emigration = float(general_config['population']['emigration'])
    immigration = float(general_config['population']['immigration'])
    immigration_phenotype = general_config['population']['immigration_phenotype'].replace(' ','').split(',')
    reproduction = int(general_config['population']['reproduction'])
    selection_group_size = int(general_config['population']['selection_group_size'])
    survival_range = [float(perc) for perc in re.split(',\s*', general_config['population']['survival_range'])]
    altruism_cost_benefit = [float(general_config['population']['altruism_cost']), float(general_config['population']['altruism_benefit'])]
    altruism_probability = float(general_config['population']['altruism_probability'])
    model = general_config['simulation']['model']

    simulation = Simulation(population_size, generations, lifespan,
                            emigration, immigration, immigration_phenotype,
                            model, reproduction, selection_group_size, survival_range,
                            altruism_cost_benefit, altruism_probability)
    simulation.populate()

    # Progress bar
    bar_msg = 'Simulation progress: '
    cols = get_terminal_size().columns - len(bar_msg)
    bar_char = '█'
    bar_end_chars = ' ▏▎▍▌▋▊▉'
    for i in range(generations):
        simulation.pass_generation()
        progress = cols*i/generations
        # print('\033[K\r' + bar_msg + bar_char*int(progress) + bar_end_chars[int((progress-int(progress))*8)] + ' ' * (cols - int(progress) - 1), end='')

    # Rounds the values of the proportions of individuals per phenotype
    for phenotype_index in range(len(simulation.simulation_summary[1])):
        for generation_index in range(len(simulation.simulation_summary[1][phenotype_index])):
            simulation.simulation_summary[1][phenotype_index][generation_index] = round(simulation.simulation_summary[1][phenotype_index][generation_index], 5)
    return simulation.simulation_summary, simulation.alleles_combinations_indexes, simulation.dict_allele_options

if __name__ == '__main__':
    simulator_main()