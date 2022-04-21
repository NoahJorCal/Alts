#!/usr/bin/python3

from os import path
import re
import random
from itertools import combinations
from itertools import product

from configparser import ConfigParser

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

#Import model configuration
model_config = ConfigParser()
model_config.read(path.join('models', general_config['simulation']['model'] + '.ini'))
py_module = 'models.' + model_config['module']['name']
model = __import__(py_module, fromlist = [''])

def delete_random_elems(input_list, n):
    to_delete = set(random.sample(range(len(input_list)), n))
    return [x for i,x in enumerate(input_list) if not i in to_delete]

class Simulation:
    def __init__(
        self,
        population_size,
        generations,
        lifespan,
        emigration,
        inmigration,
        inmigration_genotype,
        model,
        reproduction,
        selection_group_size,
        survival_range,
        altruism_cost_benefit,
        altruism_probability,
        plot_genes):
            self.__population_size = population_size
            self.__generations = generations
            self.__lifespan = lifespan
            self.__emigration = emigration
            self.__inmigration = inmigration
            self.__inmigration_genotype = inmigration_genotype
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
            #First list of the summary is survivors per generation and second the total
            self.__simulation_summary = [[],[]]
            self.__plot_genes = plot_genes
            self.__alleles_combinations_indexes = []

            for gene in model_config.sections():
                if gene != 'module':
                    incorrect_character = re.search('[^\w\s=>]', model_config[gene]['alleles'])
                    if incorrect_character:
                        raise Exception(f'Character "{incorrect_character.group(0)}" in allele not valid')
                    allele_options = re.split('\s*>\s*|\s*=\s*', model_config[gene]['alleles'])
                    noncodominant_allele_options = allele_options.copy()
                    codominant_alleles_raw = re.findall('[\w=]*=[\w=]*', model_config[gene]['alleles'].replace(' ',''))
                    for element in codominant_alleles_raw:
                        for combination in combinations(element.split('='), 2):
                            codominant_phenotype = combination[0]+'_'+combination[1]
                            allele_options.append(codominant_phenotype)
                    for allele in allele_options:
                        if re.search('\s', allele):
                            raise Exception(f'Missing inheritance symbol in "{allele}", allele names can\'t have spaces')
                    initial_frequencies = re.split('\s*,\s*', model_config[gene]['initial_frequencies'])
                    if len(initial_frequencies) != len(noncodominant_allele_options) and len(initial_frequencies) != len(noncodominant_allele_options)-1:
                        raise Exception(f'Incorrect number of allele frequencies given in gene "{gene}", there must be same number as alleles or one less')
                    #First frequence range is set manually as it needn't be calculated
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

            alleles_combinations_indexes = {}
            all_alleles = []
            for gene in self.genes:
                #[gene][0] is equal to gene's alleles
                all_alleles.append(self.__genes_properties[gene][0])
            all_alleles_combinations = list(product(*all_alleles))
            for i in range(len(all_alleles_combinations)):
                allele_combination_name = ''
                for element in all_alleles_combinations[i]:
                    allele_combination_name += element+'&'
                alleles_combinations_indexes[allele_combination_name[:-1]] = i
            self.__alleles_combinations_indexes = alleles_combinations_indexes

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
            #[gene][0] is equal to gene's alleles
            alleles[gene] = self.__genes_properties[gene][0]
        return alleles

    @property
    def alleles_combinations_indexes(self):
        return self.__alleles_combinations_indexes

    def populate(self):
        for i in range(self.__population_size):
            individual = Individual(self)
            individual.genotype = self.generate_individual_genotype()
            self.__population.append(individual)

    def generate_individual_genotype(self):
        allele_pairs = []
        for gene in self.genes:
            alleles = []
            for i in range(2):
                random_picker = random.random()
                #[gene][0] is equal to gene's alleles ranges
                for i in range(len(self.__genes_properties[gene][1])):
                    if random_picker < self.__genes_properties[gene][1][i]:
                        #[gene][0] is equal to gene's alleles
                        alleles.append(self.__genes_properties[gene][0][i])
                        break
            allele_pairs.append(alleles)
        return allele_pairs

    def group_individuals(self):
        population_copy = self.__population.copy()
        groups = []
        if len(population_copy) <= self.__selection_group_size:
            groups.append(population_copy)
            self.__groups = groups
        else:
            while len(population_copy) > self.__selection_group_size:
                group = []
                for i in range(self.__selection_group_size):
                    group.append(population_copy.pop(random.randint(0, len(population_copy) - 1)))
                groups.append(group)
            while len(population_copy) != 0:
                groups[random.randint(0, len(groups) - 1)].append(population_copy.pop(0))
            self.__groups = groups

    def assing_individuals_survival_probability(self):
        for individual in self.__population:
            individual.survival_probability = random.uniform(self.__survival_range[0], self.__survival_range[1])

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
        if self.__lifespan == 1:
            missing_population = self.__population_size
        else:
            survivors_population = [individual for individual in self.__population if individual.age < self.__lifespan]
            missing_population = self.__population_size - len(survivors_population)
        if self.__inmigration != 0:
            inmigrants = round(self.__population_size*self.__inmigration)
            if inmigrants > missing_population:
                raise Exception(f'Inmigration rate cannot exceed {missing_population/self.__population_size} with this configuration')
            for inmigran_index in range(inmigrants):
                individual = Individual(self)
                individual.genotype = self.generate_individual_genotype()
                for allele in self.__inmigration_genotype:
                    for gene_index in range(len(self.genes)):
                        if allele in self.__genes_properties[self.genes[gene_index]][0]:
                            individual.genotype[gene_index] = [allele, allele]
                            break
                new_population.append(individual)
        if self.__reproduction == 0:
            while len(new_population) < missing_population:
                reproductors = random.sample(self.__population, 2)
                new_individual_genotype = []
                #For each gene there's a probability
                for i in range(len(reproductors[0].genotype)):
                    mutation_rate = float(model_config[self.genes[i]]['mutation_rate'])
                    gene_alleles = []
                    for reproductor_index in range(2):
                        reproductor_allele = random.choice(reproductors[reproductor_index].genotype[i])
                        if random.random() < mutation_rate:
                            #[gene][0] is equal to gene's alleles
                            allele_posibilities = self.__genes_properties[self.genes[i]][0].copy()
                            allele_posibilities.remove(reproductor_allele)
                            reproductor_allele = random.choice(allele_posibilities)
                        gene_alleles.append(reproductor_allele)
                    new_individual_genotype.append(gene_alleles)
                new_individual = Individual(self)
                new_individual.genotype = new_individual_genotype
                new_population.append(new_individual)
        else:
            if self.__reproduction > len(self.__population[0].phenotype):
                raise Exception('Number of identical alleles needed for reproduction cannot exceed number of genes')
            count = 0
            while len(new_population) < missing_population:
                possible_reproductors = self.__population.copy()
                first_reproductor = random.choice(possible_reproductors)
                possible_reproductors.remove(first_reproductor)
                same_phenotype_count = 0
                while same_phenotype_count < self.__reproduction:
                    if len(possible_reproductors) == 0:
                        break
                    same_phenotype_count = 0
                    second_reproductor = random.choice(possible_reproductors)
                    for phenotype_index in range(len(first_reproductor.phenotype)):
                        if first_reproductor.phenotype[phenotype_index] == second_reproductor.phenotype[phenotype_index]:
                            same_phenotype_count += 1
                    possible_reproductors.remove(second_reproductor)
                if same_phenotype_count >= self.__reproduction:
                    new_individual_genotype = []
                    reproductors = [first_reproductor, second_reproductor]
                    for i in range(len(reproductors[0].genotype)):
                        mutation_rate = float(model_config[self.genes[i]]['mutation_rate'])
                        gene_alleles = []
                        for reproductor_index in range(2):
                            reproductor_allele = random.choice(reproductors[reproductor_index].genotype[i])
                            if random.random() < mutation_rate:
                                #[gene][0] is equal to gene's alleles
                                allele_posibilities = self.__genes_properties[self.genes[i]][0].copy()
                                allele_posibilities.remove(reproductor_allele)
                                reproductor_allele = random.choice(allele_posibilities)
                            gene_alleles.append(reproductor_allele)
                        new_individual_genotype.append(gene_alleles)
                    new_individual = Individual(self)
                    new_individual.genotype = new_individual_genotype
                    new_population.append(new_individual)
                count += 1
        if self.__lifespan == 1:
            self.__population = new_population
        else:
            for individual in survivors_population:
                individual.age_individual()
            survivors_population.extend(new_population)
            self.__population = survivors_population

    def save_generation_data(self, summary_number):
        #Posible alleles are added in a list, if more than 1 gene is considered
        #all allele combinations will also be added
        combined_phenotypes_list = list(self.__alleles_combinations_indexes.keys())
        if self.current_generation == 0:
            summary = [0 for i in range(self.generations)]
            for i in range(len(combined_phenotypes_list)):
                self.__simulation_summary[0].append(summary.copy())
                self.__simulation_summary[1].append(summary.copy())

        #Code for saving all individuals per generation
        if summary_number == 0:
            for individual in self.__population:
                match_indexes = list(range(len(combined_phenotypes_list)))
                for phenotype in individual.phenotype:
                    new_match_indexes = []
                    for index in match_indexes:
                        combination = combined_phenotypes_list[index]
                        if re.search('(?:&|^)('+phenotype+')(?:&|$)', combination):
                            new_match_indexes.append(index)
                    match_indexes = new_match_indexes
                self.__simulation_summary[1][match_indexes[0]][self.current_generation] += round(1/self.__population_size,6)

        if summary_number == 1:
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
        #start_progress('Progress bar')
        #progress(0)
        self.save_generation_data(0)
        #progress(10)
        #print('Proportions of individuals data saved')
        self.assing_individuals_survival_probability()
        #progress(25)
        #print('Survival probabilities assigned')
        self.group_individuals()
        #progress(40)
        #print('Individuals grouped')
        model.selection(self.groups)
        #progress(55)
        #print('Altruism event finished')
        self.selection_event()
        #progress(70)
        #print('Individuals have been selected')
        self.reproduce()
        #progress(85)
        #print('Individuales have reproduced')
        self.save_generation_data(1)
        #progress(100)
        #print('Survivors data saved')
        #print(self.current_generation)
        self.__generation += 1

class Individual:
    def __init__(self, simulation):
        self.__simulation = simulation
        self.__age = 0
        self.__genotype = []
        self.__genes_properties = {}
        self.__phenotype = []
        self.__survival_probability = 0

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

    #Calculates the phenotype considering the alleles and inheritance pattern
    def choose_phenotype(self):
        phenotype = []
        for i in range(len(self.__simulation.genes)):
            gene = self.__simulation.genes[i]
            #If both allele are the same the phenotype is that allele
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
                    #[gene][0] is equal to gene's alleles
                    if chosen_phenotype in self.__simulation.genes_properties[gene][0]:
                        phenotype.append(chosen_phenotype)
                    else:
                        phenotype.append(self.__genotype[i][1]+'_'+self.__genotype[i][0])
        self.phenotype = phenotype

    def age_individual(self):
        self.__age += 1

def alts_main():
    generations = int(general_config['simulation']['generations'])+1
    population_size = int(general_config['population']['size'])
    lifespan = int(general_config['population']['lifespan'])
    emigration = float(general_config['population']['emigration'])
    inmigration = float(general_config['population']['inmigration'])
    inmigration_genotype = general_config['population']['inmigration_genotype'].replace(' ','').split(',')
    reproduction = int(general_config['population']['reproduction'])
    selection_group_size = int(general_config['population']['selection_group_size'])
    survival_range = [float(perc) for perc in re.split(',\s*', general_config['population']['survival_range'])]
    altruism_cost_benefit = [float(general_config['population']['altruism_cost']), float(general_config['population']['altruism_benefit'])]
    altruism_probability = float(general_config['population']['altruism_probability'])
    model = general_config['simulation']['model']
    plot_genes = general_config['plot']['genes'].replace(' ','').split(',')

    p = Simulation(population_size, generations, lifespan, emigration, inmigration, inmigration_genotype, model, reproduction, selection_group_size, survival_range, altruism_cost_benefit, altruism_probability, plot_genes)

    p.populate()
    for i in range(generations):
        p.pass_generation()
    return p.simulation_summary, p.alleles_combinations_indexes, p.dict_allele_options

if __name__ == '__main__':
    main()
