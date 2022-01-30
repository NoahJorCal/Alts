#!/usr/bin/python3

from os import path
from configparser import ConfigParser
import re
import random
import alts_plot
from collections import defaultdict

ALLELES = 0
RANGES = 1
INHERITANCE_DIVIDER = '(?:>|=|\s|$|^)'

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

#Import model configuration
model_name = general_config['simulation']['model']
model_config = ConfigParser()
model_config.read(path.join('models', model_name + '.ini'))
py_module = 'models.' + model_config['module']['name']
model = __import__(py_module, fromlist = [''])

class Simulation:
    def __init__(self, population_size, lifespan, model, reproduction, selection_group_size, survival_range, altruism_cost_benefit, altruism_probability):
        self.__population_size = population_size
        self.__last_population_size = 0
        #self.__size_replenish = size_replenish
        self.__lifespan = lifespan
        self.__generation = 0
        self.__population = []
        self.__model = model
        self.__genes_properties = {}
        self.__reproduction = reproduction
        self.__selection_group_size = selection_group_size
        self.__groups = []
        self.__survival_range = survival_range
        self.__min_max_survival_probability = []
        self.__survival_rate_mean = 0
        self.__altruism_cost_benefit = altruism_cost_benefit
        self.__altruism_probability = altruism_probability
        self.__simulation_summary = self.generate_simulation_summary()

        for gene in model_config.sections():
            if gene != 'module' and gene != 'plot':
                incorrect_character = re.search('[^\w\s=>]', model_config[gene]['alleles'])
                if incorrect_character:
                    raise Exception(f'Character "{incorrect_character.group(0)}" in allele not valid')
                allele_options = re.split('\s*>\s*|\s*=\s*', model_config[gene]['alleles'])
                for allele in allele_options:
                    if re.search('\s', allele):
                        raise Exception(f'Missing inheritance symbol in "{allele}", allele names can\'t have spaces')
                initial_frequencies = re.split('\s*,\s*', model_config[gene]['initial_frequencies'])
                #First frequence range is set manually as it needn't be calculated
                initial_frequencies_ranges = [float(initial_frequencies[0])]
                for i in range(len(initial_frequencies))[1:]:
                    initial_frequencies_ranges.append(float(initial_frequencies[i])+initial_frequencies_ranges[i-1])
                if len(initial_frequencies_ranges) != len(allele_options):
                    initial_frequencies_ranges.append(1)
                self.__genes_properties[gene] = (allele_options, initial_frequencies_ranges)

    #@property
    #def population_size(self):
        #return self.__population_size

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
        return sorted(self.__genes_properties.keys())

    @property
    def survival_range(self):
        return self.__survival_range

    @property
    def groups(self):
        return self.__groups

    @property
    def simulation_summary(self):
        return self.__simulation_summary

    def populate(self):
        for i in range(self.__population_size):
            individual = Individual(self)
            self.__population.append(individual)
        self.save_generation_data()

    def generate_individual_genotype(self):
        allele_pairs = []
        for gene in self.genes:
            alleles = []
            for i in range(2):
                random_picker = random.random()
                for i in range(len(self.__genes_properties[gene][RANGES])):
                    if random_picker < self.__genes_properties[gene][RANGES][i]:
                        alleles.append(self.__genes_properties[gene][ALLELES][i])
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
        survivors_population = []
        survival_probabilities = []
        for individual in self.__population:
            survival_probabilities.append(individual.survival_probability)
            picker = random.random()
            if picker < individual.survival_probability:
                survivors_population.append(individual)
        self.population = survivors_population
        self.__survival_rate_mean = sum(survival_probabilities)/len(survival_probabilities)

    #def selection_event(self):
        #survivors_population = []
        ##print('==================================================')
        ##print(self.__min_max_survival_probability)
        ##print(self.__min_max_survival_probability[0], min(1,self.__min_max_survival_probability[1]))
        #random_picker = random.uniform(self.__min_max_survival_probability[0], min(1,self.__min_max_survival_probability[1]))
        #print(random_picker)
        ##print('==================================================')
        #for individual in self.__population:
            #if random_picker < individual.survival_probability:
                #survivors_population.append(individual)
        #self.population = survivors_population

    def reproduce(self):
        #print(len(self.__population))
        if len(self.__population) < 2:
            alts_plot.plot(self.simulation_summary)
            exit()
        else:
            if self.__reproduction == 'panmixia':
                new_population = []
                #if self.__size_replenish == 'False':
                    #minimum_survival_probability = self.survival_range[0] - self.__altruism_cost_benefit[0]
                    #maximum_survival_probability = min(1, (self.survival_range[1] + self.__altruism_cost_benefit[1]))
                    #multiplication_factor = (((self.__survival_rate_mean - minimum_survival_probability) * (1.05 - 0.95)) / (maximum_survival_probability - minimum_survival_probability)) + 0.95
                    #print(minimum_survival_probability, maximum_survival_probability)
                    #if self.__lifespan == 1:
                        ##print(round(len(self.__population)*
                        ##((self.__last_population_size+(self.__last_population_size - len(self.__population)))/self.__last_population_size)))
                        #print(self.__last_population_size, self.__survival_rate_mean)
                        #missing_population = self.__last_population_size * multiplication_factor
                        #print(missing_population)
                        ##random.randrange(round(len(self.__population)), round(len(self.__population) + (self.__population_size * 0.2)))
                        ##*((self.__last_population_size+(self.__last_population_size - len(self.__population)))/self.__last_population_size)))
                    #else:
                        #missing_population = random.randrange(len(self.__population), round(len(self.__population)*1.5))
                        #survivors_population = [individual for individual in self.__population if individual.age < self.__lifespan]
                #elif self.__size_replenish == 'True':
                    #if self.__lifespan == 1:
                        #missing_population = self.__population_size
                    #else:
                        #survivors_population = [individual for individual in self.__population if individual.age < self.__lifespan]
                        #missing_population = self.__population_size - len(survivors_population)
                if self.__lifespan == 1:
                    missing_population = self.__population_size
                else:
                    survivors_population = [individual for individual in self.__population if individual.age < self.__lifespan]
                    missing_population = self.__population_size - len(survivors_population)
                while len(new_population) < missing_population:
                    reproductors = random.sample(self.__population, 2)
                    new_individual_genotype = []
                    for i in range(len(reproductors[0].genotype)):
                        new_individual_genotype.append([random.choice(reproductors[0].genotype[i]), random.choice(reproductors[1].genotype[i])])
                    new_individual = Individual(self)
                    new_individual.genotype = new_individual_genotype
                    new_individual.phenotype = new_individual.choose_phenotype()
                    new_population.append(new_individual)
                if self.__lifespan == 1:
                    self.__population = new_population
                else:
                    for individual in survivors_population:
                        individual.age_individual()
                    survivors_population.extend(new_population)
                    self.__population = survivors_population

    def generate_simulation_summary(self):
        gene = model_config['plot']['gene']
        allele_options = re.split('\s*>\s*|\s*=\s*', model_config[gene]['alleles'])
        summary_dictionary = {}
        for allele in allele_options:
            summary_dictionary[allele] = []
        return summary_dictionary

    def save_generation_data(self):
        gene = model_config['plot']['gene']
        config_sections = model_config.sections()
        config_sections.remove('module')
        config_sections.remove('plot')
        gene_index = config_sections.index(gene)
        groups = defaultdict(list)
        for individual in self.__population:
            groups[individual.phenotype[gene_index]].append(individual)
        for phenotype in groups.keys():
            if len(groups[phenotype]) == len(self.__population):
                for phenotype_drift in self.__simulation_summary.keys():
                    self.__simulation_summary[phenotype_drift].append(0)
                self.__simulation_summary[phenotype].pop()
                self.__simulation_summary[phenotype].append(len(groups[phenotype]))
            else:
                self.__simulation_summary[phenotype].append(len(groups[phenotype]))


    def pass_generation(self):
        #print('EMPIEZA UNA GENERACION')
        #print(len(self.__population))
        self.assing_individuals_survival_probability()
        #print('SE HA ASIGNADO LA PROBABILIDAD DE SOBREVIVIR')
        self.group_individuals()
        #print('SE HAN AGRUPADO LOS INDIVIDUOS')
        self.__last_population_size = len(self.__population)
        self.__min_max_survival_probability = model.selection(self.groups)
        #print('SE HA HECHO AL MOVIDA DEL ALTRUISMO')
        self.selection_event()
        #print('SE HAN MUERTO LOS QUE TENÍAN POCA PROBABILIDAD')
        self.reproduce()
        #print('SE HAN REPRODUCIDO')
        self.__generation += 1
        self.save_generation_data()
        #print('HA PASADO UNA GENERACION')
        print(self.__generation)


class Individual:
    def __init__(self, simulation):
        self.__simulation = simulation
        self.__age = 0
        self.__genotype = simulation.generate_individual_genotype()
        self.__genes_properties = {}
        self.__phenotype = self.choose_phenotype()
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
            if self.__genotype[i][0] == self.__genotype[i][1]:
                phenotype.append(self.__genotype[i][0])
            else:
                characters_between_alleles = re.search(
                    f'{INHERITANCE_DIVIDER}{self.__genotype[i][0]}{INHERITANCE_DIVIDER}(.*){INHERITANCE_DIVIDER}{self.__genotype[i][1]}{INHERITANCE_DIVIDER}',
                    model_config[gene]['alleles'])
                reverse = False
                if not characters_between_alleles:
                    characters_between_alleles = re.search(
                        f'{INHERITANCE_DIVIDER}{self.__genotype[i][1]}{INHERITANCE_DIVIDER}(.*){INHERITANCE_DIVIDER}{self.__genotype[i][0]}{INHERITANCE_DIVIDER}',
                        model_config[gene]['alleles'])
                    reverse = True
                if '>' in characters_between_alleles.group(1):
                    if reverse:
                        phenotype.append(self.__genotype[i][1])
                    else:
                        phenotype.append(self.__genotype[i][0])
                else:
                    phenotype.append(self.__genotype[i][0]+'_'+self.__genotype[i][1])
        return phenotype

    def age_individual(self):
        self.__age += 1

def main():
    generations = int(general_config['simulation']['generations'])
    population_size = int(general_config['population']['size'])
    #size_replenish = general_config['population']['size_replenish']
    lifespan = int(general_config['population']['lifespan'])
    reproduction = general_config['population']['reproduction']
    selection_group_size = int(general_config['population']['selection_group_size'])
    survival_range = [float(perc) for perc in re.split(',\s*', general_config['population']['survival_range'])]
    altruism_cost_benefit = [float(general_config['population']['altruism_cost']), float(general_config['population']['altruism_benefit'])]
    altruism_probability = float(general_config['population']['altruism_probability'])
    model = general_config['simulation']['model']
    p = Simulation(population_size, lifespan, model, reproduction, selection_group_size, survival_range, altruism_cost_benefit, altruism_probability)

    #ind = Individual(p)
    #print(ind.genotype)
    #print(ind.phenotype)

    p.populate()
    for i in range(generations):
        p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #print()
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #print()
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #print()
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #print()
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #p.pass_generation()
    #for i in range(len(p.population)):
        #print(p.population[i], p.population[i].phenotype[0], p.population[i].age, p.population[i].genotype[0])
    #print()
    #p.pass_generation()
    #p.pass_generation()
    #p.pass_generation()
    #print('despues de la generacion', p.population)
    alts_plot.plot(p.simulation_summary)

if __name__ == '__main__':
    main()
