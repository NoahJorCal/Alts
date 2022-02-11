#!/usr/bin/python3

from os import path
from configparser import ConfigParser
import re
import random
#import alts_plot
from collections import defaultdict
from itertools import combinations
from itertools import product

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
    def __init__(
        self,
        population_size,
        generations,
        lifespan, model,
        reproduction,
        selection_group_size,
        survival_range,
        altruism_cost_benefit,
        altruism_probability,
        plot_genes):
            self.__population_size = population_size
            self.__generations = generations
            self.__last_population_size = 0
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
            self.__simulation_summary = []
            self.__plot_genes = plot_genes
            self.__alleles_combinations_indexes = []

            for gene in model_config.sections():
                if gene != 'module':
                    incorrect_character = re.search('[^\w\s=>]', model_config[gene]['alleles'])
                    if incorrect_character:
                        raise Exception(f'Character "{incorrect_character.group(0)}" in allele not valid')
                    allele_options = re.split('\s*>\s*|\s*=\s*', model_config[gene]['alleles'])
                    codominant_alleles_raw = re.findall('[\w=]*=[\w=]*', model_config[gene]['alleles'].replace(' ',''))
                    #print(codominant_alleles_raw)
                    for element in codominant_alleles_raw:
                        #print(element.split('='))
                        for combination in combinations(element.split('='), 2):
                            codominant_phenotype = combination[0]+'_'+combination[1]
                            allele_options.append(codominant_phenotype)
                            #print(combination)

                            #print(codominant_phenotype)
                        #list())
                    #print(allele_options)
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

            alleles_combinations_indexes = {}
            all_alleles = []
            for gene in self.genes:
                all_alleles.append(self.__genes_properties[gene][ALLELES])
            all_alleles_combinations = list(product(*all_alleles))
            for i in range(len(all_alleles_combinations)):
                allele_combination_name = ''
                for element in all_alleles_combinations[i]:
                    allele_combination_name += element+'&'
                alleles_combinations_indexes[allele_combination_name[:-1]] = i
            self.__alleles_combinations_indexes = alleles_combinations_indexes

            #config_genes = model_config.sections()
            #config_genes.remove('module')
            ##print(config_genes)
            #genes_combinations = []
            #for i in range(2, len(config_genes)+1):
                #genes_combinations += combinations(config_genes, i)
            ##print(genes_combinations)
            #for combination in genes_combinations:
                #combination_name = ''
                #for element in combination:
                    #combination_name += element + '&'
                ##print(combination_name[0:-2])
                #config_genes.append(combination_name[0:-2])
            #print(config_genes)
            #for gene in config_genes:
                #print(gene)
                #incorrect_character = re.search('[^\w\s=>]', model_config[gene]['alleles'])
                #if incorrect_character:
                    #raise Exception(f'Character "{incorrect_character.group(0)}" in allele not valid')
                #allele_options = re.split('\s*>\s*|\s*=\s*', model_config[gene]['alleles'])
                #codominant_alleles_raw = re.findall('[\w=]*=[\w=]*', model_config[gene]['alleles'].replace(' ',''))
                ##print(codominant_alleles_raw)


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
    def survival_range(self):
        return self.__survival_range

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
            alleles[gene] = self.__genes_properties[gene][ALLELES]
        return alleles

    @property
    def alleles_combinations_indexes(self):
        return self.__alleles_combinations_indexes

    def populate(self):
        for i in range(self.__population_size):
            individual = Individual(self)
            self.__population.append(individual)

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
        #if len(self.__population) < 2:
            #alts_plot.plot(self.simulation_summary)
            #exit()
        #else:
        if self.__reproduction == 'panmixia':
            new_population = []
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



    def save_generation_data(self):
        #gene = model_config['plot']['gene']
        #allele_options = re.split('\s*>\s*|\s*=\s*', model_config[gene]['alleles'])
        #summary_dictionary = {}
        #for allele in allele_options:
            #summary_dictionary[allele] = []
        #return summary_dictionary
        ##The name of the gene that will be plotted is saved
        #gene = model_config['plot']['gene']
        #config_sections = model_config.sections()
        #config_sections.remove('module')
        #config_sections.remove('plot')
        #gene_index = config_sections.index(gene)
        combined_alleles_list = list(self.__alleles_combinations_indexes.keys())
        if self.current_generation == 0:
            summary = []
            for i in range(len(combined_alleles_list)):
                summary.append([0]*self.generations)
            self.__simulation_summary = summary
        for individual in self.__population:
            match_indexes = list(range(len(combined_alleles_list)))
            for phenotype in individual.phenotype:
                new_match_indexes = []
                for index in match_indexes:
                    combination = combined_alleles_list[index]
                    if re.search('(?:&|^)('+phenotype+')(?:&|$)', combination):
                        new_match_indexes.append(index)
                match_indexes = new_match_indexes
            self.__simulation_summary[match_indexes[0]][self.current_generation] += 1
            #combined_phenotype = list(self.simulation_summary['combinations'].keys())[match_indexes[0]]
            #self.simulation_summary['combinations'][combined_phenotype][self.current_generation] += 1


            #for gene in self.genes:
                #genes = []
                #for i in range(len(self.__genes_properties[gene][ALLELES])):
                    #genes.append([0]*self.generations)
                #summary.append(genes)
            #self.simulation_summary = summary
            #print('==============')
            #print(self.simulation_summary)
            #print('==============')
            #if len(self.__plot_genes) > 1:
                #all_alleles = []
                #for gene in self.genes:
                    #all_alleles.append(self.__genes_properties[gene][ALLELES])
                #for combination in product(*all_alleles):
                    #combinations_summary[combination[0]+'__'+combination[1]] = [0]*self.generations
                #self.simulation_summary['combinations'] = combinations_summary
        #for individual in self.__population:
            #for i in range(len(self.genes)):
                ##Each individual is added by their phenotype
                #self.simulation_summary[self.genes[i]][individual.phenotype[i]][self.current_generation] += 1
            #if len(self.__plot_genes) > 1:
                #match_indexes = list(range(len(self.simulation_summary['combinations'].keys())))
                #for phenotype in individual.phenotype:
                    #new_match_indexes = []
                    #for index in match_indexes:
                        #combination = list(self.simulation_summary['combinations'].keys())[index]
                        #if re.search('(?:__|^)('+phenotype+')(?:__|$)', combination):
                            #new_match_indexes.append(index)
                    #match_indexes = new_match_indexes
                #combined_phenotype = list(self.simulation_summary['combinations'].keys())[match_indexes[0]]
                #self.simulation_summary['combinations'][combined_phenotype][self.current_generation] += 1



                #for phenotype in individual.phenotype:
                #for gene in self.__plot_genes:

                #combinations_summary = {}
                #all_alleles = []
                #for gene in self.genes:
                    #all_alleles.append(self.__genes_properties[gene][ALLELES])
                #for combination in product(*all_alleles):+
                    #combinations_summary[combination[0]+'_'+combination[1]] = [0]*self.generations
                #for individual in self.__population:
                    #for i in range(len(self.genes)):
                        ##Each individual is added by their phenotype
                        #self.simulation_summary[self.genes[i]][individual.phenotype[i]][self.current_generation] += 1
        #print(self.current_generation)
        #print(self.simulation_summary)


    #def save_generation_data(self):
        ##The name of the gene that will be plotted is saved
        #gene = model_config['plot']['gene']
        #config_sections = model_config.sections()
        #config_sections.remove('module')
        #config_sections.remove('plot')
        #gene_index = config_sections.index(gene)
        ##A dictionary with lists as values is created
        #groups = {}
        ##print('======================')
        ##print(self.__genes_properties[gene][ALLELES])
        #for phenotype in self.__genes_properties[gene][ALLELES]:
            #print(phenotype)
            #groups[phenotype] = 0
            #print(groups)
        #for individual in self.__population:
            ##Each individual is added by their phenotype
            #groups[individual.phenotype[gene_index]] += 1
        ##print(self.__simulation_summary)
        #for phenotype in groups.keys():
            #self.__simulation_summary[phenotype].append(groups[phenotype])
        #pass
        #for phenotype in groups.keys():
            #if groups[phenotype] == len(self.__population):
                #print('000000000000000000000')
                #print(self.__simulation_summary.keys())
                #print('000000000000000000000')

                #for phenotype_drift in self.__simulation_summary.keys():
                    #self.__simulation_summary[phenotype_drift].append(0)
                #self.__simulation_summary[phenotype].pop()
                #self.__simulation_summary[phenotype].append(len(groups[phenotype]))
            #else:
                #self.__simulation_summary[phenotype].append(len(groups[phenotype]))


    def pass_generation(self):
        #print('EMPIEZA UNA GENERACION')
        #print(len(self.__population))
        self.save_generation_data()
        #for individual in self.__population:
            #print(individual.phenotype)
        #print(self.simulation_summary)
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
        #print('HA PASADO UNA GENERACION')
        #print(self.__generation)


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
            #print(self.__genotype)
            if self.__genotype[i][0] == self.__genotype[i][1]:
                phenotype.append(self.__genotype[i][0])
            else:
                characters_between_alleles = re.search(
                    f'{INHERITANCE_DIVIDER}{self.__genotype[i][0]}{INHERITANCE_DIVIDER}(.*){INHERITANCE_DIVIDER}{self.__genotype[i][1]}{INHERITANCE_DIVIDER}',
                    model_config[gene]['alleles'])
                reverse = False
                if not characters_between_alleles:
                    #print('BASHJDBAOSIHDSA')
                    characters_between_alleles = re.search(
                        f'{INHERITANCE_DIVIDER}{self.__genotype[i][1]}{INHERITANCE_DIVIDER}(.*){INHERITANCE_DIVIDER}{self.__genotype[i][0]}{INHERITANCE_DIVIDER}',
                        model_config[gene]['alleles'])
                    reverse = True
                #print(characters_between_alleles)
                if '>' in characters_between_alleles.group(1):
                    if reverse:
                        phenotype.append(self.__genotype[i][1])
                    else:
                        phenotype.append(self.__genotype[i][0])
                else:
                    chosen_phenotype = self.__genotype[i][0]+'_'+self.__genotype[i][1]
                    if chosen_phenotype in self.__simulation.genes_properties[gene][ALLELES]:
                        phenotype.append(chosen_phenotype)
                    else:
                        phenotype.append(self.__genotype[i][1]+'_'+self.__genotype[i][0])
        return phenotype

    def age_individual(self):
        self.__age += 1

def alts_main():
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
    plot_genes = general_config['plot']['genes'].replace(' ','').split(',')
    p = Simulation(population_size, generations, lifespan, model, reproduction, selection_group_size, survival_range, altruism_cost_benefit, altruism_probability, plot_genes)

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


    #alts_plot.plot(p.simulation_summary, p.population_size)
    return p.simulation_summary, p.alleles_combinations_indexes, p.dict_allele_options

if __name__ == '__main__':
    main()
