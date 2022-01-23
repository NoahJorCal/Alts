#!/usr/bin/python3

from os import path
from configparser import ConfigParser
import re
import random

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

#Import model configuration
model_name = general_config['simulation']['model']
model_config = ConfigParser()
model_config.read(path.join('models', model_name + '.ini'))
py_module = 'models.' + model_config['module']['name']
model = __import__(py_module, fromlist=[''])

class Simulation:
    def __init__(self, population_size, lifespan, model):
        self.__population_size = population_size
        self.__lifespan = lifespan
        self.__generation = 0
        self.__population = []
        self.__model = model

    #@property
    #def population_size(self):
        #return self.__population_size

    #@property
    #def current_generation(self):
        #return self.__generation

    def populate(self):
        for i in range(self.__population_size):
            self.__population.append(Individual())

    @property
    def population(self):
        return self.__population

    def pass_generation(self):
        #Do stuff
        self.__generation += 1

class Individual():
    def __init__(self):
        self.__age = 0
        self.__genotype = self.set_individual_genotype()
        self.__phenotype = self.choose_phenotype()

    #Returns a list of random alleles of each gene considering initial frequencies
    def set_individual_genotype(self):
        allele_pairs = []
        for gene in model_config.sections()[1:]:
            alleles = []
            for i in range(2):
                #print(gene)
                allele_options = re.split('\s*>\s*|\s*<\s*|\s*=\s*',model_config[gene]['alleles'])
                initial_frequencies = re.split(',\s*',model_config[gene]['initial_frequencies'])
                #First frequence range is set manually as it needn't be calculated
                initial_frequencies_ranges = [float(initial_frequencies[0])]
                for i in range(len(initial_frequencies))[1:]:
                    initial_frequencies_ranges.append(float(initial_frequencies[i])+initial_frequencies_ranges[i-1])
                if len(initial_frequencies_ranges) != len(allele_options):
                    initial_frequencies_ranges.append(1)
                random_picker = random.random()
                for i in range(len(initial_frequencies_ranges)):
                    if random_picker < initial_frequencies_ranges[i]:
                        alleles.append(allele_options[i])
                        break
            allele_pairs.append(alleles)
        return allele_pairs

    #Calculates the phenotype considering the alleles and inheritance pattern
    def choose_phenotype(self):
        #print(self.genotype)
        genotype = []
        for i in range(len(self.__genotype)):
            gene = model_config.sections()[i+1]
            #print(gene)
            #print(model_config[gene]['alleles'])
            if self.__genotype[i][0] == self.__genotype[i][1]:
                genotype.append(self.__genotype[i][0])
            else:
                between_alleles = re.search(self.__genotype[i][0]+'(.*)'+self.__genotype[i][1], model_config[gene]['alleles'])
                reverse = False
                if not between_alleles:
                    between_alleles = re.search(self.__genotype[i][1]+'(.*)'+self.__genotype[i][0], model_config[gene]['alleles'])
                    reverse = True
                raw_symbol = re.search('>|<',between_alleles.group(1)[::-1])
                if not raw_symbol:
                    config_symbol = '='
                else:
                    config_symbol = raw_symbol.group()
                if config_symbol == '>':
                    if reverse:
                        genotype.append(self.__genotype[i][1])
                    else:
                        genotype.append(self.__genotype[i][0])
                elif config_symbol == '<':
                    if reverse:
                        genotype.append(self.__genotype[i][0])
                    else:
                        genotype.append(self.__genotype[i][1])
                else:
                        genotype.append(self.__genotype[i][0]+'_'+self.__genotype[i][1])
        return genotype

    @property
    def genotype(self):
        return self.__genotype

    @property
    def phenotype(self):
        return self.__phenotype

    #@property
    #def age(self):
        #return self.__age

ind = Individual()
print(ind.genotype)
print(ind.phenotype)

def main():
    population_size = int(general_config['population']['size'])
    lifespan = int(general_config['population']['lifespan'])
    model = general_config['simulation']['model']
    p = Simulation(population_size, lifespan, model)
    #p.populate()
    #for i in range(len(p.population)):
        #print(p.population[i].genotype)
    p.pass_generation()
    p.pass_generation()
    p.pass_generation()
    p.pass_generation()
    p.pass_generation()
    p.pass_generation()

if __name__ == '__main__':
    main()
