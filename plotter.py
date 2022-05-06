import re
import pickle
import argparse

from configparser import ConfigParser
from matplotlib import pyplot as plt
#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

parser = argparse.ArgumentParser(description='Altruism simulations')
parser.add_argument('-i', '--input', default = 'result.txt', help = 'input file where data will be extracted')
args = parser.parse_args()

def plot():
    population_size = int(general_config['population']['size'])
    generation_x = range(int(general_config['simulation']['generations'])+1)

    with open(vars(args)['input'], 'rb') as config_dictionary_file:
        simulations_summary = pickle.load(config_dictionary_file)
    #print(simulations_summary)
    dict_phenotypes_combinations_indexes = simulations_summary[0]
    dict_phenotype_options = simulations_summary[1]
    survivors_means = simulations_summary[4]
    proportions_means = simulations_summary[5]
    all_simulations_summary = simulations_summary[6]

    legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())

    individual_survivors_plot = plt.figure(3)
    for simulation in all_simulations_summary:
        plt.plot(simulation)
    plt.margins(0)
    plt.ylim(0,1)
    plt.title('Selfish individuals per simulation')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')

    legend_population_size = 'N = '+str(population_size)
    legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    proportions_plot = plt.figure(1)
    plt.stackplot(generation_x, proportions_means, labels = legend_phenotypes)
    plt.margins(0)
    plt.ylim(0,1)
    plt.title('Proportion of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')
    plt.legend()

    total_survivors = []
    for generation in generation_x:
        individuals_sum = 0
        for phenotype_index in range(len(legend_phenotypes)):
            individuals_sum += survivors_means[phenotype_index][generation]
        total_survivors.append(individuals_sum)

    survivors_means.append(total_survivors)

    survivors_plot = plt.figure(2)
    legend_phenotypes.append('total')
    for count in range(len(survivors_means)):
        plt.plot(survivors_means[count], label = legend_phenotypes[count])

    plt.margins(0)
    plt.ylim((0,population_size))
    plt.title('Survivors per generation by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Number of individuals')
    plt.legend()

    plt.show()

plot()
