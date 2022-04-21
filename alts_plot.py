
import re
from itertools import product

from configparser import ConfigParser
from matplotlib import pyplot as plt

from alts import alts_main

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

def plot():
    #Runs the main function of the alts script and generates a summary
    simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = alts_main()
    population_size = int(general_config['population']['size'])
    genes = general_config['plot']['genes'].replace(' ','').split(',')
    generation_x = list(range(int(general_config['simulation']['generations'])+1))
    combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    plot_info = {}
    phenotype_options = []
    for gene in genes:
        phenotype_options.append(dict_phenotype_options[gene])

    if len(genes) == 1:
        for phenotype_index in range(len(phenotype_options[0])):
            survivors_summary = ([0] * (generation_x[-1] + 1))
            total_summary = ([0] * (generation_x[-1] + 1))
            for combination_index in range(len(combined_phenotypes)):
                if re.search('(?:&|^)('+phenotype_options[0][phenotype_index]+')(?:&|$)', combined_phenotypes[combination_index]):
                    for generation in generation_x:
                        survivors_summary[generation] += simulation_summary[0][combination_index][generation]
                        total_summary[generation] += (simulation_summary[1][combination_index][generation])
                        plot_info[phenotype_options[0][phenotype_index]] = []
                        plot_info[phenotype_options[0][phenotype_index]].append(survivors_summary)
                        plot_info[phenotype_options[0][phenotype_index]].append(total_summary)

    else:
        for combination in product(*phenotype_options):
            match_indexes = list(range(len(simulation_summary[0])))
            phenotype_name = ''
            for element in combination:
                new_match_indexes = []
                for index in match_indexes:
                    phenotype = list(dict_phenotypes_combinations_indexes.keys())[index]
                    if re.search('(?:&|^)('+element+')(?:&|$)', phenotype):
                        new_match_indexes.append(index)
                match_indexes = new_match_indexes
                phenotype_name += element+'&'

            survivors_summary = [0] * (generation_x[-1] + 1)
            total_summary = [0] * (generation_x[-1] + 1)
            for generation in generation_x:
                for index in match_indexes:
                    survivors_summary[generation] += simulation_summary[0][index][generation]
                    total_summary[generation] += (simulation_summary[1][index][generation])/population_size
            plot_info[phenotype_name[:-1]] = []
            plot_info[phenotype_name[:-1]].append(survivors_summary)
            plot_info[phenotype_name[:-1]].append(total_summary)

    total_survivors = []
    for generation in generation_x:
        individuals_sum = 0
        for phenotype_index in range(len(plot_info)):
            individuals_sum += plot_info[list(plot_info.keys())[phenotype_index]][0][generation]
        total_survivors.append(individuals_sum)

    legend_population_size = 'N = '+str(population_size)
    proportions_data_y = [values[1] for values in list(plot_info.values())]
    survivors_data_y = [values[0] for values in list(plot_info.values())]
    survivors_data_y.append(total_survivors)
    legend_phenotypes = list(plot_info.keys())

    proportions_plot = plt.figure(1)
    plt.stackplot(generation_x, proportions_data_y, labels = legend_phenotypes)
    plt.margins(0)
    plt.title('Proportion of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')
    plt.legend()

    survivors_plot = plt.figure(2)
    legend_phenotypes.append('total')
    for count in range(len(survivors_data_y)):
        plt.plot(survivors_data_y[count], label = legend_phenotypes[count])
    plt.margins(0)
    plt.ylim((0,population_size))
    plt.title('Survivors per generation by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Number of individuals')
    plt.legend()

    plt.show()

plot()
