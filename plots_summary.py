
import re
from itertools import product
from time import perf_counter

from configparser import ConfigParser
from matplotlib import pyplot as plt

from alts import alts_main

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')
with open('config2.ini', 'w') as f:
    general_config.write(f, False)

start_time = perf_counter()

def create_simulation_results():
    start_counter = perf_counter()
    number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    population_size = int(general_config['population']['size'])
    generation_x = range(int(general_config['simulation']['generations'])+1)
    survivors_simulations_summary = []
    total_simulations_summary = []
    simulation_count = 0
    individual_survivors_plot = plt.figure(3)

    simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = alts_main()
    print(f'\rsimulation number {simulation_count} run in {round(perf_counter() - start_counter,2)} seconds', end = '')
    #sys.stdout.flush()
    mean_simulation_duration = perf_counter() - start_counter
    start_counter = perf_counter()
    simulation_count += 1

    legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    survivors_simulations_summary.append(simulation_summary[0])
    total_simulations_summary.append(simulation_summary[1])
    altruism_indexes = []

    for phentoype_index in range(len(legend_phenotypes)):
        if re.search('(?:&|^)(selfish)(?:&|$)', legend_phenotypes[phentoype_index]):
            altruism_indexes.append(phentoype_index)

    for phenotype in altruism_indexes:
        plt.plot(simulation_summary[1][phenotype], label = legend_phenotypes[phenotype])

    #The given number of simulations are run and their data is saved
    for i in range(number_of_simulations):
        simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = alts_main()
        survivors_simulations_summary.append(simulation_summary[0])
        total_simulations_summary.append(simulation_summary[1])
        single_simulation_summary = [0 for generation in generation_x]
        for phenotype in altruism_indexes:
            single_simulation_summary = [sum(x) for x in zip(single_simulation_summary, simulation_summary[1][phenotype])]
        plt.plot(single_simulation_summary)
        print(f'\rsimulation number {simulation_count} run in {round(perf_counter() - start_counter,2)} seconds', end = '')
        #sys.stdout.flush()
        mean_simulation_duration += perf_counter() - start_counter
        simulation_count += 1
        start_counter = perf_counter()

    mean_simulation_duration = mean_simulation_duration/number_of_simulations
    print(f'\rEach simulation took {round(mean_simulation_duration,2)} seconds on average and {round((perf_counter() - start_time)/60,2)} minutes in total')
    plt.margins(0)
    plt.ylim(0,1)
    plt.title('Selfish individuals per simulation')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')
    plt.legend()

    combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    survivors_data_y = []
    proportions_data_y = []
    for phenotype in range(len(dict_phenotypes_combinations_indexes.keys())):
        survivors_means_per_phenotype = []
        total_means_per_phenotype = []
        for generation in generation_x:
            survivors_mean = 0
            total_mean = 0
            for simulation in range(number_of_simulations):
                survivors_mean += survivors_simulations_summary[simulation][phenotype][generation]
                total_mean += total_simulations_summary[simulation][phenotype][generation]
            survivors_means_per_phenotype.append(survivors_mean/number_of_simulations)
            total_means_per_phenotype.append(total_mean/number_of_simulations)
        survivors_data_y.append(survivors_means_per_phenotype)
        proportions_data_y.append(total_means_per_phenotype)

    total_survivors = []
    for generation in generation_x:
        individuals_sum = 0
        for phenotype_index in range(len(combined_phenotypes)):
            individuals_sum += survivors_data_y[phenotype_index][generation]
        total_survivors.append(individuals_sum)

    legend_population_size = 'N = '+str(population_size)
    legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    survivors_data_y.append(total_survivors)
    proportions_plot = plt.figure(1)
    plt.stackplot(generation_x, proportions_data_y, labels = legend_phenotypes)
    plt.margins(0)
    plt.ylim(0,1)
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

create_simulation_results()
