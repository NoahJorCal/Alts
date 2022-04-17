from configparser import ConfigParser
from matplotlib import pyplot as plt
from alts import alts_main
import re
from itertools import product

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')
with open('config2.ini', 'w') as f:
    general_config.write(f, False)


def create_simulation_results():
    number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    population_size = int(general_config['population']['size'])
    generation_x = range(int(general_config['simulation']['generations'])+1)
    survivors_simulations_summary = []
    total_simulations_summary = []
    #The given number of simulations are run and their data is saved
    count = 0
    for i in range(number_of_simulations):
        simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = alts_main()
        print(count)
        #print(simulation_summary)
        survivors_simulations_summary.append(simulation_summary[0])
        total_simulations_summary.append(simulation_summary[1])
        count += 1

    combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    #print(combined_phenotypes)
    #print('RESUMEN DE SUPERVIVIENTES', survivors_simulations_summary)
    #print('RESUMEN DE TOTALES', total_simulations_summary)
    survivors_data_y = []
    proportions_data_y = []
    #data_y = []
    print('Calculating mean values')
    for phenotype in range(len(dict_phenotypes_combinations_indexes.keys())):
        #print('Con el fenotipo numero', dict_phenotypes_combinations_indexes.keys())
        survivors_means_per_phenotype = []
        total_means_per_phenotype = []
        for generation in generation_x:
            #print('Con la generacion numero', generation)
            survivors_mean = 0
            total_mean = 0
            for simulation in range(number_of_simulations):
                #print('Con la simulacion numero', simulation)
                #print('el valor es', simulations_summary[simulation][phenotype][generation])
                survivors_mean += survivors_simulations_summary[simulation][phenotype][generation]
                total_mean += total_simulations_summary[simulation][phenotype][generation]
            #print('MEDIAS',survivors_mean, '\t', total_mean)
            survivors_means_per_phenotype.append(survivors_mean/number_of_simulations)
            total_means_per_phenotype.append((total_mean/number_of_simulations)/population_size)
        survivors_data_y.append(survivors_means_per_phenotype)
        proportions_data_y.append(total_means_per_phenotype)
        #data_y.append(phenotype_data)
    #print('MEDIAS DE SUPERVIVIENTES', survivors_data_y)
    #print('MEDIAS DE TOTALES', proportions_data_y)
    #legend_population_size = 'N = '+str(population_size)
    #legend_phenotypes =dict_phenotypes_combinations_indexes.keys()

    total_survivors = []
    for generation in generation_x:
        individuals_sum = 0
        for phenotype_index in range(len(combined_phenotypes)):
            individuals_sum += survivors_data_y[phenotype_index][generation]
        total_survivors.append(individuals_sum)
    print(total_survivors)

    #plt.stackplot(generation_x, data_y, labels = legend_phenotypes)
    #plt.plot([], [], ' ', label = legend_population_size)

    #plt.margins(0)
    #plt.title('Proportion of individuals by phenotype')
    #plt.xlabel('Generation')
    #plt.ylabel('Proportion of individuals')
    #plt.legend()

    #plt.show()

    legend_population_size = 'N = '+str(population_size)
    legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    survivors_data_y.append(total_survivors)
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

create_simulation_results()
