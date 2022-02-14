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
    generation_x = range(int(general_config['simulation']['generations']))
    simulations_summary = []
    #The given number of simulations are run and their data is saved
    count = 0
    for i in range(number_of_simulations):
        print(count)
        simulation_summary, alleles_combinations_indexes, dict_allele_options = alts_main()
        simulations_summary.append(simulation_summary)
        count += 1
    mean_values = []
    data_y = []
    print('Calculating mean values')
    for phenotype in range(len(alleles_combinations_indexes.keys())):
        #print('Con el fenotipo numero', phenotype)
        phenotype_data = []
        for generation in generation_x:
            #print('Con la generacion numero', generation)
            mean = 0
            for simulation in range(number_of_simulations):
                #print('Con la simulacion numero', simulation)
                #print('el valor es', simulations_summary[simulation][phenotype][generation])
                mean += simulations_summary[simulation][phenotype][generation]
            phenotype_data.append((mean/number_of_simulations)/population_size)
        data_y.append(phenotype_data)

    legend_population_size = 'N = '+str(population_size)
    legend_phenotypes = alleles_combinations_indexes.keys()

    plt.stackplot(generation_x, data_y, labels = legend_phenotypes)
    plt.plot([], [], ' ', label = legend_population_size)

    plt.margins(0)
    plt.title('Proportion of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')
    plt.legend()

    plt.show()

create_simulation_results()
