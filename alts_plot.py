from configparser import ConfigParser
from matplotlib import pyplot as plt
from alts import alts_main
import re
from itertools import product

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

def plot():
    #Run the main function of the alts script, runs a simulation and generates a summary
    simulation_summary, alleles_combinations_indexes, dict_allele_options = alts_main()
    population_size = int(general_config['population']['size'])
    genes = general_config['plot']['genes'].replace(' ','').split(',')
    generation_x = list(range(int(general_config['simulation']['generations'])))
    plot_info = {}
    allele_options = []
    for gene in genes:
        allele_options.append(dict_allele_options[gene])
    for combination in product(*allele_options):
        match_indexes = list(range(len(simulation_summary)))
        phenotype_name = ''
        for element in combination:
            new_match_indexes = []
            for index in match_indexes:
                phenotype = list(alleles_combinations_indexes.keys())[index]
                if re.search('(?:&|^)('+element+')(?:&|$)', phenotype):
                    new_match_indexes.append(index)
            match_indexes = new_match_indexes
            phenotype_name += element+'&'
        count_list = [0]*int(general_config['simulation']['generations'])
        for generation in generation_x:
            for index in match_indexes:
                count_list[generation] += simulation_summary[index][generation]/population_size
            plot_info[phenotype_name[:-1]] = count_list
        #print(plot_info)

    data_y = plot_info.values()
    legend_phenotypes = plot_info.keys()

    plt.stackplot(generation_x, data_y, labels = legend_phenotypes)

    plt.margins(0)
    plt.title('Proportion of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')
    plt.legend(legend_phenotypes)

    plt.show()

plot()
