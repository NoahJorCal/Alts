import re
import pickle
import argparse

from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='Altruism simulations')
parser.add_argument('-i', '--input', default = 'result.alt', help = 'input file where data will be extracted')
args = parser.parse_args()

def simulation_plot():
    with open(vars(args)['input'], 'rb') as config_dictionary_file:
        simulations_summary = pickle.load(config_dictionary_file)

    population_size = simulations_summary[0]
    generation_x = range(len(simulations_summary[3][0][0]))

    dict_phenotypes_combinations_indexes = simulations_summary[1]
    dict_phenotype_options = simulations_summary[2]
    survivors_means = simulations_summary[5]
    #print(dict_phenotypes_combinations_indexes)

    if len(simulations_summary[6]) == 4:
        proportions_means = [simulations_summary[6][2], simulations_summary[6][3], simulations_summary[6][0], simulations_summary[6][1]]
    else:
        proportions_means = simulations_summary[6]

    all_simulations_summary = simulations_summary[7]

    if len(simulations_summary[6]) == 4:
        combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
        legend_phenotypes = [combined_phenotypes[2], combined_phenotypes[3], combined_phenotypes[0], combined_phenotypes[1]]
    else:
        legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())

    legend_population_size = 'N = '+str(population_size)
    proportions_plot = plt.figure(1)
    if len(simulations_summary[6]) == 4:
        color_map = ["#e74848ff", "#32be65ff", "#51bfe8ff", "#ee833c"]
        plt.stackplot(generation_x, proportions_means, labels = legend_phenotypes, colors = color_map)
    else:
        color_map = ["#51bfe8ff", "#e74848ff"]
        plt.stackplot(generation_x, proportions_means, labels = legend_phenotypes, colors = color_map)
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
    color_map.append('#000000ff')

    survivors_plot = plt.figure(2)
    legend_phenotypes.append('total')
    for count in range(len(survivors_means)):
        plt.plot(survivors_means[count], label = legend_phenotypes[count], color = color_map[count])
    plt.margins(0)
    plt.ylim((0,population_size))
    plt.title('Survivors per generation by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Number of individuals')
    plt.legend(loc = 'upper right')

    individual_survivors_plot = plt.figure(3)
    for simulation in all_simulations_summary:
        plt.plot(simulation, color = 'black', linewidth = 0.5)
    plt.margins(0)
    plt.ylim(0,1)
    plt.title('Selfish individuals per simulation')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')



    plt.show()

simulation_plot()
