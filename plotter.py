#!/usr/bin/python3

import pickle
import argparse

from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description = 'Plot results from alts.py script')
parser.add_argument('-i', '--input', default = 'result.alt', help = 'Input file where data will be extracted')
parser.add_argument('-d', '--display', action = 'store_true', help = 'Display the plots in popup windows, they will be '
                                                                     'saved .png if -q or --quiet is not used')
parser.add_argument('-q', '--quiet', action = 'store_false', help = 'Avoid saving the plots as .png files')
args = parser.parse_args()


def simulation_plot():
    with open(vars(args)['input'], 'rb') as input_file:
        simulations_summary = pickle.load(input_file)
    ''' simulations_summary is a list with:
    * population_size = The assigned size of the population
    * dict_phenotypes_combinations_indexes = The dictionary with the phenotypes and their indexes
    * Dictionary of gene's alleles
    * Survivors data by phenotype per generation per simulation
    * Proportions data by phenotype per generation per simulation
    * survivors_means = Mean survivors data by phenotype per generation of all simulations
    * proportions_means = Mean proportions data by phenotype per generation of all simulations
    * all_simulations_summary = Proportion of selfish individuals per generation in each simulation
    '''
    population_size = simulations_summary[0]
    dict_phenotypes_combinations_indexes = simulations_summary[1]
    generation_x = range(len(simulations_summary[3][0][0]))
    survivors_means = simulations_summary[5]
    if len(simulations_summary[6]) == 4:
        # Order of phenotypes is adjusted to be plotted in a more visual way
        proportions_means = [simulations_summary[6][2], simulations_summary[6][3], simulations_summary[6][0], simulations_summary[6][1]]
    else:
        proportions_means = simulations_summary[6]
    all_simulations_summary = simulations_summary[7]

    if len(simulations_summary[6]) == 4:
        # Order of phenotypes is adjusted to be plotted in a more visual way
        combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
        legend_phenotypes = [combined_phenotypes[2], combined_phenotypes[3], combined_phenotypes[0], combined_phenotypes[1]]
    else:
        legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())

    # Stackplot of the proportion of each phenotype at the beginning of the simulation
    proportions_plot = plt.figure(1)
    if len(simulations_summary[6]) == 4:
        # Order of phenotypes is adjusted to be plotted in a more visual way
        colour_map = ["#e74848ff", "#32be65ff", "#51bfe8ff", "#ee833c"]
        plt.stackplot(generation_x, proportions_means, labels = legend_phenotypes, colors = colour_map)
    else:
        colour_map = ["#51bfe8ff", "#e74848ff"]
        plt.stackplot(generation_x, proportions_means, labels = legend_phenotypes, colors = colour_map)
    plt.margins(0)
    plt.ylim(0,1)
    plt.title('Proportion of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')
    plt.legend()

    if args.quiet:
        plt.savefig('proportions_plot.png')

    # Total number of survivors is calculated to be added to the plot
    total_survivors = []
    for generation in generation_x:
        individuals_sum = 0
        for phenotype_index in range(len(legend_phenotypes)):
            individuals_sum += survivors_means[phenotype_index][generation]
        total_survivors.append(individuals_sum)

    survivors_means.append(total_survivors)
    colour_map.append('#000000ff')

    # Line plot of survivors per generation by phenotype
    survivors_plot = plt.figure(2)
    legend_phenotypes.append('total')
    for count in range(len(survivors_means)):
        plt.plot(survivors_means[count], label = legend_phenotypes[count], color = colour_map[count])
    plt.margins(0)
    plt.ylim((0,population_size))
    plt.title('Survivors per generation by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Number of individuals')
    plt.legend(loc = 'upper right')

    if args.quiet:
        plt.savefig('survivors_plot.png')

    # Line plot of the proportion of selfish individuals in every simulation
    individual_survivors_plot = plt.figure(3)
    for simulation in all_simulations_summary:
        plt.plot(simulation, color = 'black', linewidth = 0.5)
    plt.margins(0)
    plt.ylim(0, 1)
    plt.title('Selfish individuals per simulation')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')

    if args.quiet:
        plt.savefig('individual_survivors_plot.png')

    if args.display:
        plt.show()


simulation_plot()
