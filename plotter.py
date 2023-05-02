#!/usr/bin/python3
import os
from configparser import ConfigParser
import argparse
import h5py
import numpy as np
from matplotlib import pyplot as plt

# Import general configuration
general_config = ConfigParser()
config_path = os.path.join(os.path.dirname(__file__), 'config.ini')
general_config.read(config_path)

parser = argparse.ArgumentParser(description='Plot results from alts.py script')
parser.add_argument('-i', '--input', default='simulation_output.h5', help='Input file where data will be extracted')
parser.add_argument('-s', '--show', action='store_true', help='Show the plots in popup windows, they will be '
                                                              'saved .png if -ns or --ns is not used')
parser.add_argument('-ns', '--no-save', dest='save_plots', action='store_false',
                    help='Avoid saving the plots as .png files')
args = parser.parse_args()


def print_name_type(name, obj):
    print(name, type(obj))


# Proportion of individuals of each phenotype at the beginning of each generation
phenotypes_pgen_bool = general_config['output']['phenotypes_pgen'] == 'True'
# Proportion of alleles of each locus at the beginning of each generation
allele_pgen_bool = general_config['output']['allele_pgen'] == 'True'
# Proportion of selfish individuals at the beginning of each generation in each simulation
selfish_phenotypes_pgen_asim_bool = general_config['output']['selfish_phenotypes_pgen_asim'] == 'True'
# Proportion of selfish alleles at the beginning of each generation in each simulation
selfish_alleles_pgen_asim_bool = general_config['output']['selfish_alleles_pgen_asim'] == 'True'
# Altruist to selfish ratio in groups at the beginning of each generation
altruist_selfish_ratio_pgen_bool = general_config['output']['altruist_selfish_ratio_pgen'] == 'True'
# Number of total survivors by the end of the generation
survivors_pgen_bool = general_config['output']['survivors_pgen'] == 'True'
# Number of total survivors by the end of the generation based on altruist/selfish group ratio
survivors_pgr_pgen_bool = general_config['output']['survivors_pgr_pgen'] == 'True'
# Haplotypes of each locus
haplotypes_bool = general_config['output']['haplotypes'] == 'True'


def save_phenotypes_pgen(phenotype, n_loci, generation, simulation, phenotypes_list):
    simulation = int(simulation.split('_')[-1])
    generation = int(generation.split('_')[-1])
    for i in range(n_loci):
        phenotypes_list[simulation][generation][i] = np.count_nonzero(phenotype == i)


def save_allele_pgen_bool(alleles, n_loci, generation, simulation, alleles_names, alleles_list):
    simulation = int(simulation.split('_')[-1])
    generation = int(generation.split('_')[-1])
    for locus in range(n_loci):
        for i in range(len(alleles_names[locus])):
            alleles_list[locus][simulation][generation][i] = np.count_nonzero(alleles[locus] == i)


def save_altruist_selfish_ratio_pgen(group, phenotype, generation, simulation, as_ratio_list):
    simulation = int(simulation.split('_')[-1])
    generation = int(generation.split('_')[-1])
    # print(as_ratio_list, end='\n\n')
    # Find indices where a change in group occurs
    indices = np.where(group[:-1] != group[1:])[0] + 1
    indices = np.concatenate((np.array([0]), indices, np.array([len(group)])))
    # print(indices)
    for i in range(len(indices) - 1):
        altruists = np.count_nonzero(phenotype[indices[i]:indices[i + 1]])
        selfish = len(phenotype[indices[i]:indices[i + 1]]) - altruists
        as_ratio_list[simulation][generation][i] = selfish/altruists




def simulation_plot():
    with h5py.File(vars(args)['input'], 'r') as f:
        simulations = f.attrs['simulations']
        generations = f.attrs['generations']
        n_loci = f.attrs['n_loci']
        groups = f.attrs['groups']
        loci = f.attrs['loci']
        alleles_names = []
        for alleles in f.attrs['alleles_names']:
            alleles_names.append(np.flip(alleles, axis=0))
        alleles_names = np.array(alleles_names)
        phenotype_names = np.flip(f.attrs['phenotype_names'], axis=0)
        phenotypes_all = np.zeros((simulations, generations, n_loci))
        as_ratio_all = np.zeros((simulations, generations, groups))
        # ASSUMING THAT EACH LOCUS HAS TWO ALLELES
        alleles_all = np.zeros((n_loci, simulations, generations, 2))
        # print(phenotypes_all)
        # print(list(f.keys()))
        for simulation in list(f.keys()):
            # phenotypes_sim = [[] for _ in f[simulation]['generation_0']['phenotype'].attrs['phenotype_names']]
            for generation in (generation for generation in list(f[simulation].keys()) if 'generation' in generation):
                # print(list(f[simulation][generation].keys()))
                phenotype = f[simulation][generation]['phenotype'][:]
                group = f[simulation][generation]['group'][:]
                # print(phenotype)
                alleles = []
                for locus in loci:
                    alleles.append(f[simulation][generation][f'locus_{locus}'][:])
                # print(alleles)
                # print(simulation, generation)
                if phenotypes_pgen_bool:
                    save_phenotypes_pgen(phenotype, n_loci, generation, simulation, phenotypes_all)
                if allele_pgen_bool:
                    save_allele_pgen_bool(alleles, n_loci, generation, simulation, alleles_names, alleles_all)
                if altruist_selfish_ratio_pgen_bool:
                    save_altruist_selfish_ratio_pgen(group, phenotype, generation, simulation, as_ratio_all)

        if phenotypes_pgen_bool:
            phenotypes_mean = np.mean(phenotypes_all, axis=0)
            phenotypes_proportions = phenotypes_mean / np.sum(phenotypes_mean, axis=1, keepdims=True)
            phenotypes_proportions = np.flip(phenotypes_proportions, axis=1).T

        if allele_pgen_bool:
            alleles_mean = np.mean(alleles_all, axis=1)
            alleles_proportions = alleles_mean / np.sum(alleles_mean, axis=2, keepdims=True)
            alleles_proportions = np.flip(alleles_proportions, axis=2)
            trans_alleles_proportions = []
            for locus_proportions in alleles_proportions:
                trans_alleles_proportions.append(locus_proportions.T)
            alleles_proportions = np.array(trans_alleles_proportions)

        if altruist_selfish_ratio_pgen_bool:
            as_ratio_mean = np.mean(np.mean(as_ratio_all, axis=0), axis=1)
            # print(as_ratio_mean, end='\n\n')
            # # phenotypes_proportions = phenotypes_mean / np.sum(phenotypes_mean, axis=1, keepdims=True)
            # # phenotypes_proportions = np.flip(phenotypes_proportions, axis=1).T

    generation_x = np.arange(generations)

    if phenotypes_pgen_bool:
        # Stack-plot of the proportion of each phenotype at the beginning of each generation
        phenotypes_plot = plt.figure(1)
        colour_map = ["#51bfe8ff", "#e74848ff"]
        plt.stackplot(generation_x, phenotypes_proportions, labels=phenotype_names, colors=colour_map)
        plt.margins(0)
        plt.title('Proportion of individuals by phenotype')
        plt.xlabel('Generation')
        plt.ylabel('Proportion of individuals')
        plt.legend()

    if allele_pgen_bool:
        # Stack-plot of the proportion of each allele at the beginning of each generation
        alleles_plot = plt.figure(2)
        colour_map = ["#51bfe8ff", "#e74848ff"]
        plt.stackplot(generation_x, alleles_proportions[0], labels=alleles_names[0], colors=colour_map)
        plt.margins(0)
        plt.title('Proportion alleles')
        plt.xlabel('Generation')
        plt.ylabel('Proportion of alleles')
        plt.legend()

    if altruist_selfish_ratio_pgen_bool:
        # Stack-plot of the ratio selfish/altruist per generation
        as_ratio_plot = plt.figure(3)
        colour_map = ["#51bfe8ff", "#e74848ff"]
        print(generation_x)
        print(as_ratio_mean)
        plt.plot(generation_x, as_ratio_mean)
        plt.margins(0)
        plt.title('Proportion alleles')
        plt.xlabel('Generation')
        plt.ylabel('Proportion of alleles')
        plt.ylim(0)
        plt.legend()


    if args.save_plots:
        plt.savefig('proportions_plot.png')
    if args.show:
        plt.show()









#     population_size = simulations_summary[0]
#     dict_phenotypes_combinations_indexes = simulations_summary[1]
#     generation_x = range(len(simulations_summary[3][0][0]))
#     survivors_means = simulations_summary[5]
#     if len(simulations_summary[6]) == 4:
#         # Order of phenotypes is adjusted to be plotted in a more visual way
#         proportions_means = [simulations_summary[6][2], simulations_summary[6][3], simulations_summary[6][0], simulations_summary[6][1]]
#     else:
#         proportions_means = simulations_summary[6]
#     all_simulations_summary = simulations_summary[7]
#
#     if len(simulations_summary[6]) == 4:
#         # Order of phenotypes is adjusted to be plotted in a more visual way
#         combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
#         legend_phenotypes = [combined_phenotypes[2], combined_phenotypes[3], combined_phenotypes[0], combined_phenotypes[1]]
#     else:
#         legend_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
#
#     # Stackplot of the proportion of each phenotype at the beginning of the simulation
#     proportions_plot = plt.figure(1)
#     if len(simulations_summary[6]) == 4:
#         # Order of phenotypes is adjusted to be plotted in a more visual way
#         colour_map = ["#e74848ff", "#32be65ff", "#51bfe8ff", "#ee833c"]
#         plt.stackplot(generation_x, proportions_means, labels = legend_phenotypes, colors = colour_map)
#     else:
#         colour_map = ["#51bfe8ff", "#e74848ff"]
#         plt.stackplot(generation_x, proportions_means, labels = legend_phenotypes, colors = colour_map)
#     plt.margins(0)
#     plt.ylim(0, 1)
#     plt.title('Proportion of individuals by phenotype')
#     plt.xlabel('Generation')
#     plt.ylabel('Proportion of individuals')
#     plt.legend()
#
#     if not args.save_plots:
#         plt.savefig('proportions_plot.png')
#
#     # Total number of survivors is calculated to be added to the plot
#     total_survivors = []
#     for generation in generation_x:
#         individuals_sum = 0
#         for phenotype_index in range(len(legend_phenotypes)):
#             individuals_sum += survivors_means[phenotype_index][generation]
#         total_survivors.append(individuals_sum)
#
#     survivors_means.append(total_survivors)
#     colour_map.append('#000000ff')
#
#     # Line plot of survivors per generation by phenotype
#     survivors_plot = plt.figure(2)
#     legend_phenotypes.append('total')
#     for count in range(len(survivors_means)):
#         plt.plot(survivors_means[count], label=legend_phenotypes[count], color=colour_map[count])
#     plt.margins(0)
#     plt.ylim((0, population_size))
#     plt.title('Survivors per generation by phenotype')
#     plt.xlabel('Generation')
#     plt.ylabel('Number of individuals')
#     plt.legend(loc='upper right')
#
#     if not args.save_plots:
#         plt.savefig('survivors_plot.png')
#
#     # Line plot of the proportion of selfish individuals in every simulation
#     individual_survivors_plot = plt.figure(3)
#     for simulation in all_simulations_summary:
#         plt.plot(simulation, color='black', linewidth=0.5)
#     plt.margins(0)
#     plt.ylim(0, 1)
#     plt.title('Selfish individuals per simulation')
#     plt.xlabel('Generation')
#     plt.ylabel('Proportion of individuals')
#
#     if not args.save_plots:
#         plt.savefig('individual_survivors_plot.png')
#

#
#
simulation_plot()
