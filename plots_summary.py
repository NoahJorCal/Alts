import re
from time import perf_counter
import threading
import argparse
import sys
import pickle

from configparser import ConfigParser

from alts import alts_main

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

start_time = perf_counter()

parser = argparse.ArgumentParser(description='Altruism simulations')
parser.add_argument('-o', '--outfile', default = 'result.txt', help = 'output file where data will be stored')
args = parser.parse_args()

def run_simulation(generation_x, survivors_simulations_summary, total_simulations_summary):
    simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = alts_main()
    survivors_simulations_summary.append(simulation_summary[0])
    total_simulations_summary.append(simulation_summary[1])
    return simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options

def create_simulation_results():
    number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    generation_x = range(int(general_config['simulation']['generations'])+1)
    simulation_count = 0
    print('Running first simulation...')

    mean_simulation_duration = 0
    #The given number of simulations are run and their data is saved
    total_simulations_summary = []
    survivors_simulations_summary = []
    all_simulations_summary = []
    simulation_duration = 0
    for i in range(number_of_simulations):
        start_counter = perf_counter()
        simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = run_simulation(generation_x, survivors_simulations_summary, total_simulations_summary)

        phenotypes = list(dict_phenotypes_combinations_indexes.keys())
        selfish_indexes = []
        single_simulations_summary = [0 for generation in generation_x]
        for phentoype_index in range(len(phenotypes)):
            if re.search('(?:&|^)(selfish)(?:&|$)', phenotypes[phentoype_index]):
                selfish_indexes.append(phentoype_index)
        for phenotype in selfish_indexes:
            single_simulations_summary = [sum(x) for x in zip(single_simulations_summary, simulation_summary[1][phenotype])]
        all_simulations_summary.append(single_simulations_summary)

        print(f'\033[K\033[FSimulation number {simulation_count} run in {round(perf_counter() - start_counter, 2)} seconds. Running simulation {simulation_count + 1}...')
        simulation_count += 1
        simulation_duration += perf_counter() - start_counter
        mean_simulation_duration += simulation_duration
        start_counter = perf_counter()

    combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    survivors_means = []
    proportions_means = []
    for phenotype in range(len(combined_phenotypes)):
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
        survivors_means.append(survivors_means_per_phenotype)
        proportions_means.append(total_means_per_phenotype)

    print('\033[K\033[F\033[K\033[F\r')

    mean_simulation_duration = mean_simulation_duration/number_of_simulations
    print(f'\rEach simulation took {round(mean_simulation_duration,2)} seconds on average and {round((perf_counter() - start_time)/60,2)} minutes in total')

    return dict_phenotypes_combinations_indexes, dict_phenotype_options, survivors_simulations_summary, total_simulations_summary, survivors_means, proportions_means, all_simulations_summary

with open(vars(args)['outfile'], 'wb') as config_dictionary_file:
  pickle.dump(create_simulation_results(), config_dictionary_file)
