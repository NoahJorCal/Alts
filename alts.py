import re
from time import perf_counter
import threading
import argparse
import sys
import pickle
import os

from configparser import ConfigParser

from simulator import simulator_main

class CPUError(Exception):
    def __init__(self):
        super().__init__()

    def __str__(self):
        return f'Set {args.cpu} cpus when only {os.cpu_count()} are available'

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

start_time = perf_counter()

parser = argparse.ArgumentParser(description='Altruism simulations')
parser.add_argument('-o', '--outfile', default = 'result.txt', help = 'Output file where data will be stored')
parser.add_argument('-c', '--cpu', default = 1, type = int, help = 'Number of simultanious workers')
#parser.add_argument('-p', '--plot', default = True, type = bool, help = 'If plot will be produced')
args = parser.parse_args()


def run_simulation(
        generation_x,
        round_count,
        returns,
        worker_index,
        survivors_simulations_summary,
        total_simulations_summary,
        all_simulations_summary,
        mean_simulation_duration):

    start_counter = perf_counter()
    simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = simulator_main()

    phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    selfish_indexes = []
    single_simulations_summary = [0 for generation in generation_x]
    for phentoype_index in range(len(phenotypes)):
        if re.search('(?:&|^)(selfish)(?:&|$)', phenotypes[phentoype_index]):
            selfish_indexes.append(phentoype_index)
    for phenotype in selfish_indexes:
        single_simulations_summary = [sum(x) for x in zip(single_simulations_summary, simulation_summary[1][phenotype])]

    simulation_duration = perf_counter() - start_counter
    print(f'\033[K\033[FRound number {round_count} run in {round(simulation_duration, 2)} seconds. Running round {round_count + 1} with {args.cpu} simulations...')

    returns[worker_index] = (simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options, single_simulations_summary, simulation_duration)

def create_simulation_results():
    number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    generation_x = range(int(general_config['simulation']['generations'])+1)
    round_count = 0
    print(f'Running first round of {args.cpu} simulations...')

    mean_simulation_duration = 0
    #The given number of simulations are run and their data is saved
    total_simulations_summary = []
    survivors_simulations_summary = []
    all_simulations_summary = []
    simulation_duration = 0

    if args.cpu > os.cpu_count():
        raise CPUError()

    for i in range(int(number_of_simulations/args.cpu+0.5)):
        workers = [None] * args.cpu
        returns = [None] * args.cpu
        for worker_index in range(args.cpu):
            worker = threading.Thread(
                target = run_simulation,
                args = (
                    generation_x,
                    round_count,
                    returns,
                    worker_index,
                    survivors_simulations_summary,
                    total_simulations_summary,
                    all_simulations_summary,
                    mean_simulation_duration))

            workers[worker_index] = worker
            worker.start()
        for worker in workers:
            worker.join()

        for worker_return in returns:
            simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options, single_simulations_summary, simulation_duration = worker_return

            survivors_simulations_summary.append(simulation_summary[0])
            total_simulations_summary.append(simulation_summary[1])
            all_simulations_summary.append(single_simulations_summary)

            mean_simulation_duration += simulation_duration

        round_count += 1

    number_of_simulations = (int(number_of_simulations/args.cpu+0.5))*args.cpu
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

with open(args.outfile, 'wb') as config_dictionary_file:
  pickle.dump(create_simulation_results(), config_dictionary_file)

#if args.plot:
    #from plotter import plot
    #simulation_plot()

