#!/usr/bin/python3

import argparse
import itertools
import os
import threading
from configparser import ConfigParser
from time import perf_counter

import h5py

from simulator import simulator_main


class CPUError(Exception):
    def __init__(self):
        super().__init__()

    def __str__(self):
        return f'Set {args.cpu} cpus when only {os.cpu_count()} are available'


# Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

start_time = perf_counter()

parser = argparse.ArgumentParser(description='Altruism simulations')
parser.add_argument('-o', '--output', default='simulation_output.h5', help='Output file where data will be stored')
parser.add_argument('-c', '--cpu', default=1, type=int, help='Number of simultaneous workers')
args = parser.parse_args()


def run_simulation(
        generation_x,
        round_count,
        returns,
        worker_index,
        simulation_index):

    start_counter = perf_counter()
    data_file = simulator_main(args.output, simulation_index)
    simulation_duration = perf_counter() - start_counter

    # phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    # selfish_indexes = []
    # # Only selfish individuals, if other genes are present they are ignored
    # single_simulations_summary = [0 for generation in generation_x]
    # for phenotype_index in range(len(phenotypes)):
    #     if re.search('(?:&|^)(selfish)(?:&|$)', phenotypes[phenotype_index]):
    #         selfish_indexes.append(phenotype_index)
    # # Counts of selfish vs non-selfish individuals
    # for phenotype in selfish_indexes:
    #     single_simulations_summary = [sum(x) for x in zip(single_simulations_summary, simulation_summary[1][phenotype])]

    print(f'\033[K\033[FRound number {round_count} run in {round(simulation_duration, 2)} seconds.'
          f'Running round {round_count + 1} with {args.cpu} simulations...')

    returns[worker_index] = (data_file, simulation_duration)


def create_simulation_results():
    number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    generation_x = range(int(general_config['simulation']['generations']) + 1)
    simulation_iter = itertools.count()
    round_count = 0
    print(f'Running first round of {args.cpu} simulations...')
    mean_simulation_duration = 0

    if args.cpu > os.cpu_count():
        raise CPUError()

    for i in range(int(number_of_simulations/args.cpu + 0.5)):
        workers = [None] * args.cpu
        returns = [None] * args.cpu
        # The workers are initialized to run one simulation each
        for worker_index in range(args.cpu):
            simulation_index = next(simulation_iter)
            worker = threading.Thread(
                target=run_simulation,
                args=(generation_x,
                      round_count,
                      returns,
                      worker_index,
                      simulation_index))

            # noinspection PyTypeChecker
            workers[worker_index] = worker
            worker.start()

        for worker in workers:
            worker.join()

        # After all the workers have finished, the data from each simulation is saved
        for worker_return in returns:
            data_file, simulation_duration = worker_return
            mean_simulation_duration += simulation_duration

        round_count += 1

    ''' The minimum number of simulations will be the number given, but if there will be free CPUs,
    more simulations will be run '''
    number_of_simulations = (int(number_of_simulations/args.cpu+0.5))*args.cpu

    print('\033[K\033[F\033[K\033[F\r')

    mean_simulation_duration = mean_simulation_duration/number_of_simulations
    print(f'\rEach simulation took {round(mean_simulation_duration,2)} seconds on average and '
          f'{round((perf_counter() - start_time)/60,2)} minutes in total')
    return data_file


if __name__ == '__main__':
    create_simulation_results()
    # with h5py.File('simulation_output.h5', 'r') as f:
    #     print(list(f.keys()))
