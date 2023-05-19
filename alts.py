#!/usr/bin/python3

import argparse
import os
import random
import multiprocessing
from configparser import ConfigParser
from time import perf_counter

import numpy as np
from numpy import mean

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
parser.add_argument('-d', '--directory', default='output', help='Output directory where individual'
                                                                'simulation results will be stored')
parser.add_argument('-o', '--output', default='simulation.h5', help='Output file where the'
                                                                    'simulation results will be stored')
parser.add_argument('-c', '--cpu', default=1, type=int, help='Number of simultaneous workers')
parser.add_argument('-q', '--quiet', action='store_true', help='Number of simultaneous workers')
parser.add_argument('-s', '--seed', required=False, type=int, help='Set a seed for the simulations')
args = parser.parse_args()


def run_simulation(output_dir, output_file, seed, quiet):
    # start_counter = perf_counter()
    aborted_sim = simulator_main(output_dir, output_file, seed, quiet)
    # simulation_duration = perf_counter() - start_counter
    if aborted_sim:
        if not args.quiet:
            print(f'\033[K\033[FSimulation ended because altruism went extinct')
        return None
    else:
        return True
    # else:
    #     if not args.quiet:
    #         print(f'\033[K\033[FSimulations up to {simulation_index + 1} run in {round(simulation_duration, 2)} seconds. '
    #               f'Running next round with {args.cpu} simulations...')
    #     return simulation_duration


def create_simulation_results():
    outputs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
    # if os.path.exists(outputs_dir):
    #     shutil.rmtree(outputs_dir)

    # number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    if not args.quiet:
        print(f'Running first round of {args.cpu} simulations...')

    if args.cpu > os.cpu_count():
        raise CPUError()
    with multiprocessing.Pool(processes=args.cpu) as pool:
        # results = [pool.apply_async(run_simulation, args=(x, )) for x in
        #            range(number_of_simulations)]
        results = [pool.apply_async(run_simulation, args=(args.directory, args.output, args.seed, args.quiet))]
        output = [p.get() for p in results]
    output = [i for i in output if i is not None]
    output_file = args.output
    if '.h5' not in output_file and '.hdf5' not in output_file and \
       '.h5p' not in output_file and '.he5' not in output_file and \
       '.h5m' not in output_file and '.h5z' not in output_file:
        output_file += '.h5'
    simulation_results = os.path.join(os.path.dirname(__file__), args.directory, output_file)
    if output:
        if not args.quiet:
            print(f'\033[K\033[F\033[KEach simulation took {round(mean(output), 2)} seconds on average and '
                  f'{round((perf_counter() - start_time)/60,2)} minutes in total')
    else:
        os.remove(simulation_results)
        if not args.quiet:
            # print('\033[1A', end='\x1b[2K')
            # print('\033[1A', end='\x1b[2K')
            # print('\033[1A', end='\x1b[2K')
            print('Altruism went extinct in all simulations, no data generated')
            print('\033[1A', end='\x1b[2K')
            print('\033[1A', end='\x1b[2K')
            print('\033[1A', end='\x1b[2K')
        return True

    if not args.quiet:
        print(f'\x1b[2KThe results of all the simulations have been saved in {simulation_results}')
    return False

    def print_name_type(name, obj):
        print(name, type(obj))

    # print(f'Saving all simulations results in {args.output}...')
    # with h5py.File(args.output, 'w') as out_file:
    #
    #     # Loop through the list of input files
    #     for i in range(number_of_simulations):
    #
    #         # Open the input file in "read" mode
    #         with h5py.File(f'outputs/simulation_{i}.h5', 'r') as in_file:
    #             # Create a new group for the input file in the output file
    #             out_group = out_file.create_group(f'simulation_{i}')
    #
    #             # Loop through each group in the input file
    #             for group_name in in_file:
    #
    #                 # Create a new group in the output file with the same name as the input group
    #                 out_subgroup = out_group.create_group(group_name)
    #
    #                 # Loop through each dataset in the input group and copy it to the output subgroup
    #                 for dataset_name in in_file[group_name]:
    #                     in_dataset = in_file[group_name][dataset_name]
    #                     out_subgroup.create_dataset(dataset_name, data=in_dataset)
    #             out_file.attrs['simulations'] = in_file.attrs['simulations']
    #             out_file.attrs['generations'] = in_file.attrs['generations']
    #             out_file.attrs['n_loci'] = in_file.attrs['n_loci']
    #             out_file.attrs['groups'] = in_file.attrs['groups']
    #             out_file.attrs['loci'] = in_file.attrs['loci']
    #             out_file.attrs['phenotype_names'] = in_file.attrs['phenotype_names']
    #             out_file.attrs['alleles_names'] = in_file.attrs['alleles_names']
    # if os.path.exists(outputs_dir):
    #     shutil.rmtree(outputs_dir)
    # print(f'\r\033[F\rThe results of all the simulations have been saved in {args.output}\n'
    #       f'All intermediate files have been removed')


if __name__ == '__main__':
    aborted = True
    while aborted:
        aborted = create_simulation_results()
        # print(aborted)

# def create_simulation_results():
#     number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
#     generation_x = range(int(general_config['simulation']['generations']) + 1)
#     simulation_iter = itertools.count()
#     round_count = 0
#     print(f'Running first round of {args.cpu} simulations...')
#     mean_simulation_duration = 0
#
#     if args.cpu > os.cpu_count():
#         raise CPUError()
#
#     for i in range(int(number_of_simulations/args.cpu + 0.5)):
#         workers = [None] * args.cpu
#         returns = [None] * args.cpu
#         # The workers are initialized to run one simulation each
#         for worker_index in range(args.cpu):
#             simulation_index = next(simulation_iter)
#             worker = threading.Thread(
#                 target=run_simulation,
#                 args=(generation_x,
#                       round_count,
#                       returns,
#                       worker_index,
#                       simulation_index))
#
#             # noinspection PyTypeChecker
#             workers[worker_index] = worker
#             worker.start()
#
#         for worker in workers:
#             worker.join()
#
#         # After all the workers have finished, the data from each simulation is saved
#         for worker_return in returns:
#             data_file, simulation_duration = worker_return
#             mean_simulation_duration += simulation_duration
#
#         round_count += 1
#
#     ''' The minimum number of simulations will be the number given, but if there will be free CPUs,
#     more simulations will be run '''
#     number_of_simulations = (int(number_of_simulations/args.cpu+0.5))*args.cpu
#
#     print('\033[K\033[F\033[K\033[F\r')
#
#     mean_simulation_duration = mean_simulation_duration/number_of_simulations
#     print(f'\rEach simulation took {round(mean_simulation_duration,2)} seconds on average and '
#           f'{round((perf_counter() - start_time)/60,2)} minutes in total')
#     return data_file
#
#
# if __name__ == '__main__':
#     create_simulation_results()
#     # with h5py.File('simulation_output.h5', 'r') as f:
#     #     print(list(f.keys()))
