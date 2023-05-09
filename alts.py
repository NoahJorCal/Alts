#!/usr/bin/python3

import argparse
import os
import shutil
import multiprocessing
from configparser import ConfigParser
from time import perf_counter
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
parser.add_argument('-o', '--output', default='simulation_output.h5', help='Output file where data will be stored')
parser.add_argument('-c', '--cpu', default=1, type=int, help='Number of simultaneous workers')
args = parser.parse_args()


def run_simulation(simulation_index):
    start_counter = perf_counter()
    simulator_main(args.output, simulation_index)
    simulation_duration = perf_counter() - start_counter
    print(f'\033[K\033[FSimulations up to {simulation_index + 1} run in {round(simulation_duration, 2)} seconds.'
          f'Running next round with {args.cpu} simulations...')
    return simulation_duration


def create_simulation_results():
    number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    print(f'Running first round of {args.cpu} simulations...')

    if args.cpu > os.cpu_count():
        raise CPUError()
    with multiprocessing.Pool(processes=args.cpu) as pool:
        results = [pool.apply_async(run_simulation, args=(x, )) for x in
                   range(number_of_simulations)]
        output = [p.get() for p in results]
    print(f'\033[K\033[F\033[KEach simulation took {round(mean(output), 2)} seconds on average and '
          f'{round((perf_counter() - start_time)/60,2)} minutes in total')

    def print_name_type(name, obj):
        print(name, type(obj))
    # output = [p.get() for p in results]
    # with multiprocessing.Pool(processes=args.cpus) as pool:
    #     results = [pool.apply_async(run_simulation, args=(x,)) for x in range(5)]
    #     output = [p.get() for p in results]
    # with h5py.File('outputs/simulation_0', 'r') as f:
    #     f.visititems(print_name_type)
    # with h5py.File('outputs/simulation_1', 'r') as f:
    #     f.visititems(print_name_type)
    # with h5py.File('outputs/simulation_2', 'r') as f:
    #     f.visititems(print_name_type)
    # with h5py.File('outputs/simulation_3', 'r') as f:
    #     f.visititems(print_name_type)
    # with h5py.File('outputs/simulation_4', 'r') as f:
    #     f.visititems(print_name_type)
    # print(args.output)

    print(f'Saving all simulations results in {args.output}...')
    with h5py.File(args.output, 'w') as out_file:

        # Loop through the list of input files
        for i in range(number_of_simulations):

            # Open the input file in "read" mode
            with h5py.File(f'outputs/simulation_{i}.h5', 'r') as in_file:
                # Create a new group for the input file in the output file
                out_group = out_file.create_group(f'simulation_{i}')

                # Loop through each group in the input file
                for group_name in in_file:

                    # Create a new group in the output file with the same name as the input group
                    out_subgroup = out_group.create_group(group_name)

                    # Loop through each dataset in the input group and copy it to the output subgroup
                    for dataset_name in in_file[group_name]:
                        in_dataset = in_file[group_name][dataset_name]
                        out_subgroup.create_dataset(dataset_name, data=in_dataset)
    outputs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outputs')
    if os.path.exists(outputs_dir):
        shutil.rmtree(outputs_dir)
    print(f'\033[K\033[F\033[KThe results of all the simulations have been saved in {args.output}\n'
          f'All intermediate files have been removed')

if __name__ == '__main__':
    create_simulation_results()

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
