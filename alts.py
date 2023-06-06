#!/usr/bin/python3

from configparser import ConfigParser
import os
import argparse

from time import perf_counter

import multiprocessing
from numpy import mean

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
parser.add_argument('-ns', '--no-save', dest='save_data', action='store_false',
                    help='The data will not be collected')
parser.add_argument('-d', '--directory', default='output', help='Output directory where individual'
                                                                'simulation results will be stored')
parser.add_argument('-o', '--output', default='simulation.h5', help='Output file where the'
                                                                    'simulation results will be stored')
parser.add_argument('-g', '--genome', dest='simulate_genome', action='store_true',
                    help='The genomes will be simulated')
parser.add_argument('-c', '--cpu', default=1, type=int, help='Number of simultaneous workers')
parser.add_argument('-q', '--quiet', action='store_true', help='Number of simultaneous workers')
parser.add_argument('-s', '--seed', required=False, type=int, help='Set a seed for the simulations')
args = parser.parse_args()


def run_simulation(save_data, output_dir, output_file, simulate_genome, seed, quiet):
    """
    Runs one simulation from simulator.py script.
    :param bool save_data: If True, saves the data of each simulation in the output file.
    :param str output_dir: Output directory name where the H5PY files will be stored.
    :param str output_file: Output H5PY file name where the data will be stored.
    :param bool simulate_genome: If True, the genomes with SNVs will be simulated.
    :param int seed: Seed used to fix simulation's result, defaults to None.
    :param bool quiet: Flag to avoid printing the program's feedback, defaults to False.
    :return: A tuple with None if the simulation ended because of lack of altruists or individuals or the duration of
    the simulation if it ended without any problems, and the name of the output file,
    as it could have changed because of uniquify.
    """
    start_counter = perf_counter()
    stop, output = simulator_main(save_data, output_dir, output_file, simulate_genome, seed, quiet)
    simulation_duration = perf_counter() - start_counter
    if stop:
        if not args.quiet:
            print(f'\033[K\033[FSimulation ended because all individuals died or altruism went extinct')
        return None, output
    else:
        return simulation_duration, output


def create_simulation_results():
    """
    Manages multiprocessing of parallel calls of the simulation.
    :return:
    """
    if not args.quiet:
        print(f'Running first round of {args.cpu} simulations...')
    # If the configured CPUs to use is higher than the machine's CPUs
    if args.cpu > os.cpu_count():
        raise CPUError()
    with multiprocessing.Pool(processes=args.cpu) as pool:
        results = [pool.apply_async(run_simulation, args=(args.save_data, args.directory, args.output,
                                                          args.simulate_genome, args.seed, args.quiet))]
        output = [p.get() for p in results]
    if args.save_data:
        out_file = os.path.join(os.path.dirname(__file__), args.directory, output[0][1])
    if output[0][0]:
        with open('output.csv', 'a') as out_file:
            out_file.write(str(output[0][1])[1:-1].replace(' ', '') + '\n')
        if not args.quiet:
            minutes_avg, seconds_avg = divmod(mean(output[0][0]), 60)
            hours_avg, minutes_avg = divmod(minutes_avg, 60)
            minutes, seconds = divmod(perf_counter() - start_time, 60)
            hours, minutes = divmod(minutes, 60)
            print(f'\033[K\033[F\033[KEach simulation took {round(hours_avg):d}:{round(minutes_avg):02d}:'
                  f'{round(seconds):02d} on average and {round(hours):d}:{round(minutes):02d}:{round(seconds):02d} '
                  f'in total')

    else:
        # If the first element of the result is None, the simulation didn't end, the uncompleted results file is deleted
        # and a new simulation will start over
        if args.save_data:
            os.remove(out_file)
        if not args.quiet:
            print('All individuals died or altruism went extinct in all simulations, no data generated')
            print('\033[1A', end='\x1b[2K')
            print('\033[1A', end='\x1b[2K')
            print('\033[1A', end='\x1b[2K')
        return True

    if not args.quiet:
        print(f'\x1b[2KThe results of the simulation have been saved in {out_file}')
    return False


if __name__ == '__main__':
    if not os.path.exists('output.csv'):
        with open('output.csv', 'w') as f:
            f.write('group_size,group_size_limit,population_size_limit,group_migration,survival_probability_mean,'
                    'survival_probability_sd,benefit_relatedness_exp_factor,cost_benefit_ratio,minimum_benefit,'
                    'maximum_cost,help_higher_sp_probability,help_lower_sp_probability,gained_lost_ratio,'
                    'gained_per_competition,maximum_gained,compete_higher_sp_probability,compete_lower_sp_probability,'
                    'altruism_initial_freq,altruist_percentage,last_generation\n')
    aborted = True
    # If the simulation was aborted because of lack of individuals or altruists, a new simulation will start
    while aborted:
        aborted = create_simulation_results()
