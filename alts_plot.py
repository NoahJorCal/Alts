from configparser import ConfigParser
from os import path
from matplotlib import pyplot as plt

def plot(simulation_summary):
    generation_x = range(len(list(simulation_summary.values())[0]))
    phenotypes_y = []
    legend_phenotypes = []
    for phenotype, count in simulation_summary.items():
        phenotypes_y.append(count)
        legend_phenotypes.append(phenotype)
    for phenotype in phenotypes_y:
        plt.plot(generation_x, phenotype)

    plt.title('Number of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Number of individuals')
    plt.legend(legend_phenotypes)

    plt.show()













