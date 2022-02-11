from configparser import ConfigParser
from alts import main
from matplotlib import pyplot as plt

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

def create_simulation_results():
    number_of_simulations = int(general_config['simulation']['simulations_per_summary'])
    population_size = int(general_config['population']['size'])
    generation_x = range(int(general_config['simulation']['generations'])+1)
    simulations_summary = []
    count = 0
    for i in range(number_of_simulations):
        simulation_summary = main()
        #print(list(simulation_summary.values())[0])
        simulations_summary.append(list(simulation_summary.values())[0])
        print(count)
        count += 1
    mean_values = []
    other_list = []
    #print(simulations_summary)
    for i in range(len(simulations_summary[0])):
        mean = 0
        for summary in simulations_summary:
            mean += summary[i]
        print(mean)
        mean = (mean/len(simulations_summary))/population_size
        print(mean)
        mean_values.append(mean)
        other_list.append(1-mean)
    print(mean_values)
    print(other_list)
    plt.stackplot(generation_x, mean_values, other_list)

    plt.title('Number of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Number of individuals')
    plt.legend(['selfish', 'altruism'])

    plt.show()




create_simulation_results()

