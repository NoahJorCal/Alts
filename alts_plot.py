from configparser import ConfigParser
from matplotlib import pyplot as plt
from alts import alts_main
import re
from itertools import product

#Import general configuration
general_config = ConfigParser()
general_config.read('config.ini')

def plot():
    #Run the main function of the alts script, runs a simulation and generates a summary
    simulation_summary, dict_phenotypes_combinations_indexes, dict_phenotype_options = alts_main()
    #print(dict_phenotype_options, dict_phenotypes_combinations_indexes)
    population_size = int(general_config['population']['size'])
    genes = general_config['plot']['genes'].replace(' ','').split(',')
    generation_x = list(range(int(general_config['simulation']['generations'])+1))
    combined_phenotypes = list(dict_phenotypes_combinations_indexes.keys())
    plot_info = {}
    phenotype_options = []
    #print('SUMMARY', simulation_summary)
    #print('GENES A GRAFICAR', genes)
    for gene in genes:
        #print('OPCIONES DE ESE GEN',dict_phenotype_options[gene])
        phenotype_options.append(dict_phenotype_options[gene])

    #print(phenotype_options)
    if len(genes) == 1:
        #print('SOLO SE GRAFICA UN GEN')
        #match_indexes = list(range(len(simulation_summary[0])))
        #for allele in phenotype_options:
            #print('EL PRIMER ALELO ES', allele)
            #new_match_indexes = []
            #for index in match_indexes:
                #print('LAS COMBINACIONES DE ALELOS SON', list(dict_phenotypes_combinations_indexes.keys()))
                #phenotype = list(dict_phenotypes_combinations_indexes.keys())[index]
                #print(phenotype)
                #if re.search('(?:&|^)('+allele+')(?:&|$)', phenotype):
                    #new_match_indexes.append(index)
                    #print('VA AQUI Y LOS INDICES SON', new_match_indexes)
        #A new list of lists is created to save the data of the selected genes
        #survivors_summary = []
        #total_summary = []
        for phenotype_index in range(len(phenotype_options[0])):
           #print('EL ALELO ES',phenotype_options[0][phenotype_index])
            survivors_summary = ([0] * (generation_x[-1] + 1))
            total_summary = ([0] * (generation_x[-1] + 1))
           #print('PLOT DATA SUMMARY', survivors_summary)
            for combination_index in range(len(combined_phenotypes)):
               #print('PHENOTYPE', combined_phenotypes[combination_index])
                if re.search('(?:&|^)('+phenotype_options[0][phenotype_index]+')(?:&|$)', combined_phenotypes[combination_index]):
                   #print(phenotype_options[0][phenotype_index], combined_phenotypes[combination_index])
                    for generation in generation_x:
                        survivors_summary[generation] += simulation_summary[0][combination_index][generation]
                        total_summary[generation] += (simulation_summary[1][combination_index][generation])/population_size
                        plot_info[phenotype_options[0][phenotype_index]] = []
                        plot_info[phenotype_options[0][phenotype_index]].append(survivors_summary)
                        plot_info[phenotype_options[0][phenotype_index]].append(total_summary)
       #print(plot_info)

        #for generation in generation_x:
            #print(generation)
            #for index in match_indexes:
                #count_list[generation] += simulation_summary[0][index][generation]/population_size
            #plot_info[phenotype_name[:-1]] = [count_list, simulation_summary[1][index]]


    else:
        for combination in product(*phenotype_options):
           #print('hasgbdnk',combination)
            match_indexes = list(range(len(simulation_summary[0])))
            phenotype_name = ''
            for element in combination:
               #print(element)
                new_match_indexes = []
                for index in match_indexes:
                    phenotype = list(dict_phenotypes_combinations_indexes.keys())[index]
                   #print(phenotype)
                    if re.search('(?:&|^)('+element+')(?:&|$)', phenotype):
                       #print('COINCIDE', element, 'CON', phenotype)
                        new_match_indexes.append(index)
                match_indexes = new_match_indexes
                phenotype_name += element+'&'
               #print('EL NOMBRE DEL FENOTIPO SERA',phenotype_name)
           #print('LOS INDICES QUE HACEN MATCH SON', match_indexes)

            survivors_summary = [0] * (generation_x[-1] + 1)
            total_summary = [0] * (generation_x[-1] + 1)
            for generation in generation_x:
                for index in match_indexes:
                    survivors_summary[generation] += simulation_summary[0][index][generation]
                    total_summary[generation] += (simulation_summary[1][index][generation])/population_size
               #print('LOS SUPERVIVIENTES SE VAN RELLENANDO',survivors_summary)
               #print('LOS TOTALES SE VAN RELLENANDO',total_summary)
            plot_info[phenotype_name[:-1]] = []
            plot_info[phenotype_name[:-1]].append(survivors_summary)
            plot_info[phenotype_name[:-1]].append(total_summary)

       #print('SUPERVIVIENTES DEL SUMMARY',simulation_summary[0])
       #print('TOTALES DEL SUMMARY',simulation_summary[1])
       #print('LA INFO DEL PLOT HA QUEDADO', plot_info)
   #print('LONGITUD DE PLOT INFO',len(plot_info))
    total_survivors = []
    #for generation in generation_x:
        #individuals_sum = 0
        #for phenotype in range(len(survivors_summary)):
            #individuals_sum += survivors_summary[phenotype][generation]
        #total_survivors.append(individuals_sum)

    for generation in generation_x:
        individuals_sum = 0
        for phenotype_index in range(len(plot_info)):
            individuals_sum += plot_info[list(plot_info.keys())[phenotype_index]][0][generation]
        total_survivors.append(individuals_sum)
   #print('TOTAL SURVIVORS', total_survivors)
   #print('RESUMEN',survivors_summary, '\t',total_summary)

    legend_population_size = 'N = '+str(population_size)
    #print(list(plot_info.values()))
    proportions_data_y = [values[1] for values in list(plot_info.values())]
    survivors_data_y = [values[0] for values in list(plot_info.values())]
    survivors_data_y.append(total_survivors)
   #print('PROPORCIONES',proportions_data_y)
   #print('SUPERVIVIENTES',survivors_data_y)
    legend_phenotypes = list(plot_info.keys())
   #print('FENOTIPOS PARA LA LEYENDA', legend_phenotypes)
   #print('PRINCIPIO 1er plot')
    proportions_plot = plt.figure(1)
    plt.stackplot(generation_x, proportions_data_y, labels = legend_phenotypes)
    plt.margins(0)
    plt.title('Proportion of individuals by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Proportion of individuals')
    plt.legend()

   #print('FIN 1er plot')

    survivors_plot = plt.figure(2)


    #print('PRINCIPIO 2do plot')
    legend_phenotypes.append('total')
    for count in range(len(survivors_data_y)):
        plt.plot(survivors_data_y[count], label = legend_phenotypes[count])

    #plt.plot([], [], ' ', label = legend_population_size)

    plt.margins(0)
    plt.ylim((0,population_size))
    plt.title('Survivors per generation by phenotype')
    plt.xlabel('Generation')
    plt.ylabel('Number of individuals')
    plt.legend()

   #print('FIN 2do plot')
    plt.show()

plot()
