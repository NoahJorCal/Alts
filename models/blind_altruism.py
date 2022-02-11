#!/usr/bin/python3

import random
from configparser import ConfigParser
import os
from os import path
#import alts

general_config = ConfigParser()
general_config.read('config.ini')

#model_name = os.path.basename(__file__)[:-3]
#model_config = ConfigParser()
#model_config.read(path.join('models', model_name + '.ini'))

altruism_probability = float(general_config['population']['altruism_probability'])
altruism_cost = float(general_config['population']['altruism_cost'])
altruism_benefit = float(general_config['population']['altruism_benefit'])

def selection(groups):
    survival_probabilities = []
    for group in groups:
        for individual in group:
            if individual.phenotype[0] == 'altruistic':
                if random.random() < altruism_probability:
                    possible_beneficiaries = group.copy()
                    possible_beneficiaries.pop(group.index(individual))
                    individual.survival_probability -= altruism_cost
                    benefactor = random.choice(possible_beneficiaries)
                    benefactor.survival_probability += altruism_benefit
            elif '_' in individual.phenotype[0]:
                if random.random() < altruism_probability/2:
                    possible_beneficiaries = group.copy()
                    possible_beneficiaries.pop(group.index(individual))
                    individual.survival_probability -= altruism_cost
                    benefactor = random.choice(possible_beneficiaries)
                    benefactor.survival_probability += altruism_benefit
            survival_probabilities.append(individual.survival_probability)
    #print('===========================')
    #print(survival_probabilities)
    #print([min(survival_probabilities), max(survival_probabilities)])
    #print('===========================')
    return [min(survival_probabilities), max(survival_probabilities)]
