#!/usr/bin/python3

import random
from configparser import ConfigParser
import os
from os import path

general_config = ConfigParser()
general_config.read('config.ini')

altruism_probability = float(general_config['population']['altruism_probability'])
altruism_cost = float(general_config['population']['altruism_cost'])
altruism_benefit = float(general_config['population']['altruism_benefit'])

def selection(groups):
    survival_probabilities = []
    for group in groups:
        for individual in group:
            if individual.phenotype[0] == 'altruistic':
                if random.random() < altruism_probability:
                    possible_recipients = group.copy()
                    #If the altruistic individual is going to act altruistic, that
                    #individual is deleted from the posible benefile beneficiaries
                    possible_recipients.pop(group.index(individual))
                    individual.survival_probability -= altruism_cost
                    recipient = random.choice(possible_recipients)
                    recipient.survival_probability += altruism_benefit
            #If the individual has the combined allele of selfish and altruistic
            #the probability to act altruistic is half
            elif '_' in individual.phenotype[0]:
                if random.random() < altruism_probability/2:
                    possible_recipients = group.copy()
                    possible_recipients.pop(group.index(individual))
                    individual.survival_probability -= altruism_cost
                    recipient = random.choice(possible_recipients)
                    recipient.survival_probability += altruism_benefit
            #survival_probabilities.append(individual.survival_probability)
    #return [min(survival_probabilities), max(survival_probabilities)]
