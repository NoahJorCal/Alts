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
                    group_wo_altruistic = group.copy()
                    possible_recipients = group.copy()
                    group_wo_altruistic.pop(group.index(individual))
                    possible_recipients.pop(group.index(individual))
                    #Non altruistic individuals are removed from the list
                    for possible_recipient in group_wo_altruistic[::-1]:
                        if possible_recipient.phenotype[0] != 'altruistic':
                            possible_recipients.pop(group_wo_altruistic.index(possible_recipient))
                    #If the altruistic individual is going to act altruistic, that
                    #individual is deleted from the posible benefile beneficiaries
                    #possible_recipients.pop(group.index(individual))
                    if len(possible_recipients) != 0:
                        individual.survival_probability -= altruism_cost
                        recipient = random.choice(possible_recipients)
                        recipient.survival_probability += altruism_benefit
            #If the individual has the combined allele of selfish and altruistic
            #the probability to act altruistic is half
            elif '_' in individual.phenotype[0]:
                if random.random() < altruism_probability/2:
                    group_wo_altruistic = group.copy()
                    possible_recipients = group.copy()
                    group_wo_altruistic.pop(group.index(individual))
                    possible_recipients.pop(group.index(individual))
                    for possible_recipient in group_wo_altruistic[::-1]:
                        if possible_recipient.phenotype[0] != 'altruistic':
                            possible_recipients.pop(group_wo_altruistic.index(possible_recipient))
                    if len(possible_recipients) != 0:
                        individual.survival_probability -= altruism_cost
                        recipient = random.choice(possible_recipients)
                        recipient.survival_probability += altruism_benefit


