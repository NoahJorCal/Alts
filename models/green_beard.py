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
                    group_wo_altruistic.pop(group.index(individual))
                    possible_recipients = group_wo_altruistic.copy()
                    # Non altruistic individuals are removed from the list
                    for possible_recipient in group_wo_altruistic[::-1]:
                        if possible_recipient.phenotype[1] != 'beard':
                            possible_recipients.pop(group_wo_altruistic.index(possible_recipient))
                    ''' If the altruistic individual is going to act altruistic, that
                    individual is deleted from the posible benefile beneficiaries '''
                    if len(possible_recipients) != 0:
                        individual.survival_probability -= altruism_cost
                        recipient = random.choice(possible_recipients)
                        recipient.survival_probability += altruism_benefit
            elif '_' in individual.phenotype[0]:
                ''' If the individual has the combined allele of selfish and altruistic,
                the probability to act altruistic is half '''
                if random.random() < altruism_probability:
                    group_wo_altruistic = group.copy()
                    group_wo_altruistic.pop(group.index(individual))
                    possible_recipients = group_wo_altruistic.copy()
                    for possible_recipient in group_wo_altruistic[::-1]:
                        if possible_recipient.phenotype[1] != 'beard':
                            possible_recipients.pop(group_wo_altruistic.index(possible_recipient))
                    if len(possible_recipients) != 0:
                        individual.survival_probability -= altruism_cost
                        recipient = random.choice(possible_recipients)
                        recipient.survival_probability += altruism_benefit
