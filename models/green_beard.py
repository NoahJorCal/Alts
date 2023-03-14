#!/usr/bin/python3

import random
from configparser import ConfigParser

general_config = ConfigParser()
general_config.read('config.ini')

altruism_probability = float(general_config['population']['altruism_probability'])
relatedness = float(general_config['population']['relatedness'])
altruism_cost = float(general_config['population']['altruism_cost'])
altruism_benefit = float(general_config['population']['altruism_benefit'])


def selection(groups, pedigree = None):
    for group in groups:
        for individual in group:
            if individual.phenotype[0] == 'altruistic':
                if random.random() < altruism_probability:
                    ''' If the altruistic individual is going to act altruistic, that
                    individual is deleted from the possible beneficiaries '''
                    group_wo_altruist = group.copy()
                    group_wo_altruist.pop(group.index(individual))
                    possible_recipients = group_wo_altruist.copy()
                    # Non-altruistic individuals are removed from the list
                    if pedigree:
                        ''' If some amount of relatedness is needed to act altruistic,
                        the possible recipients are also filtered by relatedness '''
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[1] != 'beard' or \
                                    pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    else:
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[1] != 'beard':
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))

                    if len(possible_recipients) != 0:
                        individual.survival_probability -= altruism_cost
                        recipient = random.choice(possible_recipients)
                        recipient.survival_probability += altruism_benefit
            elif '_' in individual.phenotype[0]:
                ''' If the individual has the combined allele of selfish and altruistic,
                the probability to act altruistic is half '''
                if random.random() < altruism_probability:
                    group_wo_altruist = group.copy()
                    group_wo_altruist.pop(group.index(individual))
                    possible_recipients = group_wo_altruist.copy()
                    if pedigree:
                        ''' If some amount of relatedness is needed to act altruistic,
                        the possible recipients are also filtered by relatedness '''
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[1] != 'beard' or \
                                    pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    else:
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[1] != 'beard':
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))

                    if len(possible_recipients) != 0:
                        individual.survival_probability -= altruism_cost
                        recipient = random.choice(possible_recipients)
                        recipient.survival_probability += altruism_benefit
