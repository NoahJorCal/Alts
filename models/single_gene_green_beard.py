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
            if individual.phenotype[0] == 'altruist':
                if random.random() < altruism_probability:
                    ''' If the altruist individual is going to act altruist, that
                    individual is deleted from the possible beneficiaries '''
                    group_wo_altruist = group.copy()
                    group_wo_altruist.pop(group.index(individual))
                    possible_recipients = group_wo_altruist.copy()
                    # Non-altruistic individuals are removed from the list
                    if pedigree:
                        ''' If some amount of relatedness is needed to act altruistic,
                        the possible recipients are also filtered by relatedness '''
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[0] != 'altruist' or \
                                    pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    else:
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[0] != 'altruist' or \
                                    pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))

                    if len(possible_recipients) != 0:
                        recipient = random.choice(possible_recipients)
                        individual.survival_probability -= altruism_cost
                        recipient.survival_probability += altruism_benefit

            elif '_' in individual.phenotype[0]:
                ''' If the individual has the combined allele of selfish and altruist
                the probability to act altruist is half '''
                if random.random() < altruism_probability/2:
                    group_wo_altruist = group.copy()
                    possible_recipients = group.copy()
                    group_wo_altruist.pop(group.index(individual))
                    possible_recipients.pop(group.index(individual))
                    if pedigree:
                        ''' If some amount of relatedness is needed to act altruistic,
                        the possible recipients are also filtered by relatedness '''
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[0] != 'altruist' or \
                                    pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    else:
                        for possible_recipient in group_wo_altruist[::-1]:
                            if possible_recipient.phenotype[0] != 'altruist' or \
                                    pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))

                    if len(possible_recipients) != 0:
                        individual.survival_probability -= altruism_cost
                        recipient = random.choice(possible_recipients)
                        recipient.survival_probability += altruism_benefit
