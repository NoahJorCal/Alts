#!/usr/bin/python3
import random
from configparser import ConfigParser

from relatedness import relatedness_coefficient

general_config = ConfigParser()
general_config.read('config.ini')

altruism_probability = float(general_config['population']['altruism_probability'])
relatedness = float(general_config['population']['relatedness'])
altruism_cost = float(general_config['population']['altruism_cost'])
altruism_benefit = float(general_config['population']['altruism_benefit'])


def selection(groups):
    for group in groups:
        for individual in group:
            if individual.phenotype[0] == 'altruistic':
                if random.random() < altruism_probability:
                    ''' If the altruistic individual is going to act altruistic, that
                    individual is deleted from the possible  beneficiaries '''
                    group_wo_altruist = group.copy()
                    group_wo_altruist.pop(group.index(individual))
                    possible_recipients = group_wo_altruist.copy()
                    # Individuals with less relatedness with the altruist than the minimum required are deleted
                    if relatedness > 0:
                        ''' If some amount of relatedness is needed to act altruistic,
                        the possible recipients are also filtered by relatedness '''
                        for possible_recipient in group_wo_altruist[::-1]:
                            ancestry_1 = [[possible_recipient.id]]
                            ancestry_1.extend(possible_recipient.ancestry)
                            ancestry_2 = [[individual.id]]
                            ancestry_2.extend(individual.ancestry)
                            relatedness_coef = relatedness_coefficient(ancestry_1, ancestry_2)
                            if relatedness_coef > relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    individual.survival_probability -= altruism_cost
                    recipient = random.choice(possible_recipients)
                    recipient.survival_probability += altruism_benefit
            elif '_' in individual.phenotype[0]:
                ''' If the individual has the combined allele of selfish and altruistic,
                the probability to act altruistic is half '''
                if random.random() < altruism_probability/2:
                    group_wo_altruist = group.copy()
                    group_wo_altruist.pop(group.index(individual))
                    possible_recipients = group_wo_altruist.copy()
                    if relatedness > 0:
                        ''' If some amount of relatedness is needed to act altruistic,
                        the possible recipients are also filtered by relatedness '''
                        for possible_recipient in group_wo_altruist[::-1]:
                            ancestry_1 = [[possible_recipient.id]]
                            ancestry_1.extend(possible_recipient.ancestry)
                            ancestry_2 = [[individual.id]]
                            ancestry_2.extend(individual.ancestry)
                            relatedness_coef = relatedness_coefficient(ancestry_1, ancestry_2)
                            if relatedness_coef > relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    individual.survival_probability -= altruism_cost
                    recipient = random.choice(possible_recipients)
                    recipient.survival_probability += altruism_benefit
