#!/usr/bin/python3
import copy
import random
from configparser import ConfigParser

import numpy as np
from math import sqrt
import networkx as nx

from pedigree import recursive_relatedness

general_config = ConfigParser()
general_config.read('config.ini')

altruism_probability = float(general_config['population']['altruism_probability'])
relatedness = float(general_config['population']['relatedness'])
altruism_cost = float(general_config['population']['altruism_cost'])
altruism_benefit = float(general_config['population']['altruism_benefit'])


def selection(groups, pedigree=None):
    for group in groups:
        for individual in group:
            # print('BEGINNING', individual.id, individual.ancestry)
            beg = str(individual.id) + str(individual.ancestry)
            if individual.phenotype[0] == 'altruistic':
                if random.random() < altruism_probability:
                    ''' If the altruistic individual is going to act altruistic, that
                    individual is deleted from the possible  beneficiaries '''
                    group_wo_altruist = group.copy()
                    group_wo_altruist.pop(group.index(individual))
                    possible_recipients = group_wo_altruist.copy()
                    # Individuals with less relatedness with the altruist than the minimum required are deleted
                    if pedigree:
                        for possible_recipient in group_wo_altruist[::-1]:
                            # print('POSSIBLE RECIPIENTS:')
                            # print('FIRST', possible_recipient.id, possible_recipient.ancestry)

                            ancestry_1 = [[possible_recipient.id]]
                            ancestry_1.extend(possible_recipient.ancestry)
                            ancestry_2 = [[individual.id]]
                            ancestry_2.extend(individual.ancestry)
                            temp_a = ancestry_1
                            temp_b = ancestry_2
                            matrix = pedigree.relatedness([possible_recipient.id, individual.id])
                            ped = recursive_relatedness(ancestry_1, ancestry_2, 0)
                            if matrix != ped and matrix - ped > 0.1:
                                print(f'----------- Altruist {individual.id} with {possible_recipient.id} -----------')
                                print(f'{temp_a} with age {possible_recipient.age}')
                                print(f'{temp_b} with age {possible_recipient.age}')
                                # recursive_relatedness(ancestry_1, ancestry_2)
                                print(f'Matrix:   {matrix}')
                                print(f'Pedigree: {recursive_relatedness(ancestry_1, ancestry_2, 1)}')
                                # print(pedigree.removed)
                                # with open('genetic_relationship_matrix.txt', 'w') as f:
                                #     f.write(str(np.round(pedigree.kinship.todense(), 3)))

                                np.savetxt('genetic_relationship_matrix.csv', np.round(pedigree.kinship.todense(), 3), delimiter=",")
                            if pedigree.relatedness([possible_recipient.id, individual.id]) > relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                            # print('THEN ', possible_recipient.id, possible_recipient.ancestry)
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
                    if pedigree:
                        for possible_recipient in group_wo_altruist[::-1]:
                            if pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    individual.survival_probability -= altruism_cost
                    recipient = random.choice(possible_recipients)
                    recipient.survival_probability += altruism_benefit

            # print('END      ', individual.id, individual.ancestry)
            end = str(individual.id) + str(individual.ancestry)
            if beg != end:
                print('BEGINNING', beg)
                print('END      ', end)
                raise Exception('The ancestry has changed')
