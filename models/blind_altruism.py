#!/usr/bin/python3
import copy
import random
from configparser import ConfigParser

import numpy as np
from math import sqrt
import networkx as nx

from pedigree import inbreeding

general_config = ConfigParser()
general_config.read('config.ini')

altruism_probability = float(general_config['population']['altruism_probability'])
relatedness = float(general_config['population']['relatedness'])
altruism_cost = float(general_config['population']['altruism_cost'])
altruism_benefit = float(general_config['population']['altruism_benefit'])



def inbreeding_coefficient(ancestry_1, ancestry_2, temp):
    common = []
    inbreeding_coef = 0
    if temp == 1:
        print(f'Ancestry 1 :{ancestry_1}')
        print(f'Ancestry 2 :{ancestry_2}')

    for l1_index in range(len(ancestry_1)):
        for l2_index in range(len(ancestry_2)):
            for g1_index in range(len(ancestry_1[l1_index])):
                for g2_index in range(len(ancestry_2[l2_index])):
                    if ancestry_1[l1_index][g1_index] != 0 and \
                            ancestry_1[l1_index][g1_index] == ancestry_2[l2_index][g2_index]:
                        common.append(ancestry_1[l1_index][g1_index])
                        edges = l1_index + l2_index + 1
                        inbreeding_coef += (1 / 2) ** edges

                        if temp == 1:
                            print(f'They have ancestors {common} in common')
                            print(f'There are {edges} edges for the last ancestor')
                            print(f'{(1 / 2) ** edges} is added to the inbreeding coefficient -> {inbreeding_coef}')

                        for l_i in range(l1_index + 1, len(ancestry_1)):
                            from_index = int(g1_index * (len(ancestry_1[l_i]) / len(ancestry_1[l1_index])))
                            to_index = int(from_index + (1 / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
                            for i in range(from_index, to_index):
                                ancestry_1[l_i][i] = 0
                            # ancestry_1[l_i] = np.delete(ancestry_1[l_i], range(from_index, to_index))
                        for l_i in range(l2_index + 1, len(ancestry_2)):
                            from_index = int(g2_index * (len(ancestry_2[l_i]) / len(ancestry_2[l2_index])))
                            to_index = int(from_index + (1 / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
                            for i in range(from_index, to_index):
                                ancestry_2[l_i][i] = 0
    if temp == 1:
        print(f'Relatedness is {inbreeding_coef}')
    return inbreeding_coef


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
                            ancestry_1 = [[possible_recipient.id]]
                            ancestry_1.extend(possible_recipient.ancestry)
                            ancestry_2 = [[individual.id]]
                            ancestry_2.extend(individual.ancestry)
                            rel = inbreeding_coefficient(ancestry_1, ancestry_2, 0)
                            # ped = pedigree_from_ancestry(ancestry_1, ancestry_2)
                            # matrix = pedigree.relatedness([possible_recipient.id, individual.id])
                            # small_matrix = ped.relatedness([possible_recipient.id, individual.id])
                            # if matrix != small_matrix and matrix - small_matrix > 0.1:
                            #     print(f'----------- Altruist {individual.id} with {possible_recipient.id} -----------')
                            #     print(f'{ancestry_1} with age {possible_recipient.age}')
                            #     print(f'{ancestry_2} with age {possible_recipient.age}')
                            #     # recursive_relatedness(ancestry_1, ancestry_2)
                            #     print(f'Matrix:   {matrix}')
                            #     print(f'Pedigree: {small_matrix}')
                            #     # print(pedigree.removed)
                            #     # with open('genetic_relationship_matrix.txt', 'w') as f:
                            #     #     f.write(str(np.round(pedigree.kinship.todense(), 3)))
                            #
                            #     np.savetxt('genetic_relationship_matrix.csv', np.round(pedigree.kinship.todense(), 3), delimiter=",")
                            if rel > relatedness:
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
