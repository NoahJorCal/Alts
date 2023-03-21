#!/usr/bin/python3

import random
from configparser import ConfigParser

import numpy as np
from math import sqrt

general_config = ConfigParser()
general_config.read('config.ini')

altruism_probability = float(general_config['population']['altruism_probability'])
relatedness = float(general_config['population']['relatedness'])
altruism_cost = float(general_config['population']['altruism_cost'])
altruism_benefit = float(general_config['population']['altruism_benefit'])


# Relatedness between two individuals without taking into account inbreeding
def inbreeding(ancestry_1, ancestry_2, temp):
    common = []
    relatedness = 0
    if temp == 0:
        for l1_index in range(len(ancestry_1)):
            for l2_index in range(len(ancestry_2)):
                for g1_index in range(len(ancestry_1[l1_index])):
                    for g2_index in range(len(ancestry_2[l2_index])):
                        if ancestry_1[l1_index][g1_index] != 0 and \
                                ancestry_1[l1_index][g1_index] == ancestry_2[l2_index][g2_index]:
                            common.append(ancestry_1[l1_index][g1_index])
                            edges = l1_index + l2_index
                            relatedness += (1 / 2) ** edges
                            for l_i in range(l1_index + 1, len(ancestry_1)):
                                from_index = int((g1_index / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
                                to_index = int(from_index + (1 / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
                                ancestry_1[l_i] = np.delete(ancestry_1[l_i], range(from_index, to_index))
                            for l_i in range(l2_index + 1, len(ancestry_2)):
                                from_index = int((g2_index / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
                                to_index = int(from_index + (1 / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
                                ancestry_2[l_i] = np.delete(ancestry_2[l_i], range(from_index, to_index))
        # print(relatedness)
        # print(float(relatedness))
        return float(relatedness)

    else:
        print('Calculating relatedness for:')
        print(str(ancestry_1)+'\n'+str(ancestry_2))
        common = []
        relatedness = 0
        for l1_index in range(len(ancestry_1)):
            for l2_index in range(len(ancestry_2)):
                for g1_index in range(len(ancestry_1[l1_index])):
                    for g2_index in range(len(ancestry_2[l2_index])):
                        if ancestry_1[l1_index][g1_index] != 0 and \
                                ancestry_1[l1_index][g1_index] == ancestry_2[l2_index][g2_index]:
                            print(
                                f'From the first ind, in generation {l1_index + 1} individual {g1_index + 1} is {ancestry_1[l1_index][g1_index]} and '
                                f'from the second ind, in generation {l2_index + 1} individual {g2_index + 1} is {ancestry_2[l2_index][g2_index]}')
                            common.append(ancestry_1[l1_index][g1_index])
                            edges = l1_index + l2_index
                            print(f'Number of edges between them is {edges}')
                            relatedness += (1 / 2) ** edges
                            for l_i in range(l1_index + 1, len(ancestry_1)):
                                from_index = int((g1_index / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
                                to_index = int(from_index + (1 / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
                                ancestry_1[l_i] = np.delete(ancestry_1[l_i], range(from_index, to_index))
                            for l_i in range(l2_index + 1, len(ancestry_2)):
                                from_index = int((g2_index / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
                                to_index = int(from_index + (1 / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
                                ancestry_2[l_i] = np.delete(ancestry_2[l_i], range(from_index, to_index))
        string = str(ancestry_1) + '\n' + str(ancestry_2)
        print('Common components are:', common)
        print(string)
        print('Relatedess:', relatedness)


def calc_kinship(ind_1, ind_2, temp):
    inbreedings = []
    for individual in (ind_1, ind_2):
        ancestry_1 = []
        ancestry_2 = []
        for i in range(len(individual.ancestry)):
            split = 2 ** i
            ancestry_1.append(individual.ancestry[i][:split])
            ancestry_2.append(individual.ancestry[i][split:])
        # print(individual.ancestry)
        # print(ancestry_1)
        # print(ancestry_2)
        inbreedings.append(inbreeding(ancestry_1, ancestry_2, temp))
    # print(inbreeding)

    ancestry_1 = [[ind_1.id]]
    ancestry_1.extend(ind_1.ancestry)
    ancestry_2 = [[ind_2.id]]
    ancestry_2.extend(ind_2.ancestry)
    inbreedings.append(inbreeding(ancestry_1, ancestry_2, temp))

    # print(f'inbreedings before change:', inbreedings)
    inbreedings = [0 if elem is None else elem for elem in inbreedings]
    # print(f'inbreedings after change:', inbreedings)
    ''' rsd = 2fo / √(1 + fs)(1 + fd) where fo is the inbreeding coefficient of the descendants of the two individuals, 
    and fs and fd are the inbreeding coefficients of the two individuals (sire and dam)'''
    if temp == 0:
        relatedness = inbreedings[2] / sqrt((1 + inbreedings[0]) * (1 + inbreedings[1]))
        # print(f'Calculated inbreedings: {inbreedings[2]} / sqrt((1 + {inbreedings[0]}) * (1 + {inbreedings[1]})) = {relatedness}', inbreedings)
        return relatedness
    else:
        print(f'inbreedings: {inbreedings[2]} / sqrt((1 + {inbreedings[0]}) * (1 + {inbreedings[1]}))', inbreedings)


def selection(groups, pedigree = None):
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
                    if pedigree:
                        for possible_recipient in group_wo_altruist[::-1]:
                            if pedigree.relatedness([possible_recipient.id, individual.id]) != calc_kinship(possible_recipient, individual, 0) and pedigree.relatedness([possible_recipient.id, individual.id]) == 0.25:
                                print(f'-------------- Altruist {individual.id} with {possible_recipient.id} --------------')
                                print(possible_recipient.id, possible_recipient.ancestry)
                                print(individual.id, individual.ancestry)
                                calc_kinship(possible_recipient, individual, 1)
                                print(f'Matrix:   {pedigree.relatedness([possible_recipient.id, individual.id])}')
                                print(f'Pedigree: {calc_kinship(possible_recipient, individual, 0)}')
                            if pedigree.relatedness([possible_recipient.id, individual.id]) > relatedness:
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
                    if pedigree:
                        for possible_recipient in group_wo_altruist[::-1]:
                            if pedigree.relatedness([possible_recipient.id, individual.id]) < relatedness:
                                possible_recipients.pop(group_wo_altruist.index(possible_recipient))
                    individual.survival_probability -= altruism_cost
                    recipient = random.choice(possible_recipients)
                    recipient.survival_probability += altruism_benefit


# def calc_kinship(ind_1, ind_2, temp):
#     if temp == 0:
#
#         ancestry_1 = [[ind_1.id]]
#         ancestry_1.extend(ind_1.ancestry)
#         ancestry_2 = [[ind_2.id]]
#         ancestry_2.extend(ind_2.ancestry)
#
#         common = []
#         kinship = 0
#
#         for l1_index in range(len(ancestry_1)):
#             for l2_index in range(len(ancestry_2)):
#                 for g1_index in range(len(ancestry_1[l1_index])):
#                     for g2_index in range(len(ancestry_2[l2_index])):
#                         # if l1_index == 1 and l2_index == 1:
#                         if ancestry_1[l1_index][g1_index] != 0 and \
#                                 ancestry_1[l1_index][g1_index] == ancestry_2[l2_index][g2_index]:
#                             # print('AAAAAAAAAAAAAAAA', ancestry_1[l1_index][g1_index])
#                             # print(ancestry_2[l2_index])
#                             common.append(ancestry_1[l1_index][g1_index])
#                             edges = l1_index + l2_index
#                             # print(f'In generation {l1_index} of first and generation {l1_index} from second, '
#                             # f'a common ancestor is found ({i}). They are connected through it by {edges} edges')
#                             kinship += (1 / 2) ** edges
#                             # print(f'Kinship is {kinship} up to this point')'a',
#                             # print(ancestry_1)
#                             # print(ancestry_2)
#                             # print(ancestry_1[l1_index])
#                             # print(ancestry_2[l2_index])
#                             # print('EMPIEZA EL PRIMER FOR')
#                             # print(f'l1_index + 1: {l1_index + 1}')
#                             # print(f'g1_index: {g1_index}')
#                             # print(f'len(ancestry_1): {len(ancestry_1)}')
#                             for l_i in range(l1_index + 1, len(ancestry_1)):
#                                 # print(l_i)
#                                 # print('Before cut', ancestry_1[l_i])
#                                 # print(l1_index, len(ancestry_1[l1_index]))
#                                 from_index = int((g1_index / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
#                                 to_index = int(from_index + (1 / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
#                                 # print(from_index, to_index)
#                                 # print(ancestry_1[l_i])
#                                 # print('slice', ancestry_1[l_i][from_index: to_index])
#                                 ancestry_1[l_i] = np.delete(ancestry_1[l_i], range(from_index, to_index))
#                                 # print('After cut', ancestry_1[l_i])
#                             # print('EMPIEZA EL SEGUNDO FOR')
#                             # print(f'l2_index + 1: {l2_index + 1}')
#                             # print(f'g2_index: {g2_index}')
#                             # print(f'len(ancestry_2): {len(ancestry_2)}')
#                             for l_i in range(l2_index + 1, len(ancestry_2)):
#                                 # print(l_i)
#                                 # print('Before cut', ancestry_2[l_i])
#                                 # print(l2_index, len(ancestry_2[l2_index]))
#                                 from_index = int((g2_index / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
#                                 to_index = int(from_index + (1 / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
#                                 # print(from_index, to_index)
#                                 # print(ancestry_2[l_i])
#                                 # print('slice', ancestry_2[l_i][from_index: to_index])
#                                 ancestry_2[l_i] = np.delete(ancestry_2[l_i], range(from_index, to_index))
#                                 # print('After cut', ancestry_2[l_i])
#                             # print('-----------', ancestry_1)
#                             # print('-----------', ancestry_2)
#
#         return kinship
#
#     else:
#         ancestry_1 = [[ind_1.id]]
#         ancestry_1.extend(ind_1.ancestry)
#         ancestry_2 = [[ind_2.id]]
#         ancestry_2.extend(ind_2.ancestry)
#         common = []
#         kinship = 0
#         for l1_index in range(len(ancestry_1)):
#             for l2_index in range(len(ancestry_2)):
#                 for g1_index in range(len(ancestry_1[l1_index])):
#                     for g2_index in range(len(ancestry_2[l2_index])):
#                         if ancestry_1[l1_index][g1_index] != 0 and \
#                                 ancestry_1[l1_index][g1_index] == ancestry_2[l2_index][g2_index]:
#                             print(
#                                 f'From the first ind, in generation {l1_index + 1} individual {g1_index + 1} is {ancestry_1[l1_index][g1_index]} and '
#                                 f'from the second ind, in generation {l2_index + 1} individual {g2_index + 1} is {ancestry_2[l2_index][g2_index]}')
#                             common.append(ancestry_1[l1_index][g1_index])
#                             edges = l1_index + l2_index
#                             print(f'Number of edges between them is {edges}')
#                             kinship += (1 / 2) ** edges
#                             for l_i in range(l1_index + 1, len(ancestry_1)):
#                                 from_index = int((g1_index / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
#                                 to_index = int(from_index + (1 / len(ancestry_1[l1_index])) * len(ancestry_1[l_i]))
#                                 ancestry_1[l_i] = np.delete(ancestry_1[l_i], range(from_index, to_index))
#                             for l_i in range(l2_index + 1, len(ancestry_2)):
#                                 from_index = int((g2_index / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
#                                 to_index = int(from_index + (1 / len(ancestry_2[l2_index])) * len(ancestry_2[l_i]))
#                                 ancestry_2[l_i] = np.delete(ancestry_2[l_i], range(from_index, to_index))
#         string = str(ancestry_1) + '\n' + str(ancestry_2)
#         print(common)
#         return string