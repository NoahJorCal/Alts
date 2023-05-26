#!/usr/bin/python3
from math import exp
import random
from configparser import ConfigParser
import sys
from os import path
from relatedness import relatedness_coefficient

# Add the alts directory to the Python path
sys.path.append(path.join(path.dirname(__file__), '..'))

# Import general configuration
general_config = ConfigParser()
config_path = path.join(path.dirname(__file__), '..', 'config.ini')
general_config.read(config_path)

# Save altruism configuration
exp_factor_config = float(general_config['population']['benefit_relatedness_exp_factor'])
cost_benefit_ratio_config = float(general_config['population']['cost_benefit_ratio'])
minimum_benefit_config = float(general_config['population']['minimum_benefit'])
maximum_cost_config = float(general_config['population']['maximum_cost'])


def exp_f(x, factor, penetrance, min_benefit):
    """
    Exponential distribution for calculating benefit (y) based on relatedness (x)
    :param float x: Relatedness
    :param float factor: Factor that determines the slope of the exponential distribution. The higher the value,
    the steeper slope.
    :param float penetrance: Effect of altruism.
    In homozygous individuals the penetrance is 1 and in heterozygous it is 0.5.
    :param float min_benefit: Minimum benefit altruists will give regardless of the relatedness.
    :return: Benefit for the recipient based on the relatedness.
    """
    # When the factor is 0, the distribution is lineal
    if factor == 0:
        if x * penetrance < min_benefit:
            return min_benefit
        else:
            return x * penetrance
    y = exp(factor * x) - 1
    minimum = exp(factor * 0) - 1
    maximum = exp(factor * 1) - 1
    # Distribution scaled between 0 and 1
    y = (y - minimum) / (maximum - minimum) * penetrance
    # If the calculated benefit from the relatedness is less than the minimum benefit,
    # the later will be the final benefit
    if y < min_benefit:
        return min_benefit
    else:
        return y


def get_possible_recipients(altruist, group):
    """
    Discards incompatible recipients for altruism. In this model the only one is the altruist itself.
    :param simulator.Individual altruist: Altruistic individual that will be deleted from the group.
    :param list[simulator.Individual] group: Groups where the altruist belongs to
    :return: Groups without the incompatible recipients.
    """
    possible_recipients = group.copy()
    possible_recipients.pop(group.index(altruist))
    return possible_recipients


def altruistic_act(altruist, possible_recipients, penetrance, added_cost,
                   exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost):
    """
    The altruistic individual gives part of its survival probability to a random individual from the possible recipients
    based their relatedness and an exponential distribution.
    :param simulator.Individual altruist: Altruistic individual that will give part of its survival probability.
    :param list[simulator.Individual] possible_recipients: Altruist's group without it from which the recipient will
    be selected
    :param float penetrance: Benefit based on relatedness, minimum benefit and maximum cost is based on penetrance.
    In homozygous individuals the penetrance is 1 and in heterozygous it is 0.5.
    :param float added_cost: Added survival probability sacrificed by the current altruist.
    :param float exp_factor: Factor that determines the slope of the exponential distribution. The higher the value,
    the steeper slope.
    :param float cost_benefit_ratio: Cost = Benefit * Ratio
    :param float minimum_benefit: Minimum benefit altruists will give regardless of the relatedness.
    :param float maximum_cost: When the added cost reaches the maximum cost the benefit in the interaction will be the
    remaining cost to reach the maximum.
    :return: The new modified added cost.
    """
    # The benefit is multiplied by the penetrance
    # If there is no codominance the penetrance is 1, if there is, it is 0.5
    minimum_benefit = minimum_benefit * penetrance
    maximum_cost = maximum_cost * penetrance
    # A random individual from the group is selected and deleted from the group
    random.shuffle(possible_recipients)
    recipient = possible_recipients.pop()
    # The altruist will only help the selected individual if it has equal or less survival probability
    if recipient.survival_probability <= altruist.survival_probability:
        # The ancestry list is prepared for the calculation of the relatedness
        recipient_ancestry = [[recipient.id]]
        recipient_ancestry.extend(recipient.ancestry)
        altruist_ancestry = [[altruist.id]]
        altruist_ancestry.extend(altruist.ancestry)
        relatedness_coef = relatedness_coefficient(recipient_ancestry, altruist_ancestry)
        # If the two individuals are not related, the benefit will be the base number
        if relatedness_coef == 0:
            benefit = minimum_benefit
        else:
            benefit = exp_f(relatedness_coef, exp_factor, penetrance, minimum_benefit)
        # The cost is calculated from the benefit and the survival probabilities are recalculated
        cost = benefit * cost_benefit_ratio
        # If the cost exceeds the maximum cost the altruist will sacrifice up to the maximum,
        # instead of based on the relatedness
        if maximum_cost - added_cost < cost:
            cost = maximum_cost - added_cost
            benefit = cost / cost_benefit_ratio
        # The survival probabilities are modified
        recipient.survival_probability += benefit
        altruist.survival_probability -= cost
        return added_cost + cost
    else:
        return added_cost


def selection(groups,
              exp_factor=exp_factor_config, cost_benefit_ratio=cost_benefit_ratio_config,
              minimum_benefit=minimum_benefit_config, maximum_cost=maximum_cost_config):
    """
    For each individual in the group, if it is altruistic, it will help other individuals.
    :param list[list[simulator.Individual]] groups: Whole population structured in groups.
    :param float exp_factor: Factor that determines the slope of the exponential distribution. The higher the value,
    the steeper slope.
    :param float cost_benefit_ratio: Cost = Benefit * Ratio
    :param float minimum_benefit: Minimum benefit altruists will give regardless of the relatedness.
    :param float maximum_cost: When the added cost reaches the maximum cost the benefit in the interaction will be the
    remaining cost to reach the maximum.
    """
    for group in groups:
        for individual in group:
            # If the individual is altruistic homozygous
            if individual.phenotype[0] == 'altruistic':
                possible_recipients = get_possible_recipients(individual, group)
                if len(possible_recipients) != 0:
                    added_cost = 0
                    while maximum_cost > added_cost:
                        # If the new added cost ends up being higher than the maximum cost,
                        # the last interaction will not be reverted, so in some cases,
                        # the individual will have sacrificed more fitness than expected
                        added_cost = altruistic_act(individual, possible_recipients, 1, added_cost,
                                                    exp_factor, cost_benefit_ratio,
                                                    minimum_benefit, maximum_cost)
                        if len(possible_recipients) == 0:
                            break
            # If the individual is altruistic heterozygous
            elif '_' in individual.phenotype[0]:
                possible_recipients = get_possible_recipients(individual, group)
                if len(possible_recipients) != 0:
                    added_cost = 0
                    while (maximum_cost * 0.5) > added_cost:
                        # If the new added cost ends up being higher than the maximum cost,
                        # the last interaction will not be reverted, so in some cases,
                        # the individual will have sacrificed more fitness than expected
                        added_cost = altruistic_act(individual, possible_recipients, 0.5, added_cost,
                                                    exp_factor, cost_benefit_ratio,
                                                    minimum_benefit, maximum_cost)
                        if len(possible_recipients) == 0:
                            break
