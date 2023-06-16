#!/usr/bin/python3
from math import exp
import random
from configparser import ConfigParser
import sys
from os import path
from relatedness import relatedness_coefficient
from scipy.stats import truncnorm

# # Add the alts directory to the Python path
# sys.path.append(path.join(path.dirname(__file__), '..'))
#
# # Import general configuration
# general_config = ConfigParser()
# config_path = path.join(path.dirname(__file__), '..', 'config.ini')
# general_config.read(config_path)

# # Import variable parameters configuration
# var_params_config = ConfigParser()
# var_params_path = path.join(path.dirname(__file__), '..', 'variable_config.ini')
# var_params_config.read(var_params_path)

# # Save altruism configuration
# exp_factor_config = float(general_config['population']['benefit_relatedness_exp_factor'])
# cost_benefit_ratio_config = float(general_config['population']['cost_benefit_ratio'])
# minimum_benefit_config = float(general_config['population']['minimum_benefit'])
# maximum_cost_config = float(general_config['population']['maximum_cost'])
# help_higher_sp_probability_config = float(general_config['population']['help_higher_sp_probability'])
# help_lower_sp_probability_config = float(general_config['population']['help_lower_sp_probability'])
# # Save selfish configuration
# gained_lost_ratio_config = float(general_config['population']['gained_lost_ratio'])
# gained_per_competition_config = float(general_config['population']['gained_per_competition'])
# maximum_gained_config = float(general_config['population']['maximum_gained'])
# compete_higher_sp_probability_config = float(general_config['population']['compete_higher_sp_probability'])
# compete_lower_sp_probability_config = float(general_config['population']['compete_higher_sp_probability'])


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
                   exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,
                   help_higher_sp_probability, help_lower_sp_probability):
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
    :param float help_higher_sp_probability: Probability of helping an individual with higher survival probability.
    :param float help_lower_sp_probability: Probability of helping an individual with lower survival probability.
    :return: The new modified added cost.
    """
    # The benefit is multiplied by the penetrance
    # If there is no codominance the penetrance is 1, if there is, it is 0.5
    minimum_benefit *= penetrance
    maximum_cost *= penetrance
    # A random individual from the group is selected and deleted from the group
    random.shuffle(possible_recipients)
    recipient = possible_recipients.pop()
    # Decides if it's going to help based on if the recipient has more or less sp and the configured probabilities
    if altruist.survival_probability >= recipient.survival_probability:
        act_altruistic = random.random() < help_higher_sp_probability
    else:
        act_altruistic = random.random() < help_lower_sp_probability
    if act_altruistic:
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


def selfish_act(selfish, possible_recipients, penetrance, added_gained,
                gained_lost_ratio, gained_per_competition, maximum_gained,
                compete_higher_sp_probability, compete_lower_sp_probability):
    """
    The selfish individual competes with another individual, the selfish will gain less survival probability than lost
    by the competitor as the competition has a wear effect on both individuals.
    :param simulator.Individual selfish: Selfish individual that will compete with other individuals
    :param list[simulator.Individual] possible_recipients: Altruist's group without it from which the recipient will
    be selected
    :param float penetrance: Lost based on gain, gain per interaction and maximum gain is based on penetrance.
    In homozygous individuals the penetrance is 1 and in heterozygous it is 0.5.
    :param float added_gained: Added survival probability gained by the current selfish individual.
    :param float gained_lost_ratio: Lost = Gain * Ratio
    :param float gained_per_competition: Gained survival probability by the selfish in each interaction.
    :param float maximum_gained: When the added gain reaches the maximum, the gain in the interaction will be the
    remaining gain to reach the maximum.
    :param float compete_higher_sp_probability: Probability of competing with an individual with higher
    survival probability.
    :param float compete_lower_sp_probability: Probability of competing with an individual with lower
    survival probability.
    :return: The new modified added gain.
    """
    # The gained survival probability is multiplied by the penetrance
    gained_per_competition *= penetrance
    maximum_gained *= penetrance
    random.shuffle(possible_recipients)
    competitor = possible_recipients.pop()
    # Decides if it's going to compete based on if the competitor has more or less sp and the configured probabilities
    if competitor.survival_probability >= selfish.survival_probability:
        act_altruistic = random.random() < compete_higher_sp_probability
    else:
        act_altruistic = random.random() < compete_lower_sp_probability
    if act_altruistic:
        # If the gained sp in this interaction is going to exceed the maximum, the remainder is gained instead
        if maximum_gained - added_gained < gained_per_competition:
            gained = maximum_gained - added_gained
        else:
            gained = gained_per_competition
        # The lost is calculated from the gained sp and the ratio
        lost = gained * gained_lost_ratio
        selfish.survival_probability += gained
        competitor.survival_probability -= lost
        return added_gained + gained
    else:
        return added_gained


def selection(groups, altruism_config, selfishness_config):
    """
    For each individual in the group, if it is altruistic, it will help other individuals.
    :param list[list[simulator.Individual]] groups: Whole population structured in groups.
    :param list[float] altruism_config: List with altruism configuration parameters:\n
    exp_factor: Factor that determines the slope of the exponential distribution.
    The higher the value, the steeper slope.\n
    cost_benefit_ratio: Cost = Benefit * Ratio\n
    minimum_benefit: Minimum benefit altruists will give regardless of the relatedness.\n
    maximum_cost: When the added cost reaches the maximum cost the benefit in the interaction
    will be the remaining cost to reach the maximum.\n
    help_higher_sp_probability: Probability of helping an individual with higher survival probability.\n
    help_lower_sp_probability: Probability of helping an individual with lower survival probability.\n
    :param list[float] selfishness_config: List with selfishness configuration parameters:\n
    gained_lost_ratio: Lost = Gained * Ratio\n
    gained_per_competition: Survival probability gained in each competition interaction.\n
    maximum_gained: When the added gained survival probability reaches the maximum, the lost survival probability
    in the interaction will be the remaining gain to reach the maximum.\n
    compete_higher_sp_probability: Probability of competing with an individual with higher survival probability.\n
    compete_lower_sp_probability: Probability of competing with an individual with lower survival probability.\n
    """
    exp_factor = altruism_config[0]
    cost_benefit_ratio = altruism_config[1]
    minimum_benefit = altruism_config[2]
    maximum_cost = altruism_config[3]
    help_higher_sp_probability = altruism_config[4]
    help_lower_sp_probability = altruism_config[5]
    gained_lost_ratio = selfishness_config[0]
    gained_per_competition = selfishness_config[1]
    maximum_gained = selfishness_config[2]
    compete_higher_sp_probability = selfishness_config[3]
    compete_lower_sp_probability = selfishness_config[4]
    for group in groups:
        for individual in group:
            # If the individual is altruistic homozygous
            if individual.phenotype[0] == 'altruistic':
                possible_recipients = get_possible_recipients(individual, group)
                if len(possible_recipients) != 0:
                    added_cost = 0
                    while maximum_cost > added_cost:
                        added_cost = altruistic_act(individual, possible_recipients, 1, added_cost,
                                                    exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,
                                                    help_higher_sp_probability, help_lower_sp_probability)
                        if len(possible_recipients) == 0:
                            break
            # If the individual is selfish homozygous
            elif individual.phenotype[0] == 'selfish':
                possible_recipients = get_possible_recipients(individual, group)
                if len(possible_recipients) != 0:
                    added_gained = 0
                    while maximum_gained > added_gained:
                        added_gained = selfish_act(individual, possible_recipients, 1, added_gained,
                                                   gained_lost_ratio, gained_per_competition, maximum_gained,
                                                   compete_higher_sp_probability, compete_lower_sp_probability)
                        if len(possible_recipients) == 0:
                            break
            # If the individual is heterozygous
            else:
                # Chose if the individual will act selfish or altruistically
                if random.random() < 0.5:
                    # Altruistic
                    possible_recipients = get_possible_recipients(individual, group)
                    if len(possible_recipients) != 0:
                        added_cost = 0
                        while maximum_cost > added_cost:
                            added_cost = altruistic_act(individual, possible_recipients, 0.5, added_cost,
                                                        exp_factor, cost_benefit_ratio, minimum_benefit, maximum_cost,
                                                        help_higher_sp_probability, help_lower_sp_probability)
                            if len(possible_recipients) == 0:
                                break
                else:
                    # Selfish
                    possible_recipients = get_possible_recipients(individual, group)
                    if len(possible_recipients) != 0:
                        added_gained = 0
                        while maximum_gained > added_gained:
                            added_gained = selfish_act(individual, possible_recipients, 0.5, added_gained,
                                                       gained_lost_ratio, gained_per_competition, maximum_gained,
                                                       compete_higher_sp_probability, compete_lower_sp_probability)
                            if len(possible_recipients) == 0:
                                break
