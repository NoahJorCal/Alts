def relatedness_coefficient(ancestry_1, ancestry_2):
    """
    Calculated relatedness from the ancestry of two individuals.
    :param list[list[int]] ancestry_1: IDs of the ancestors of an individual. Each list within contains the IDs of the
    individuals of one tree's generation.
    :param list[list[int]] ancestry_2: IDs of the ancestors of an individual. Each list within contains the IDs of the
    individuals of one tree's generation.
    :return: Relatedness coefficient between the two individuals.
    """
    relatedness = 0

    indexes_to_skip = {0: {0: [], 1: [], 2: [], 3: []},
                       1: {0: [], 1: [], 2: [], 3: []},
                       2: {0: [], 1: [], 2: [], 3: []},
                       3: {0: [], 1: [], 2: [], 3: []}}
    # For generation in the genealogy
    for gen_1_index in range(len(ancestry_1)):
        for gen_2_index in range(len(ancestry_2)):
            # For ancestor in the generation
            for ind_1_index in range(len(ancestry_1[gen_1_index])):
                for ind_2_index in range(len(ancestry_2[gen_2_index])):
                    # If the ancestor is not unknown, it matches any other ancestor from the other ancestry and was not
                    # already considered by a closer ancestor
                    if ancestry_1[gen_1_index][ind_1_index] != 0 and \
                            ancestry_1[gen_1_index][ind_1_index] == ancestry_2[gen_2_index][ind_2_index] and \
                            (ind_1_index, ind_2_index) not in indexes_to_skip[gen_1_index][gen_2_index]:
                        # Generations between the common ancestor to both individuals
                        edges = gen_1_index + gen_2_index
                        relatedness += (1 / 2) ** edges
                        younger_gen_index = min(gen_1_index, gen_2_index)
                        gen_difference = abs(gen_1_index - gen_2_index)
                        # Adds indices to skip
                        for gen_i in range(younger_gen_index + 1, len(ancestry_1) - gen_difference):
                            if gen_1_index == younger_gen_index:
                                from_index_1 = int(ind_1_index * (len(ancestry_1[gen_i]) /
                                                                  len(ancestry_1[gen_1_index])))
                                to_index_1 = int(from_index_1 + (1 / len(ancestry_1[gen_1_index]))
                                                 * len(ancestry_1[gen_i]))
                                from_index_2 = int(ind_2_index * (len(ancestry_2[gen_i + gen_difference]) /
                                                                  len(ancestry_2[gen_2_index])))
                                to_index_2 = int(from_index_2 + (1 / len(ancestry_2[gen_2_index])) *
                                                 len(ancestry_2[gen_i + gen_difference]))
                                indexes_to_skip[gen_i][gen_i + gen_difference].extend(
                                    list(zip(range(from_index_1, to_index_1),
                                             range(from_index_2, to_index_2))))
                            else:
                                from_index_1 = int(ind_1_index * (len(ancestry_1[gen_i + gen_difference]) /
                                                                  len(ancestry_1[gen_1_index])))
                                to_index_1 = int(from_index_1 + (1 / len(ancestry_1[gen_1_index])) *
                                                 len(ancestry_1[gen_i + gen_difference]))
                                from_index_2 = int(ind_2_index * (len(ancestry_2[gen_i]) /
                                                                  len(ancestry_2[gen_2_index])))
                                to_index_2 = int(from_index_2 + (1 / len(ancestry_2[gen_2_index])) *
                                                 len(ancestry_2[gen_i]))
                                indexes_to_skip[gen_i + gen_difference][gen_i].extend(
                                    list(zip(range(from_index_1, to_index_1),
                                             range(from_index_2, to_index_2))))
    return relatedness
