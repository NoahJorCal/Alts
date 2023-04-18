def relatedness_coefficient(ancestry_1, ancestry_2):
    relatedness = 0
    indexes_to_skip = {0: {0: [], 1: [], 2: [], 3: []},
                       1: {0: [], 1: [], 2: [], 3: []},
                       2: {0: [], 1: [], 2: [], 3: []},
                       3: {0: [], 1: [], 2: [], 3: []}}
    for gen_1_index in range(len(ancestry_1)):
        for gen_2_index in range(len(ancestry_2)):
            for ind_1_index in range(len(ancestry_1[gen_1_index])):
                for ind_2_index in range(len(ancestry_2[gen_2_index])):
                    # temp1 = ancestry_1[gen_1_index][ind_1_index]
                    # temp2 = ancestry_2[gen_2_index][ind_2_index]
                    # temp1000 = indexes_to_skip[gen_1_index]
                    # iff = (ind_1_index, ind_2_index)
                    # inn = indexes_to_skip[gen_1_index][gen_2_index]
                    if ancestry_1[gen_1_index][ind_1_index] != 0 and \
                            ancestry_1[gen_1_index][ind_1_index] == ancestry_2[gen_2_index][ind_2_index] and \
                            (ind_1_index, ind_2_index) not in indexes_to_skip[gen_1_index][gen_2_index]:
                        edges = gen_1_index + gen_2_index
                        relatedness += (1 / 2) ** edges
                        younger_gen_index = min(gen_1_index, gen_2_index)
                        gen_difference = abs(gen_1_index - gen_2_index)
                        for gen_i in range(younger_gen_index + 1, len(ancestry_1) - gen_difference):
                            if gen_1_index == younger_gen_index:
                                from_index_1 = int(ind_1_index * (len(ancestry_1[gen_i]) / len(ancestry_1[gen_1_index])))
                                to_index_1 = int(from_index_1 + (1 / len(ancestry_1[gen_1_index])) * len(ancestry_1[gen_i]))
                                from_index_2 = int(ind_2_index * (len(ancestry_2[gen_i + gen_difference]) / len(ancestry_2[gen_2_index])))
                                to_index_2 = int(from_index_2 + (1 / len(ancestry_2[gen_2_index])) * len(ancestry_2[gen_i + gen_difference]))
                                indexes_to_skip[gen_i][gen_i + gen_difference].extend(
                                    list(zip(range(from_index_1, to_index_1),
                                             range(from_index_2, to_index_2))))
                            else:
                                from_index_1 = int(ind_1_index * (len(ancestry_1[gen_i + gen_difference]) / len(ancestry_1[gen_1_index])))
                                to_index_1 = int(from_index_1 + (1 / len(ancestry_1[gen_1_index])) * len(ancestry_1[gen_i + gen_difference]))
                                from_index_2 = int(ind_2_index * (len(ancestry_2[gen_i]) / len(ancestry_2[gen_2_index])))
                                to_index_2 = int(from_index_2 + (1 / len(ancestry_2[gen_2_index])) * len(ancestry_2[gen_i]))
                                indexes_to_skip[gen_i + gen_difference][gen_i].extend(
                                    list(zip(range(from_index_1, to_index_1),
                                             range(from_index_2, to_index_2))))
    return relatedness


# print(relatedness_coefficient([[14], [10, 7], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],
#                               [[24], [10, 14], [0, 0, 10, 7], [0, 0, 0, 0, 0, 0, 0, 0]]))

# relatedness_coefficient([[2000], [1000, 1], [1, 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]],
#                         [[1000], [1, 2], [0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]])

