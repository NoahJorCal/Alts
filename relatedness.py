def relatedness_coefficient(ancestry_1, ancestry_2):
    relatedness = 0
    indexes_to_skip = {0: [], 1: [], 2: [], 3: []}
    for gen_1_index in range(len(ancestry_1)):
        for gen_2_index in range(len(ancestry_2)):
            for ind_1_index in range(len(ancestry_1[gen_1_index])):
                for ind_2_index in range(len(ancestry_2[gen_2_index])):
                    if ancestry_1[gen_1_index][ind_1_index] != 0 and \
                            ancestry_1[gen_1_index][ind_1_index] == ancestry_2[gen_2_index][ind_2_index] and \
                            (ind_1_index, ind_2_index) not in indexes_to_skip[gen_1_index]:
                        edges = gen_1_index + gen_2_index
                        relatedness += (1 / 2) ** edges

                        for gen_i in range(gen_1_index + 1, len(ancestry_1)):
                            from_index_1 = int(ind_1_index * (len(ancestry_1[gen_i]) / len(ancestry_1[gen_1_index])))
                            to_index_1 = int(from_index_1 + (1 / len(ancestry_1[gen_1_index])) * len(ancestry_1[gen_i]))
                            from_index_2 = int(ind_2_index * (len(ancestry_2[gen_i]) / len(ancestry_2[gen_2_index])))
                            to_index_2 = int(from_index_2 + (1 / len(ancestry_2[gen_2_index])) * len(ancestry_2[gen_i]))
                            indexes_to_skip[gen_i] = list(zip(range(from_index_1, to_index_1),
                                                              range(from_index_2, to_index_2)))
    return relatedness
