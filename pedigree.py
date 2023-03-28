from itertools import product
from math import sqrt


def ancestry_to_graph(ancestry_1, ancestry_2):
    relatedness_inds = (ancestry_1[0][0], ancestry_2[0][0])
    coef_relatedness = 0
    pedigree = {}
    for ancestry in (ancestry_1, ancestry_2):
        for generation_index in range(len(ancestry) - 2, -1, -1):
            for individual_index in range(len(ancestry[generation_index])):
                from_index = int(individual_index * (len(ancestry[generation_index + 1]) / len(ancestry[generation_index])))
                to_index = int(from_index + (1 / len(ancestry_1[generation_index])) * len(ancestry_1[generation_index + 1]))
                parents = ancestry[generation_index + 1][from_index:to_index]
                for parent in parents:
                    if parent in relatedness_inds:
                        coef_relatedness += 0.5
                    elif parent != 0:
                        ind = ancestry[generation_index][individual_index]
                        if parent not in pedigree.keys():
                            pedigree[parent] = [ind]
                        else:
                            if ind not in pedigree[parent]:
                                pedigree[parent].append(ind)
    return pedigree, coef_relatedness, relatedness_inds


def flat_list(nested_list, output):
    for i in nested_list:
        if type(i) == list:
            flat_list(i, output)
        else:
            output.append(i)
    return output


def split_list(lst, val):
    result = []
    current_sublist = []
    for i, num in enumerate(lst):
        current_sublist.append(num)
        if num == val:
            result.append(current_sublist[:-1])
            current_sublist = [val]
    result.append(current_sublist)
    result = [x for x in result if x]
    return result


def clean_list(input_list, staring_node, end_nodes):
    flat_routes = flat_list(input_list, [])

    seen = set()
    repeated = []
    for x in flat_routes:
        if x in seen:
            repeated.append(x)
        else:
            seen.add(x)
    contains_end_nodes = all(item in flat_routes for item in end_nodes)
    excluded = list(end_nodes)
    excluded.append(staring_node)
    repeated_ancestors = len(set(repeated).difference(excluded))
    if contains_end_nodes and repeated_ancestors == 0:
        cleaned_routes = split_list(flat_routes, staring_node)
        return cleaned_routes


def every_dfs(graph, start_node, end_nodes):
    if start_node in end_nodes:
        return [start_node]
    else:
        route = []
        for neighbour in graph[start_node]:
            route.append([start_node])
            route[-1].extend(every_dfs(graph, neighbour, end_nodes))
    return route


# ped = {10: [13], 6: [13, 21, 12], 13: [21], 1: [12], 12: [18], 2: [18]}
# ped = {1: [6, 7], 2: [6, 7], 3: [8], 4: [8, 9], 5: [9], 6: [10], 7: [11], 8: [10], 9: [11]}
# ped = {1: [4], 2: [4], 3: [5], 4: [5, 6, 7], 5: [6], 6: [7]}
# relatedness_individuals = [6, 7]


def inbreeding(ancestry_1, ancestry_2, temp):
    pedigree, coef_relatedness, relatedness_inds = ancestry_to_graph(ancestry_1, ancestry_2)
    # print(pedigree)
    routes = []
    for node in pedigree.keys():
        # print(node)
        nested_routes = every_dfs(pedigree, node, relatedness_inds)
        node_routes = clean_list(nested_routes, node, relatedness_inds)
        if node_routes:
            routes.append(node_routes)
    # print(routes)

    for ancestor_routes in routes:
        # print(ancestor_routes)
        ind_paths = []
        for ind in relatedness_inds:
            path = []
            for route in ancestor_routes:
                if ind in route:
                    path.append(route)
            ind_paths.append(path)
            # ind_paths.append([route for route in ancestor_routes if ind in route])
        for connection in product(ind_paths[0], ind_paths[1]):
            # fo = Σ 0.5^(n + n' + 1)(1 + fa)
            edges = len(connection[0]) + len(connection[1]) - 1
            coef_relatedness += 0.5 ** edges

    if relatedness_inds[0] in pedigree.keys() and relatedness_inds[1] in pedigree[relatedness_inds[0]]:
        coef_relatedness += 0.5
    elif relatedness_inds[1] in pedigree.keys() and relatedness_inds[0] in pedigree[relatedness_inds[1]]:
        coef_relatedness += 0.5

    return coef_relatedness


def recursive_relatedness(ancestry_1, ancestry_2, temp):
    if len(ancestry_1) == 2:
        inbreeding_coef = inbreeding(ancestry_1, ancestry_2, temp)
        return inbreeding_coef
    else:
        sire_ancestry_1 = []
        sire_ancestry_2 = []
        for i in range(1, len(ancestry_1)):
            split = 2 ** (i - 1)
            sire_ancestry_1.append(ancestry_1[i][:split])
            sire_ancestry_2.append(ancestry_1[i][split:])
        sire_inbreeding = recursive_relatedness(sire_ancestry_1, sire_ancestry_2, temp)
        dam_ancestry_1 = []
        dam_ancestry_2 = []
        for i in range(1, len(ancestry_2)):
            split = 2 ** (i - 1)
            dam_ancestry_1.append(ancestry_2[i][:split])
            dam_ancestry_2.append(ancestry_2[i][split:])
        dam_inbreeding = recursive_relatedness(dam_ancestry_1, dam_ancestry_2, temp)
        desc_inbreeding = inbreeding(ancestry_1, ancestry_2, temp)
        ''' rsd = 2fo / √(1 + fs)(1 + fd) where fo is the inbreeding coefficient of the descendants of the sire and the
        dam, and fs and fd are the inbreeding coefficients of the two individuals (sire and dam)'''
        relatedness = 2 * desc_inbreeding / sqrt((1 + sire_inbreeding) * (1 + dam_inbreeding))
        if temp == 1:
            print(f'Relatedness between {ancestry_1[0][0]} and {ancestry_2[0][0]}')
            print(f'rsd = 2 · {desc_inbreeding} / √(1 + {sire_inbreeding})(1 + {dam_inbreeding}) = {relatedness}')
            print('------------------- Next level -------------------')
        return relatedness


print(recursive_relatedness([[21], [13, 6], [10, 6, 0, 0]],
                            [[18], [12, 2], [6, 1, 0, 0]], 0))

