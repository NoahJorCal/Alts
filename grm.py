#!/usr/bin/python3

from scipy import sparse
# from simulator import Individual


class Pedigree:
    def __init__(self, size):
        self.__individuals = list(range(1, size + 1))
        self.__sires = [0] * size
        self.__dams = [0] * size
        # Relatedness or kinship matrix
        self.__kinship = sparse.lil_matrix((size, size))
        # Diagonal values are set to 1 (inbreeding)
        self.__kinship.setdiag(1)
        self.__inds_removed = 0

    @property
    def kinship(self):
        return self.__kinship

    @property
    def pedigree(self):
        message = f'Inds:  {" ".join(map(str, self.__individuals))}\n' \
                  f'Sires: {" ".join(map(str, self.__sires))}\n' \
                  f'Dams:  {" ".join(map(str, self.__dams))}'
        return message

    def calculate_kinship(self, ind):
        """ When accessing the matrix, the row will be the higher
        index as only the bottom left half of the matrix is filled """
        ind_id = ind.id - self.__inds_removed
        # Individual's ids start on 1 but indexes start on 0
        ind_index = ind_id - 1
        sire_index = ind.sire - 1 - self.__inds_removed
        dam_index = ind.dam - 1 - self.__inds_removed
        self.__kinship.resize(ind_id, ind_id)
        # If the parents are unknown the inbreeding is by default 1
        if ind.sire == 0 and ind.dam == 0:
            self.__kinship[ind_index, ind_index] = 1
        else:
            # Base inbreeding + inbreeding from the parents
            self.__kinship[ind_index, ind_index] = 1 + (0.5 * self.__kinship[max(sire_index, dam_index),
                                                                             min(sire_index, dam_index)])
            for col in range(ind_index):
                # Kinship between the last individual added and the individual in the column col
                kinship = 0.5 * (self.__kinship[max(sire_index, col), min(sire_index, col)] +
                                 self.__kinship[max(dam_index, col), min(dam_index, col)])
                self.__kinship[ind_index, col] = kinship

    def add_individual(self, individual):
        self.__individuals.append(individual.id)
        self.__sires.append(individual.sire)
        self.__dams.append(individual.dam)
        self.calculate_kinship(individual)

    def trim_kinship_matrix(self, ind_id):
        ind_id = ind_id - self.__inds_removed
        self.__kinship = self.__kinship[ind_id:, ind_id:]
        self.__inds_removed += ind_id

    def relatedness(self, ind_list):
        max_ind_id = max(ind_list) - 1 - self.__inds_removed
        min_ind_id = min(ind_list) - 1 - self.__inds_removed
        return self.__kinship[max_ind_id, min_ind_id]


# ped = Pedigree(3)
#
# sim = 1
# a = Individual(sim)
# a.id = 1
# ped.add_individual(a)
#
# a = Individual(sim)
# a.id = 2
# ped.add_individual(a)
#
# a = Individual(sim)
# a.id = 6
# ped.add_individual(a)
#
# a = Individual(sim)
# a.id = 10
# ped.add_individual(a)
#
# a = Individual(sim)
# a.id = 12
# a.sire = 1
# a.dam = 6
# ped.add_individual(a)
#
# a = Individual(sim)
# a.id = 13
# a.sire = 6
# a.dam = 10
# ped.add_individual(a)
#
# a = Individual(sim)
# a.id = 18
# a.sire = 12
# a.dam = 2
# ped.add_individual(a)
#
# a = Individual(sim)
# a.id = 21
# a.sire = 6
# a.dam = 13
# ped.add_individual(a)
# print(ped.relatedness([21, 18]))
#
