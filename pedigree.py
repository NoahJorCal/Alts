#!/usr/bin/python3

from scipy import sparse
import numpy
# from simulator import Individual


class Pedigree:
    def __init__(self, size):
        self.__individuals = list(range(1, size + 1))
        self.__sires = [0] * size
        self.__dams = [0] * size
        # Relatedness or kinship matrix
        self.__kinship = sparse.lil_matrix((size, size))
        # Diagonal values are set to 1 (consanguinity)
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
        ''' When accessing the matrix, the row will be the higher
        index as only the bottom left half of the matrix is filled '''
        ind_id = ind.id - self.__inds_removed
        # Individual's ids start on 1 but indexes start on 0
        ind_index = ind_id - 1
        sire_index = ind.sire - 1 - self.__inds_removed
        dam_index = ind.dam - 1 - self.__inds_removed
        self.__kinship.resize(ind_id, ind_id)
        # If the parents are unknown the consanguinity is by default 1
        print(ind.id, ind_index)
        if ind.sire == 0 and ind.dam == 0:
            print('Parents unknown')
            self.__kinship[ind_index, ind_index] = 1
            print(numpy.around(self.__kinship.todense(), decimals=1))
        else:
            # Base consanguinity + consanguinity from the parents
            print(f'Consanguinity: {1 + (0.5 * self.__kinship[max(sire_index, dam_index), min(sire_index, dam_index)])}')
            self.__kinship[ind_index, ind_index] = 1 + (0.5 * self.__kinship[max(sire_index, dam_index),
                                                                             min(sire_index, dam_index)])
            print(f'sire_index: {sire_index}, dam_index: {dam_index}')
            for col in range(ind_index):
                # Kinship between the last individual added and the individual in the column col
                kinship = 0.5 * (self.__kinship[max(sire_index, col), min(sire_index, col)] +
                                 self.__kinship[max(dam_index, col), min(dam_index, col)])
                self.__kinship[ind_index, col] = kinship
                # print(f'0.5 * ({self.__kinship[max(sire_index, col), min(sire_index, col)]} + {self.__kinship[max(dam_index, col), min(dam_index, col)]}) = {kinship}')
        print(numpy.around(self.__kinship.todense(), decimals=1))

    def add_individual(self, individual):
        self.__individuals.append(individual.id)
        self.__sires.append(individual.sire)
        self.__dams.append(individual.dam)
        self.calculate_kinship(individual)
        # print(self.__kinship)

    def trim_kinship_matrix(self, ind_id):
        ind_id = ind_id - self.__inds_removed
        print('Pre-trim')
        print(numpy.around(self.__kinship.todense(), decimals=1))
        self.__kinship = self.__kinship[ind_id:, ind_id:]
        self.__inds_removed += ind_id
        print('Post-trim')
        print(numpy.around(self.__kinship.todense(), decimals=1))

    def relatedness(self, ind_list):
        print(ind_list)
        max_ind_id = max(ind_list) - 1 - self.__inds_removed
        min_ind_id = min(ind_list) - 1 - self.__inds_removed
        print(self.__kinship.shape)
        return self.__kinship[max_ind_id, min_ind_id]

# ped = Pedigree(2)
# km = ped.kinship
# # print(km.todense())
# a = 1
#
# indiv = Individual(a)
# indiv.id = 1
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 2
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 3
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 4
# indiv.sire = 2
# indiv.dam = 1
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 5
# indiv.sire = 1
# indiv.dam = 2
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 6
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 7
# indiv.sire = 3
# indiv.dam = 4
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 8
# indiv.sire = 4
# indiv.dam = 3
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 9
# indiv.sire = 6
# indiv.dam = 5
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 10
# ped.add_individual(indiv)
#
# km = ped.kinship
# print('FIN')
# print(km.todense())
#
# print()
# print('========================================================================================================')
# print()
# ped = Pedigree(2)
# km = ped.kinship
# # print(km.todense())
# a = 1
#
# indiv = Individual(a)
# indiv.id = 1
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 2
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 3
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 4
# indiv.sire = 1
# indiv.dam = 2
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 5
# indiv.sire = 3
# indiv.dam = 4
# ped.add_individual(indiv)
#
# indiv = Individual(a)
# indiv.id = 6
# indiv.sire = 4
# indiv.dam = 5
# ped.add_individual(indiv)
#
# print('Pre-trim')
# print(ped.kinship.todense())
# ped.trim_kinship_matrix(3)
# print('Post-trim')
# print(ped.kinship.todense())
#
# indiv = Individual(a)
# indiv.id = 7
# indiv.sire = 4
# indiv.dam = 6
# ped.add_individual(indiv)
#
# km = ped.kinship
# print('FIN')
# print(km.todense())
#
# print(ped.relatedness([5, 7]))