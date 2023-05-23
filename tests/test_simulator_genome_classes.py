import random
import numpy as np
import pytest
import sys
from os import path
# Add the alts directory to the Python path
sys.path.append(path.join(path.dirname(__file__), '..'))
from simulator import Chromosome
from simulator import Locus
from simulator import Genome


@pytest.mark.parametrize('locus_size, mutation_rate, expected_mutations', [
    (1000, 0, 0),
    (1000, 0.001, 1),
    (10000, 0.001, 9),
])
def test_chromosome_mutate(locus_size, mutation_rate, expected_mutations):
    np.random.seed(909848)
    random.seed(909848)
    chromosome = Chromosome('neutral')
    chromosome.mutate(locus_size, mutation_rate)
    is_ordered = all(chromosome.snvs[i] <= chromosome.snvs[i + 1] for i in range(len(chromosome.snvs) - 1))
    assert len(chromosome.snvs) == expected_mutations
    assert is_ordered


@pytest.mark.parametrize('snvs_list, locus_size, expected_sequence', [
    ([0.0267, 0.7647, 0.9273], 10,
     [1, 0, 0, 0, 0, 0, 0, 0, 1, 1]),
    ([0.3162, 0.3653, 0.5736, 0.8176], 20,
     [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0]),
    ([0.5258, 0.5412], 10,
     [0, 0, 0, 0, 0, 1, 1, 0, 0, 0]),
    ([1], 10,
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
])
def test_chromosome_snvs_to_sequence(snvs_list, locus_size, expected_sequence):
    chromosome = Chromosome('neutral')
    chromosome.snvs = snvs_list
    result_sequence = list(chromosome.snvs_to_sequence(locus_size))
    assert result_sequence == expected_sequence


@pytest.mark.parametrize('crossovers, snvs_lists, expected_snvs', [
    ([0.317834, 0.864235],
     [[0.145689, 0.368904, 0.398149, 0.678918, 0.918726],
     [0.125613, 0.236673, 0.286352, 0.362346, 0.743542]],
     [[0.145689, 0.362346, 0.743542, 0.918726],
     [0.125613, 0.236673, 0.286352, 0.368904, 0.398149, 0.678918]]),
    ([0.472436, 0.823462, 0.458345, 0.612478, 0.492346],
     [[0.145689, 0.152345, 0.172341, 0.216233, 0.262347],
     [0.283414, 0.294523, 0.367234, 0.398172, 0.437593]],
     [[0.145689, 0.152345, 0.172341, 0.216233, 0.262347],
      [0.283414, 0.294523, 0.367234, 0.398172, 0.437593]]),
    ([0.7435421, 0.996172, 0.146234, 0.362346, 0.916752],
     [[0.145689, 0.184245, 0.523686, 0.912552, 0.918726],
     [0.274513, 0.583425, 0.592345, 0.834525]],
     [[0.145689, 0.274513, 0.523686, 0.834525, 0.918726],
     [0.184245, 0.583425, 0.592345, 0.912552]]),
])
def test_locus_recombine(crossovers, snvs_lists, expected_snvs):
    chromosomes = []
    for snvs_list in snvs_lists:
        chrom = Chromosome('neutral')
        chrom.snvs = snvs_list
        chromosomes.append(chrom)
    locus = Locus(chromosomes, 'neutral', 1000, 0, 0.001)
    locus.recombine(crossovers=crossovers)
    result_snvs = [list(chromosome.snvs) for chromosome in locus.chromosomes]
    assert result_snvs == expected_snvs


@pytest.mark.parametrize('snvs, expected_haplotype', [
    ([[[0.0267, 0.7647, 0.9273], [0.5736, 0.8176]],
     [[0.1893, 0.4895], [0.3653, 0.7274, 0.9863]]],
     ([[1, 0, 0, 0, 1], [0, 0, 0, 1, 1]],
      [[0, 1, 1, 0, 0], [0, 0, 1, 0, 1]]))
])
def test_genome_haplotype(snvs, expected_haplotype):
    loci = []
    for locus in snvs:
        chromosomes = []
        for snvs_list in locus:
            chrom = Chromosome('neutral')
            chrom.snvs = snvs_list
            chromosomes.append(chrom)
        loci.append(Locus(chromosomes, 'neutral', 5, 0, 0))
    genome = Genome(loci)
    haplotype = genome.haplotype()
    result_haplotype = []
    for locus in haplotype:
        locus_haplotype = []
        for chromosome in locus:
            locus_haplotype.append(list(chromosome))
        result_haplotype.append(locus_haplotype)
    result_haplotype = tuple(result_haplotype)
    assert result_haplotype == expected_haplotype










