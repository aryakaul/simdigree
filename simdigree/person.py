#!/bin/env/python

import numpy as np
import random
import sys

class Person:
    def __init__(self, name, genotype=None, parents=None, childrennames=None, affected="unknown", genotype_snps=None):
        self.name = name
        self.genotype = genotype
        self.parents = parents
        self.children = childrennames
        self.affected = affected
        self.genotype_snps = genotype_snps

    def __str__(self):
        return("Person: %s\nParents: %s\nChildren: %s\nAffected?: %s\n" % (self.get_name(), self.get_parents(), self.get_children(), self.is_affected()))

    def is_founder(self):
        return self.get_name().startswith("founder")

    def get_name(self):
        return self.name

    def get_genotype(self, founder_genotype_matrix, curr_gen_genotype_matrix):
        genidx = self.genotype
        if self.is_founder(): return list(founder_genotype_matrix[genidx,:])
        elif self.is_random(): return self.get_genotype_snps()
        else: return list(curr_gen_genotype_matrix[genidx,:])

    def get_parents(self):
        if self.parents:
            return self.parents
        else:
            #print("No parents from individual")
            return None

    def get_children(self):
        return self.children

    def is_random(self):
        return self.get_name().startswith("random")

    def is_affected(self):
        if self.affected==True:
            return True
        elif self.affected==False:
            return False
        else:
            return "unknown"

    def get_genotype_snps(self):
        return self.genotype_snps

    def set_parents(self, parent1id, parent2id):
        self.parents = (parent1id, parent2id)


    def set_affected(self, affected):
        self.affected = affected

    def set_genotype(self, genotype):
        self.genotype = genotype

    def add_children(self, childrenid):
        if self.get_children() == None:
            childlist = []
        else:
            childlist = self.children
        childlist.append(childrenid)
        self.children = childlist

    def set_genotype_snps(self, snps):
        self.genotype_snps = snps

def mating(child_id, parent1, parent2, founder_mat, currgen_genmat):
    child_gen = phase_parents(parent1, parent2, founder_mat, currgen_genmat)
    parent1.add_children(child_id)
    parent2.add_children(child_id)
    child = Person(child_id)
    child.set_parents(parent1.get_name(), parent2.get_name())
    return child, child_gen

#modified from Onuralp
def random_breaks(no_loci, coldspot=0.1):
    """
    Generate random breakpoints between `no_loci` unlinked loci from $Poisson(\lambda=0.20)$.
    For `coldspot` fraction of total mutational target, set $\lambda=0.02$
    (10-fold lower rate of recombination within coldspot region).
    """
    r = 1e-08 # genome-wide recombination rate per base pair per generation (uniform)
    L = 4.5e03 # mutational target size per locus

    target_rate = r*L*no_loci
    whole = np.random.poisson(lam=target_rate, size=no_loci)
    cold_region = int(no_loci * coldspot)
    start = np.random.randint(low=1, high=no_loci-cold_region)
    end = start + cold_region
    whole[start:end] = np.random.poisson(lam=0.02, size=cold_region)

    # 1-based list of loci at which a breakpoint occurs
    breakpoints = np.nonzero(whole)[0] + 1

    return breakpoints

#modified from Onuralp
def designate_breakpoints(chr_index, no_loci, no_samples):
    """
    Generate boolean mask from chromosome index
    """
    mask = np.empty((no_samples, len(chr_index)), dtype=bool)
    bp_list = []

    for i in range(no_samples):
        bp = random_breaks(no_loci, coldspot=0.1)
        mask[i] = np.isin(chr_index, [str(x) for x in bp])
        bp_list.append(bp)
    return mask, bp_list

#modified from Onuralp
def recombine(geno_matrix, chr_index, no_loci, no_samples):
    """
    Recombine at randomly generated breakpoints.
    """
    recomb = {0: 0, 1: 2, 2: 1, 3: 3} # '0|1' <-> '1|0'
    masked, bp_list = designate_breakpoints(chr_index, no_loci, no_samples)
    z = np.copy(geno_matrix)
    if np.asarray(bp_list).size > 0:
        # this would modify the original geno_matrix too! Work with copy!
        try:
            z[masked] = np.vectorize(recomb.get)(z[masked])
        except:
            return z
    return z#, bp_list

#modified from Onuralp
def add_denovo_hms(X, s, old_index, num_loci,gamma_shape=0.31623, gamma_scale=0.01):
    '''
    Given inherited genotypes `X`, return genotype matrix with de novo mutations.
    Append selective effects of de novo mutations to the base selection vector.
    '''

    U = 2.36e-08 # genome-wide mutation rate per base pair per generation
    L = 4.5e03 # mutational target size per locus (in base pairs)
    N = num_loci # number of causal loci

    a, b = X.shape
    denovo_per_sample = np.random.poisson(lam=U*L*N, size=a)
    total_num_denovo = np.sum(denovo_per_sample)

    #denovo_per_sample = np.array([0,0,1,0,0,0,0,0,1,0])
    #total_num_denovo = np.sum(denovo_per_sample)

    # draw selective effect for each de novo mutation from Gamma
    k, theta = gamma_shape, gamma_scale
    s_denovo = np.random.gamma(shape=k, scale=theta, size=total_num_denovo)
    s_concat = np.concatenate((s, s_denovo))

    # calculate effect sizes for de novo and append
    #b_denovo = (s_denovo ** tau) * scaling_factor
    #b_concat = np.concatenate((b, b_denovo))

    # pre-allocate genotypes for de novo mutations
    X_denovo = np.tile(np.zeros(total_num_denovo), (a,1))

    # sample genotypes for de novo 
    j = 0
    for i,x in enumerate(denovo_per_sample):
        sampled = np.random.choice([1,2], size=x)
        X_denovo[i][j:j+x] = sampled
        j = j + x

    X_concat = np.concatenate((X, X_denovo), axis=1)

    # where to add these de novo mutations
    # update chrom index

    chrom_denovo = np.random.choice(range(1,N+1), size=total_num_denovo)

    new_index = np.concatenate((old_index, chrom_denovo), axis=0)

    return X_concat, s_concat, new_index #, b_concat

#modified from Onuralp
def phase_parents(parent1, parent2, founder_mat, currgen_mat):
    parent1 = parent1.get_genotype(founder_mat, currgen_mat)
    parent2 = parent2.get_genotype(founder_mat, currgen_mat)
    """
    Randomly sample gametes from each parent and generate offspring
    genotypes.
    """
    gt = {0: [0, 0], 1: [0, 1], 2: [1, 0], 3: [1, 1]}
    a, b = np.random.randint(low=0, high=2, size=2)
    p1 = np.array([gt[x][a] for x in parent1])
    p2 = np.array([gt[x][b] for x in parent2])

    ## in case one has denovo mutations. the length of their genotypes will be different. 
    ## If this is the case, add zeros to one st. it goes 
    while len(p1) > len(p2):
        np.append(p2, 0)
    while len(p2) > len(p1):
        np.append(p1, 0)

    proband = (2 * p1) + p2

    return proband

def rollDie(prob, sides):
    x = random.random()
    if 0<x<=prob[0]:
        return sides[0]
    for i in range(1, len(prob)):
        if prob[i-1]<x<=prob[i]:
            return sides[i]

def main():
   print(snp_mating(2,1))

if __name__ == "__main__":
    main()
