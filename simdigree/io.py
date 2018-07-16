#!/bin/env/python 
import numpy.random as npr
import sys
import os
try:
    from person import Person
except:
    from .person import Person
import numpy as np
try:
    from dummy import Dummy
except:
    from .dummy import Dummy
import pandas as pd

#modified from Onuralp's code
def read_vcf(path, founders_desired):

    """
    Extract selective effects and allele frequencies of causal variants from VCF generated using SLiM. Note that the
    header and non-causal variants (those outside the mutational target) were removed from the original VCF output.
    """

    # Read in the vcf file and extract a couple useful things:
    ## founder_genotype_phase_matrix
    #### A matrix where each row is an individual, and each column is a SNP
    #### Every value is either a 0, 1, 2, or 3. A 0 corresponds to a homozygous
    #### reference genotype, a 1 is a '0|1' genotype, a 2 is a '1|0' genotype
    #### and a 3 is a homozygous alternate genotype.
    ## selection_coeff
    #### A column vector with a row for each SNP. Every value denotes the
    #### selective effect of each SNP, as defined by SLiM
    ## allele_freqs
    #### Column vector with row for each SNP, and the value associated with
    #### each SNP's allele frequency.
    ## chrom_num
    #### Column vector with row for each SNP, and the value associated with
    #### the chromosomal location of each SNP
    ## all_founders
    #### a dictionary with a key, value pair where the key is 'i' followed
    #### by the number where it appears in the founder matrix and the value
    #### associated with the Person object created with that founder
    #read in vcf file

    df = pd.read_csv(path, delim_whitespace=True, header=None, comment='#', keep_default_na=False, low_memory=False, engine='c', dtype='str')
    no_samples = df.shape[1]-9

    # extract selective effects `s` and allele counts `ac`, calculate allele frequencies `x`
    sel_coeff = df[7].str.extract(r';S=(\S+);DOM', expand=False).values.astype(float)
    allele_counts = df[7].str.extract(r';AC=(\S+);DP', expand=False).values.astype(int)
    allele_freqs = allele_counts / (2 * no_samples)

    # track locus id
    chrom_num = df[0].values

    #generate random founders to subset
    founders_to_subset=npr.randint(low=9,high=no_samples+9, size=founders_desired)
    print("Founders will be the following... %s" % founders_to_subset)

    # construct genotype matrix (one-hot encoding)
    genotypes = df.loc[:, founders_to_subset].T.values
    dosage = {'0|0': 0, '0|1': 1, '1|0': 2, '1|1': 3}
    geno_dosage = np.vectorize(dosage.get)(genotypes)

    #create dictionary with everybody 
    allpeople = {}
    for i in range(len(geno_dosage)):
        allpeople["i"+str(founders_to_subset[i])]=Person("founder-i"+str(founders_to_subset[i]), i)
        x = list(geno_dosage[i,:])
        allpeople["i"+str(founders_to_subset[i])].set_genotype_snps(x)

    return geno_dosage, abs(sel_coeff), allele_freqs, chrom_num, allpeople

def read_vcf_founderliab(path):

    """
    Read whole vcf and return ONLY founder matrix
    """

    df = pd.read_csv(path, delim_whitespace=True, header=None, comment='#', keep_default_na=False, low_memory=False, engine='c', dtype='str')
    no_samples = df.shape[1]-9

    # construct genotype matrix (one-hot encoding)
    genotypes = df.loc[:, founders_to_subset].T.values
    dosage = {'0|0': 0, '0|1': 1, '1|0': 2, '1|1': 3}
    geno_dosage = np.vectorize(dosage.get)(genotypes)

    return geno_dosage

def write_fam(listofpeople, outfile):

    """
    Write the generation lists to the fam file
    """

    with open(outfile, 'w') as f:
        for people in listofpeople:
            line = []
            line.append("1") #FID
            line.append(people.get_name()) #IID
            try:
                line.append(people.get_parents()[0]) # parents
                line.append(people.get_parents()[1])
            except:
                line.append("0")
                line.append("0")
            line.append("0") #sex
            if people.is_affected() is True: #affected status
                line.append("2")
            elif people.is_affected() is False:
                line.append("1")
            else:
                line.append("0") # should never occur
            linestr = '\t'.join([str(i) for i in line])
            linestr += "\n"
            f.write(linestr)

def read_fam(famfile):

    """
    Read the fam file and create Dummy's for each person
    """

    individs = {}
    with open(famfile, 'r') as f:
        for lines in f:
            # extract line information
            line = lines.rstrip().split()
            famid = line[0]
            individ = line[1]
            father = line[2]
            mother = line[3]
            x = [father, mother]
            # sort so deterministic
            x.sort()

            # if no parents known, then it's a founder
            if father == "0" and mother == "0": founder = True
            else: founder = False

            # create dummy
            if individ not in individs:
                dummy = Dummy((x[0], x[1]), founder)
                individs[individ] = dummy

        # loop through dummies
        for dummies in individs:
            dummy = individs[dummies]
            if dummy.is_founder(): continue

            # add child to the mother and father Dummy's
            parents = dummy.get_parents()
            f = individs[parents[0]]
            m = individs[parents[1]]
            f.add_child(dummies)
            m.add_child(dummies)
    return individs

def write_effect_snps(colvec, output):

    """
    Write the column vector to the file
    """

    with open(output, 'w') as fi:
        for x in colvec:
            fi.write("%s\n" % x)

def write_genotype_matrix(generations, output):

    """
    Write the dosage matrix to a file
    """

    gen_matrix = {}
    todosage = {0:0, 1:1, 2:2, 3:2} # dictionary to convert to dosage matrix
    for i in generations: gen_matrix[i.get_name()] = i.get_genotype_snps()
    max_len = 0

    # find maximum length of snps
    for g in gen_matrix:
        if len(gen_matrix[g]) > max_len: max_len = len(gen_matrix[g])

    # add zeros to those snps who don't have max len
    for g in gen_matrix:
        diff = max_len - len(gen_matrix[g])
        if diff < 0:
            print("error diff = %s" % diff)
            sys.exit(2)
        elif diff == 0: continue
        else:
            zeros = np.zeros(diff)
            gen_matrix[g]=np.append(gen_matrix[g], zeros)
            gen_matrix[g]=np.vectorize(todosage.get)(gen_matrix[g])

    # write the string to the file described
    with open(output, 'w') as fileout:
        for g in gen_matrix:
            stringToWrite = []
            stringToWrite.append(g)
            for j in gen_matrix[g]: stringToWrite.append(str(j))
            linestr = '\t'.join(stringToWrite)
            fileout.write(linestr+"\n")

def main():
    pass

if __name__ == "__main__":
    main()
