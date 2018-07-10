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

    allpeople = {}
    for i in range(len(geno_dosage)):
        allpeople["i"+str(founders_to_subset[i])]=Person("founder-i"+str(founders_to_subset[i]), i)
        x = list(geno_dosage[i,:])
        allpeople["i"+str(founders_to_subset[i])].set_genotype_snps(x)

    return geno_dosage, abs(sel_coeff), allele_freqs, chrom_num, allpeople

def parse_inputNames(founder):
    if founder.startswith("indiv"):
        try:
            return int(founder.split("v")[1])
        except:
            print("Illegal argument found for founder: %s\nSee help message" % founder)
    elif founder.startswith("i"):
        try:
            return int(founder.split("i")[1])
        except:
            print("Illegal argument found for founder: %s\nSee help message" % founder)
    else:
        try:
            return int(founder)
        except:
            print("Illegal argument found for founder: %s\nSee help message" % founder)

def readin_snpmatrix(filein):
    allpeople = {}
    snp_matrix = []
    ctr = 1
    with open(filein, 'r') as f:
        for lines in f:
            line = lines.rstrip().split()
            genotype = line
            for i in range(len(genotype)): genotype[i] = int(genotype[i])
            id = "i"+str(ctr)
            person = Person(id, genotype)
            allpeople[person.get_name()] = person
            snp_matrix.append(genotype)
            ctr+=1
    return allpeople, snp_matrix

def readin_effectvec(filein):
    effectcolvec = []
    with open(filein, 'r') as f:
        for lines in f:
            line = lines.rstrip()
            effectcolvec.append(float(line))
    return effectcolvec

def gen_healthy_affected(effectcolvec, snp_matrix, liabilityThreshold, allpeople):
    cutoff=int(round(len(allpeople)*liabilityThreshold))
    if cutoff==0: cutoff += 1
    cutoff = len(allpeople)-cutoff
    people_affected_value = np.matmul(snp_matrix, effectcolvec)
    argsort = np.argsort(people_affected_value)
    healthy = {}
    affected = {}
    for indiv in allpeople:
        number = int(indiv.split("v")[1])-1
        if number >= cutoff:
            person = allpeople[indiv]
            person.set_affected(True)
            affected[indiv] = person
        else:
            person = allpeople[indiv]
            person.set_affected(False)
            healthy[indiv] = person
    return healthy,affected

def determine_healthy_affected(snpfilein, effectfilein, liabilityThres):
    allpeople, snp_matrix = readin_snpmatrix(snpfilein)
    effectcolvec = readin_effectvec(effectfilein)
    healthy,affected = gen_healthy_affected(np.array(effectcolvec).T, np.asmatrix(snp_matrix), liabilityThres, allpeople)
    return healthy, affected

def write_fam(listofpeople, outfile):
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
            if people.is_affected() is True:
                line.append("2")
            elif people.is_affected() is False:
                line.append("1")
            else:
                line.append("0")
            linestr = '\t'.join([str(i) for i in line])
            linestr += "\n"
            #linestr = str(*line, sep="\t")
            f.write(linestr)

def read_fam(famfile):
    individs = {}
    with open(famfile, 'r') as f:
        for lines in f:
            line = lines.rstrip().split()
            famid = line[0]
            individ = line[1]
            father = line[2]
            mother = line[3]
            x = [father, mother]
            x.sort()
            if father == "0" and mother == "0": founder = True
            else: founder = False
            if individ not in individs:
                dummy = Dummy((x[0], x[1]), founder)
                individs[individ] = dummy
        for dummies in individs:
            dummy = individs[dummies]
            if dummy.is_founder(): continue
            parents = dummy.get_parents()
            f = individs[parents[0]]
            m = individs[parents[1]]
            f.add_child(dummies)
            m.add_child(dummies)
    return individs

def main():
    pass

if __name__ == "__main__":
    main()
