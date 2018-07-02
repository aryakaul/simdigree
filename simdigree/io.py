#!/bin/env/python

import sys
import os
try:
    from person import Person
except:
    from .person import Person
import numpy as np

def readin_vcf(filein):
    for lines in filein:
        if lines.startswith("#"): continue
        line = lines.rstrip().split()
        info = line[8]
        gt = line[10:]

def readin_snpmatrix(filein):
    allpeople = {}
    snp_matrix = []
    ctr = 1
    with open(filein, 'r') as f:
        for lines in f:
            line = lines.rstrip().split()
            genotype = line
            for i in range(len(genotype)): genotype[i] = int(genotype[i])
            id = "indiv"+str(ctr)
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

def write_pedigree(listofpeople, outfile):
    with open(outfile, 'w') as f:
        for people in listofpeople:
            line = []
            line.append("1")
            line.append(people.get_name())
            try:
                line.append(people.get_parents()[0])
                line.append(people.get_parents()[1])
            except:
                line.append(".")
                line.append(".")
            line.append("3")
            for snp in people.get_genotype():
                line.append(str(snp))
            linestr = '\t'.join([str(i) for i in line])
            linestr += "\n"
            #linestr = str(*line, sep="\t")
            f.write(linestr)

def main():
    allpeople, snp_matrix = readin_snpmatrix("../tests/matrices/test-snpmatrix_4ppl-5snps")
    effectcolvec = readin_effectvec("../tests/matrices/test-effectvec_4ppl-5snps")
    healthy,affected = gen_healthy_affected(np.array(effectcolvec).T, np.asmatrix(snp_matrix), 0.01, allpeople)

if __name__ == "__main__":
    main()
