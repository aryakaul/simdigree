#!/bin/env/python

import sys
import os
from .person import Person
import numpy as np

def readin_matrix(filein):
    people = {}
    ctr = 1
    with open(filein, 'r') as f:
        for lines in f:
            line = lines.rstrip().split()
            isaffect = line[0]
            affect_decoder = {"U":False, "A":True}
            genotype = line[1:]
            for i in range(len(genotype)): genotype[i] = int(genotype[i])
            id = "indiv"+str(ctr)
            person = Person(id, genotype, affected=affect_decoder[isaffect])
            #print(person)
            people[person.get_name()] = person
            ctr+=1
            numsnps = len(genotype)
    return people, numsnps

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
            #print(person)
            allpeople[person.get_name()] = person
            snp_matrix.append(genotype)
            ctr+=1
    return allpeople, snp_matrix

def readin_effectmatrix(filein):
    effectcolvec = []
    with open(filein, 'r') as f:
        for lines in f:
            line = lines.rstrip().split()
            effectcolvec.append(int(line))
    return effectcolvec

def gen_healthy_affected(effectcolvec, snp_matrix):
    healthy = {}
    affected = {}
    #TODO 

def main():
    readin_matrix("../tests/test-matrix_4ppl-5snps")

if __name__ == "__main__":
    main()
