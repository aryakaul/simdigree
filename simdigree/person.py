#!/bin/env/python

import random
import sys

class Person:
    def __init__(self, name, genotype, parents=None, childrennames=None, affected="unknown"):
        self.name = name
        self.genotype = genotype
        self.parents = parents
        self.children = childrennames
        self.affected = affected

    def __str__(self):
        return("Person: %s\nGenotype: %s\nParents: %s\nChildren: %s\nAffected?: %s\n" % (self.get_name(), self.get_genotype(), self.get_parents(), self.get_children(), self.is_affected()))


    def get_name(self):
        return self.name

    def get_genotype(self):
        return self.genotype

    def get_parents(self):
        if self.parents:
            return self.parents
        else:
            #print("No parents from individual")
            return None

    def get_children(self):
        return self.children

    def is_affected(self):
        if self.affected==True:
            return True
        elif self.affected==False:
            return False
        else:
            return "unknown"

    def set_parents(self, parent1id, parent2id):
        self.parents = (parent1id, parent2id)

    def set_affected(self, affected):
        self.affected = affected

    def add_children(self, childrenid):
        if self.get_children() == None:
            childlist = []
        else:
            childlist = self.children
        childlist.append(childrenid)
        self.children = childlist

def mating(child_id, parent1, parent2):
    child_gen = genotype_mating(parent1, parent2)
    child = Person(child_id, child_gen)
    parent1.add_children(child_id)
    parent2.add_children(child_id)
    child.set_parents(parent1.get_name(), parent2.get_name())
    return child

def genotype_mating(person1, person2):
    person1Genotype = person1.get_genotype()
    person2Genotype = person2.get_genotype()

    if len(person1Genotype) != len(person2Genotype):
        print("ERROR - Genotypes are not the same length")
        print("%s - Person 1's # of genotypes" % len(person1Genotype))
        print("%s - Person 2's # of genotypes" % len(person2Genotype))
        sys.exit(2)

    genotype = []
    for snps in range(len(person1Genotype)):
        currSNP1 = person1Genotype[snps]
        currSNP2 = person2Genotype[snps]
        babysnp = snp_mating(currSNP1, currSNP2)
        genotype.append(babysnp)
    return genotype

def snp_mating(snp1genotype, snp2genotype):
    possiblecrossings = {
        "hetcrosshet":(1,1),
        "hetcrosshomalt":(1,2),
        "hetcrosshomref":(0,1),
        "homaltcrossref":(0,2),
        "homrefcrosshomref":(0,0),
        "homaltcrosshomalt":(2,2),
    }

    ourPairing = tuple(sorted([snp1genotype, snp2genotype]))
    realcross = None
    for gens in possiblecrossings:
        if ourPairing == possiblecrossings[gens]:
            realcross = gens
            #print(ourPairing)
            #print(realcross)

    if realcross == "hetcrosshet":
        prob = [0.25,0.75,1]
        sides = [0,1,2]

    elif realcross == "hetcrosshomalt":
        prob = [0.5,1]
        sides = [1,2]

    elif realcross == "hetcrosshomref":
        prob = [0.5,1]
        sides = [0,1]

    elif realcross == "homaltcrossref":
        prob = [1]
        sides = [1]

    elif realcross == "homrefcrosshomref":
        prob = [1]
        sides = [0]

    elif realcross == "homaltcrosshomalt":
        prob = [1]
        sides = [2]

    else:
        print("Impossible cross being studied")
        print("%s=first parent's genotype" % (snp1genotype))
        print("%s=second parent's genotype" % (snp2genotype))
        sys.exit(2)

    babyGen = rollDie(prob, sides)
    return babyGen

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
