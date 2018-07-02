#!/bin/env/python

import argparse
import sys
import os
from numpy import random as npr
try:
    from person import Person, mating
except:
    from .person import Person, mating
import random

def generate_random_normal_INT(mean):
    return int(npr.normal(loc=mean))

def generate_random_person(name, numSNPs):
    choices = [0,1,2]
    rand_gen = list(npr.choice(choices, numSNPs))
    randopeep = Person(name, rand_gen)
    return randopeep

def generate_marriage_choice(prob):
    flip = npr.random_sample()
    return (prob>flip)

def generate_shuffle(dic):
    keys = list(dic.keys())
    return (random.sample(keys, len(keys)))

def generate_healthy_or_affected():
    die = npr.rand()
    if die > 0.5: return "healthy"
    else: return "affected"

def generate_mate(healthy, affected, healthy_shuffle, affected_shuffle, numSNPs, randoctr):
    randomname = "RANDOM-indiv"
    choice = generate_healthy_or_affected()
    mate = None
    print(choice)
    if choice == "healthy" or len(affected_shuffle)==0:
        try:
            mate = healthy_shuffle.pop(0)
            print("Healthy mate: \n%s" % mate)
            mate = healthy[mate]
        except:
            pass
    elif choice == "affected" or len(healthy_shuffle)==0:
        try:
            mate = affected_shuffle.pop(0)
            print("Affected mate: \n%s" % mate)
            mate = affected[mate]
        except:
            pass
    if mate is None:
        randommatename = randomname+str(randoctr)
        mate = generate_random_person(randommatename, numSNPs)
        print("Random mate: \n%s" % mate)
        randoctr += 1
    return mate, randoctr


def generate_pedigree(healthy, affected, healthy_shuffle, affected_shuffle, founder1, founder2, reproductionRate, generationNumber, marriageRate):
    generationLists = []
    print(healthy_shuffle)
    #Either find the founders or randomly select them
    ## Ensure that we remove them from the shuffled populations
    if founder1 is not None:
        try:
            founder1 = healthy["indiv"+str(founder1)]
            healthy_shuffle.remove(founder1)
        except:
            founder1 = affected["indiv"+str(founder1)]
            affected_shuffle.remove(founder1)
    else:
        founder1 = healthy_shuffle.pop(0)
        founder1 = healthy[founder1]
    if founder2 is not None:
        try:
            founder2 = healthy["indiv"+str(founder2)]
            healthy_shuffle.remove(founder2)
        except:
            founder2 = affected["indiv"+str(founder2)]
            affected_shuffle.remove(founder2)
    else:
        founder2 = affected_shuffle.pop(0)
        founder2 = affected[founder2]
    print("Founder1 will be %s" % (founder1))
    print("Founder2 will be %s" % (founder2))
    generationLists.append(founder1)
    generationLists.append(founder2)
    #determine first round of children, make it nonzero
    noChildren = abs(generate_random_normal_INT(reproductionRate))+1
    print("They will have %s children" % noChildren)
    generation = "G1-FOUNDER-C"
    currGen = []
    for i in range(noChildren):
        childname = generation + str(i+1)
        child = mating(childname, founder1, founder2)
        print("Child #%s created" % (i+1))
        currGen.append(child)
        generationLists.append(child)

    numSNPs = len(founder1.get_genotype())
    randompersonctr = 1 #count of individuals completely randomly generated
    gen_cont = False #is this generation continuing? i.e. is the number of kids for the next gen nonzero?

    #loop through the next generations we need to generate
    for generations in range(2,generationNumber):

        #loop through every person in the current generation
        for eachperson in range(len(currGen)):
            newGen = []

            #will they marry? / do we need them to marry to continue gen?
            if generate_marriage_choice(marriageRate) or (gen_cont==False and eachperson==len(currGen)-1):
                partner, randompersonctr = generate_mate(healthy, affected, healthy_shuffle, affected_shuffle, numSNPs, randompersonctr)
                generationLists.append(partner)

                # generate the partnership 
                print("Marriage between %s and %s" % (currGen[eachperson].get_name(),partner.get_name()))

                # how many children will they have?
                noChildren = abs(generate_random_normal_INT(reproductionRate))
                print("They will have %s children" % noChildren)

                # if we need a child to sustain generation...
                if eachperson==len(currGen)-1 and noChildren == 0 and gen_cont == False:
                    noChildren += 1

                # create each child
                for i in range(noChildren):
                    childname = "G" + str(generations) + "-"+partner.get_name()+","+currGen[eachperson].get_name()+"-C"+str(i+1)
                    child = mating(childname, currGen[eachperson], partner)
                    print("Child #%s created" % (i+1))
                    newGen.append(child)
                    gen_cont = True
                    generationLists.append(child)
    return generationLists

def main():
    return None

if __name__ == "__main__":
    main()
