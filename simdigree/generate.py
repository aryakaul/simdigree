#!/bin/env/python

import argparse
import sys
import os
from numpy import random as npr
from .person import Person

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


def generate_pedigree(people, founder1, founder2, reproductionRate, generationNumber, marriageRate):
        ped_out = []
        
        people,numSNPs = readin_matrix(args.input)
        founder1 = people["indiv"+str(args.ancestor1)]
        founder2 = people["indiv"+str(args.ancestor2)]
        print("Founder1 will be %s" % (founder1))
        print("Founder2 will be %s" % (founder2))
        noChildren = abs(generate_random_normal_INT(args.reproduction))+1
        print("They will have %s children" % noChildren)
        generation = "G1-FOUNDER-C"
        currGen = []
        for i in range(noChildren):
            childname = generation + str(i+1)
            child = mating(childname, founder1, founder2)
            print("Child #%s created" % (i+1))
            currGen.append(child)
        randomname = "RANDOM-indiv"
        randompersonctr = 1
        gen_cont = False
        for generations in range(2,args.generation):
            for eachperson in range(len(currGen)):
                newGen = []
                if generate_marriage_choice(args.marriagerate) or (gen_cont==False and eachperson==len(currGen)-1):
                    partner = generate_random_person(randomname+str(randompersonctr), numSNPs)
                    print("Marriage between %s and %s" % (currGen[eachperson].get_name(),partner.get_name()))
                    randompersonctr+=1
                    noChildren = abs(generate_random_normal_INT(args.reproduction))
                    if eachperson==len(currGen)-1 and noChildren == 0 and gen_cont == False:
                        noChildren += 1
                    for i in range(noChildren):
                        childname = "G" + str(generations) + "-"+partner.get_name()+"-C"+str(i+1)
                        child = mating(childname, currGen[eachperson], partner)
                        print("Child #%s created" % (i+1))
                        newGen.append(child)
                        gen_cont = True



def main():
    args = parser_args(sys.argv[1:])

if __name__ == "__main__":
    main()
