#!/bin/env/python

import argparse
import sys
import os
from numpy import random as npr
try:
    from person import Person, mating
except:
    from .person import Person, mating, recombine, add_denovo_hms
import random
import numpy as np

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

def generate_mate(founders, numSNPs, randoctr):
    randomname = "RANDOM-indiv"
    mate = None
    if len(founders) != 0:
        mate = founders.pop(0)
    if mate is None:
        randommatename = randomname+str(randoctr)
        mate = generate_random_person(randommatename, numSNPs)
        #print("Random mate: \n%s" % mate)
        randoctr += 1
    return mate, randoctr

def generate_pedigree(founder1Name, founder2Name, founder_genotype_phase_matrix, reproductionRate, generationNumber, marriageRate, all_founders, chrom_num, num_loci, sel_coeff):
    generationLists = []
    founders_shuffled = generate_shuffle(all_founders)
    num_founders = len(all_founders)
    numSNPs = founder_genotype_phase_matrix.shape[1]

    #Either find the founders or randomly select them
    ## Ensure that we remove them from the shuffled populations
    if founder1Name is not None:
        founder1 = all_founders["i"+str(founder1Name)]
        founders_shuffled.remove("i"+str(founder1Name))
    else:
        founder1Name = founders_shuffled.pop(0)
        founder1 = all_founders[founder1Name]

    if founder2Name is not None:
        founder2 = all_founders["i"+str(founder2Name)]
        founders_shuffled.remove("i"+str(founder2Name))
    else:
        founder2Name = founders_shuffled.pop(0)
        founder2 = all_founders[founder2Name]
    print("Founder1 will be %s" % (founder1))
    print("Founder2 will be %s" % (founder2))
    generationLists.append(founder1)
    generationLists.append(founder2)

    #determine number in first round of children, make it nonzero
    noChildren = abs(generate_random_normal_INT(reproductionRate))+1
    print("They will have %s children" % noChildren)

    #print(num_founders)
    #print(num_loci)
    #print(chrom_num)

    #recombine founder genotypes, no need to add de-novo since it's the FOUNDER and they already have them
    recombine_founder_genotype_phase_matrix = recombine(founder_genotype_phase_matrix, chrom_num, num_loci, num_founders)


    generation = "g1-i"
    currGen = []
    gen_genotypes = []
    for i in range(noChildren):
        childname = generation + str(i)
        #no current gen genotype matrix hence the None
        child, child_gen = mating(childname, founder1, founder2, recombine_founder_genotype_phase_matrix, None)
        child.set_genotype(i)
        print("Child #%s created" % (i))
        currGen.append(child)
        generationLists.append(child)
        gen_genotypes.append(child_gen)
    current_gen_phase_matrix = np.vstack(gen_genotypes)

    randompersonctr = 1 #count of individuals completely randomly generated
    gen_cont = False #is this generation continuing? i.e. is the number of kids for the next gen nonzero?

    #loop through the next generations we need to generate
    for generations in range(2,generationNumber):

        gen_genotypes = []
        #generate partner and new genotype matrix of the CURRENT generation. Founder matrix stays the same
        curr_postrecomb_matrix = recombine(current_gen_phase_matrix, chrom_num, num_loci, num_founders)
        curr_postdenovo_matrix,sel_coeff,chrom_num= add_denovo_hms(curr_postrecomb_matrix, sel_coeff, chrom_num, num_loci)
        #loop through every person in the current generation. Their index corresponds to the row of them in the 
        for eachperson in range(len(currGen)):
            newGen = []

            #will they marry? / do we need them to marry to continue gen?
            if generate_marriage_choice(marriageRate) or (gen_cont==False and eachperson==len(currGen)-1):
                partner, randompersonctr = generate_mate(founders_shuffled, numSNPs, randompersonctr)
                partner = all_founders[partner]

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
                    print(partner)
                    print(currGen[eachperson])
                    childname = "G" + str(generations) + "-"+partner.get_name()+","+currGen[eachperson].get_name()+"-C"+str(i+1)
                    child, child_gen = mating(childname, currGen[eachperson], partner, recombine_founder_genotype_phase_matrix, curr_postdenovo_matrix)
                    child.set_genotype(i)
                    print("Child #%s created" % (i+1))
                    newGen.append(child)
                    gen_cont = True
                    generationLists.append(child)
                    gen_genotypes.append(child_gen)
            curr_gen_phase_matrix = np.vstack(gen_genotypes)
            currGen = newGen
    return generationLists

def main():
    return None

if __name__ == "__main__":
    main()
