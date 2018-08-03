#!/bin/env/python

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

    """
    Returns a number from a normal distribution
    with given mean
    """

    return int(npr.normal(loc=mean))

def generate_random_person(name, numSNPs):

    """
    Generate a random individual with a bunch of
    random snps.
    """

    choices = [0,1,2,3]
    rand_gen = list(npr.choice(choices, numSNPs))
    randopeep = Person(name, rand_gen)
    randopeep.set_genotype_snps(rand_gen)
    return randopeep

def generate_marriage_choice(prob):

    """
    We roll a dice and if it's greater than the
    given probability of a marriage, we return True
    """

    flip = npr.random_sample()
    return (prob>flip)

def generate_shuffle(dic):

    """
    Given a dictionary, shuffle the keys
    """

    keys = list(dic.keys())
    return (random.sample(keys, len(keys)))

def generate_mate(founders, numSNPs, randoctr):

    """
    Given a shuffled list of founders. Determine if we need to generate,
    a random mate. If we do, generate one and return it. If not, just return
    the founder
    """

    randomname = "random-indiv"
    mate = None
    if len(founders) != 0:
        mate = founders.pop(0)
    if mate is None:
        randommatename = randomname+str(randoctr)
        mate = generate_random_person(randommatename, numSNPs)
        #print("Random mate: \n%s" % mate)
        randoctr += 1
    return mate, randoctr

def generate_pedigree(founder_genotype_phase_matrix, reproductionRate, generationNumber, marriageRate, all_founders, chrom_num, num_loci, sel_coeff):

    """
    Given a lot of information, generate a novel, stochastically determined pedigree.
    """

    #generated pedigree genotype matrix
    GT_MATRIX = []

    #this list will store the necessary information for a later FAM and PED file output
    generationLists = []

    #shuffle the founders for random ordering
    founders_shuffled = generate_shuffle(all_founders)

    gt_matrix_ctr = 0

    # Randomly determine initial founders
    # Ensure that we remove them from the shuffled populations
    founder1Name = founders_shuffled.pop(0)
    founder1 = all_founders[founder1Name]
    founder1_genotype = founder1.get_genotype(founder_genotype_phase_matrix, None)
    GT_MATRIX.append(founder1_genotype)
    founder1.gt_matrix_ctr = gt_matrix_ctr
    gt_matrix_ctr += 1
    founder2Name = founders_shuffled.pop(0)
    founder2 = all_founders[founder2Name]
    founder2_genotype = founder2.get_genotype(founder_genotype_phase_matrix, None)
    GT_MATRIX.append(founder2_genotype)
    founder2.gt_matrix_ctr = gt_matrix_ctr
    gt_matrix_ctr += 1
    print("Founder1 will be %s" % (founder1))
    print("Founder2 will be %s" % (founder2))
    generationLists.append(founder1)
    generationLists.append(founder2)
    GT_MATRIX = np.vstack(GT_MATRIX)
    RECOMBINED_GT_MATRIX = recombine(GT_MATRIX, chrom_num, num_loci)

    #determine number in first round of children, make it nonzero
    noChildren = abs(generate_random_normal_INT(reproductionRate))+1
    print("They will have %s children" % noChildren)


    #create first generation of kiddies
    generation = "g1-i"
    currGen = []
    for i in range(noChildren):
        childname = generation + str(i)
        child, child_gen = mating(childname, founder1, founder2, RECOMBINED_GT_MATRIX)
        print("%s created" % (child.get_name()))
        generationLists.append(child)
        currGen.append(child)
        RECOMBINED_GT_MATRIX=np.vstack((RECOMBINED_GT_MATRIX, child_gen))
        child.gt_matrix_ctr = gt_matrix_ctr
        gt_matrix_ctr+=1

    GT_MATRIX = RECOMBINED_GT_MATRIX

    #loop through the next generations we need to generate
    for generations in range(1,generationNumber):
        newGen = []
        gen_cont = False
        for eachperson in range(len(currGen)):

            #will they marry? / do we need them to marry to continue gen?
            if generate_marriage_choice(marriageRate) or (gen_cont==False and eachperson==len(currGen)-1):
                try:
                    partnerName = founders_shuffled.pop(0)
                except:
                    print("Do not have enough founders to mate. Please run again and set a higher value for -f")
                    sys.exit(2)
                partner = all_founders[partnerName]
                partner_genotype = partner.get_genotype(founder_genotype_phase_matrix, None)
                GT_MATRIX = np.vstack((GT_MATRIX,partner_genotype))
                partner.gt_matrix_ctr = gt_matrix_ctr
                gt_matrix_ctr += 1
                generationLists.append(partner)

                # generate the partnership 
                print("Marriage between %s and %s" % (currGen[eachperson].get_name(),partner.get_name()))

                # how many children will they have?
                noChildren = abs(generate_random_normal_INT(reproductionRate))

                # if we need a child to sustain generation...
                if eachperson==len(currGen)-1 and noChildren == 0 and gen_cont == False:
                    noChildren += 1
                if noChildren == 0: continue
                print("They will have %s children" % noChildren)
                RECOMBINED_GT_MATRIX = recombine(GT_MATRIX, chrom_num, num_loci)
                CURR_GT_MATRIX = []
                # create each child
                for i in range(noChildren):
                    childname = "g" + str(generations+1) + "-c"+str(len(newGen))
                    child, child_gen = mating(childname, currGen[eachperson], partner, RECOMBINED_GT_MATRIX)
                    print("%s created" % (childname))
                    newGen.append(child)
                    generationLists.append(child)
                    CURR_GT_MATRIX.append(child_gen)
                    child.gt_matrix_ctr = gt_matrix_ctr
                    gt_matrix_ctr+=1
                    gen_cont = True
                CURR_GT_MATRIX = np.vstack(CURR_GT_MATRIX)
                DENOVO_CURR_GT_MATRIX,sel_coeff,chrom_num=add_denovo_hms(CURR_GT_MATRIX, sel_coeff, chrom_num, num_loci)

                numdenovo = DENOVO_CURR_GT_MATRIX.shape[1]-RECOMBINED_GT_MATRIX.shape[1]
                if numdenovo != 0:
                    zero_col = np.zeros((RECOMBINED_GT_MATRIX.shape[0], numdenovo))
                    RECOMBINED_GT_MATRIX = np.hstack((RECOMBINED_GT_MATRIX, zero_col))
                RECOMBINED_GT_MATRIX = np.vstack((RECOMBINED_GT_MATRIX, DENOVO_CURR_GT_MATRIX))
                GT_MATRIX = RECOMBINED_GT_MATRIX

        # overwrite the new current generation
        currGen = newGen

    return generationLists, sel_coeff, GT_MATRIX

def main():
    return None

if __name__ == "__main__":
    main()
