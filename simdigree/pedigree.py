import sys
import numpy as np
import os
from numpy import random as npr
try:
    from person import Person, mating, recombine, add_denovo_hms
except:
    from .person import Person, mating, recombine, add_denovo_hms
try:
    from generate import generate_shuffle, generate_mate
except:
    from .generate import generate_shuffle, generate_mate

def return_subset_number(pedigree):

    """
    Find out how many founders there are in the given file.
    """

    founder = 0
    with open(pedigree, 'r') as filein:
        for lines in filein:
            line = lines.rstrip().split()
            if line[2] == "0" and line[3] == "0": founder += 1

    return founder

def recreate_pedigree(individs, dummypairs, all_founders, founder_genotype_phase_matrix, num_loci, chrom_num, sel_coeff):

    """
    Given a variety of inputs recreate the pedigree
    """

    #pedigree genotype matrix
    GT_MATRIX = []

    #list that will store all the people created
    generationLists = []

    #shuffle the founders so random
    founders_shuffled = generate_shuffle(all_founders)
    num_founders = len(all_founders)
    numSNPs = founder_genotype_phase_matrix.shape[1]

    #set all founders and add them to new genotype matrix
    gt_matrix_ctr = 0
    founderdummy_to_founder = {}
    for people in individs:
        if individs[people].is_founder():
            if people not in founderdummy_to_founder:
                founderName = founders_shuffled.pop(0)
                founder = all_founders[founderName]
                founderdummy_to_founder[people] = founder
                founder_genotype = founder.get_genotype(founder_genotype_phase_matrix, None)
                GT_MATRIX.append(founder_genotype)
                founder.gt_matrix_ctr = gt_matrix_ctr
                gt_matrix_ctr += 1
    GT_MATRIX = np.vstack(GT_MATRIX)


    #loop through dummy pairs and find all founder pairs
    founder_pairs = []
    for pairname in dummypairs:
        pair = dummypairs[pairname]
        par1 = pair.get_pair()[0]
        par2 = pair.get_pair()[1]
        if individs[par1].is_founder() and individs[par2].is_founder():
            print("Founder pair found")
            print(pair)
            founder_pairs.append(pair)
    RECOMBINED_GT_MATRIX = recombine(GT_MATRIX, chrom_num, num_loci)

    # loop through founder pairs
    currGen = []
    gen_genotypes = []
    for founder_pair in founder_pairs:
        dummypairs.pop(founder_pair.get_pair())
        f1 = founder_pair.get_pair()[0]
        founder1 = founderdummy_to_founder[f1]
        f2 = founder_pair.get_pair()[1]
        founder2 = founderdummy_to_founder[f2]
        print("Marriage between %s and %s" % (founder1.get_name(), founder2.get_name()))
        generationLists.append(founder1)
        generationLists.append(founder2)

        # determine number of children and loop through them
        print("They will have %s children" % founder_pair.get_num_children())
        for i in range(founder_pair.get_num_children()):
            childname = founder_pair.get_children()[i]
            #no current gen genotype matrix hence the None below
            child, child_gen = mating(childname, founder1, founder2, RECOMBINED_GT_MATRIX)
            #child.set_genotype_snps(child_gen)
            print("%s created" % (child.get_name()))
            generationLists.append(child)
            currGen.append(child)
            RECOMBINED_GT_MATRIX=np.vstack((RECOMBINED_GT_MATRIX, child_gen))
            child.gt_matrix_ctr = gt_matrix_ctr
            gt_matrix_ctr+=1

    GT_MATRIX = RECOMBINED_GT_MATRIX

    # loop through all remaining dummypairs
    while len(dummypairs) > 0:

        newGen = [] #store kids created here

        # loop through all children
        for children in currGen:
            toRemove = []

            # loop through remaining dummy pairs
            for pairs in dummypairs:

                # if the current child is part of a marriage pair
                if children.get_name() in pairs:

                    # find the child's partner
                    # NOTE doesn't allow consanguineous rels. also implies that we have enough founders
                    lp = list(pairs)
                    lp.remove(children.get_name())
                    partner = lp[0]
                    partner = founderdummy_to_founder[partner]
                    print("Marriage between %s and %s" % (children.get_name(), partner.get_name()))
                    generationLists.append(partner)
                    noChildren = dummypairs[pairs].get_num_children()
                    print("They will have %s children" % noChildren)

                    CURR_GT_MATRIX = []
                    # Do recombination for this pair
                    RECOMBINED_GT_MATRIX = recombine(GT_MATRIX, chrom_num, num_loci)
                    # loop through that pair's no of children
                    for i in range(noChildren):
                        childname = dummypairs[pairs].get_children()[i]
                        child, child_gen = mating(childname, children, partner, RECOMBINED_GT_MATRIX)
                        print("%s created" % (child.get_name()))
                        newGen.append(child)
                        generationLists.append(child)
                        CURR_GT_MATRIX.append(child_gen)
                        child.gt_matrix_ctr = gt_matrix_ctr
                        gt_matrix_ctr+=1
                    CURR_GT_MATRIX = np.vstack(CURR_GT_MATRIX)
                    DENOVO_CURR_GT_MATRIX,sel_coeff,chrom_num=add_denovo_hms(CURR_GT_MATRIX, sel_coeff, chrom_num, num_loci)
                    numdenovo = DENOVO_CURR_GT_MATRIX.shape[1]-RECOMBINED_GT_MATRIX.shape[1]
                    if numdenovo != 0:
                        zero_col = np.zeros((RECOMBINED_GT_MATRIX.shape[0], numdenovo))
                        RECOMBINED_GT_MATRIX = np.hstack((RECOMBINED_GT_MATRIX, zero_col))
                    RECOMBINED_GT_MATRIX = np.vstack((RECOMBINED_GT_MATRIX, DENOVO_CURR_GT_MATRIX))

                    toRemove.append(pairs)
                    # reset currGen
                    currGen = newGen
                    GT_MATRIX = RECOMBINED_GT_MATRIX

            # remove the pairs we went through
            for k in toRemove: dummypairs.pop(k)

    # return generation lists and new selection coefficients
    return generationLists, sel_coeff, GT_MATRIX
