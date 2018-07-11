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
    founder = 0
    with open(pedigree, 'r') as filein:
        for lines in filein:
            line = lines.rstrip().split()
            if line[2] == "0" and line[3] == "0": founder += 1

    return founder

def recreate_pedigree(individs, dummypairs, all_founders, founder_genotype_phase_matrix, num_loci, chrom_num, sel_coeff):
    generationLists = []
    founders_shuffled = generate_shuffle(all_founders)
    num_founders = len(all_founders)
    numSNPs = founder_genotype_phase_matrix.shape[1]

    recombine_founder_genotype_phase_matrix = recombine(founder_genotype_phase_matrix, chrom_num, num_loci, num_founders)

    founder_pair = None
    for pairname in dummypairs:
        pair = dummypairs[pairname]
        par1 = pair.get_pair()[0]
        par2 = pair.get_pair()[1]
        if individs[par1].is_founder() and individs[par2].is_founder():
            print("Founder pair found")
            print(pair)
            founder_pair = pair
            dummypairs.pop(pairname)
            break
    founderName = founders_shuffled.pop(0)
    founder1Name = all_founders[founderName]
    founder1 = all_founders[founderName]
    founderName = founders_shuffled.pop(0)
    founder2Name = all_founders[founderName]
    founder2 = all_founders[founderName]
    print("Marriage between %s and %s" % (founder1.get_name(), founder2.get_name()))
    founder1.set_genotype_snps(list(recombine_founder_genotype_phase_matrix[founder1.genotype,:]))
    founder2.set_genotype_snps(list(recombine_founder_genotype_phase_matrix[founder2.genotype,:]))
    generationLists.append(founder1)
    generationLists.append(founder2)
    generation = "g1-i"
    currGen = []
    gen_genotypes = []
    print("They will have %s children" % founder_pair.get_num_children())
    for i in range(founder_pair.get_num_children()):
        childname = founder_pair.get_children()[i]
        #no current gen genotype matrix hence the None below
        child, child_gen = mating(childname, founder1, founder2, recombine_founder_genotype_phase_matrix, None)
        child.set_genotype(i)
        child.set_genotype_snps(child_gen)
        print("Child #%s created" % (i+1))
        currGen.append(child)
        generationLists.append(child)
        gen_genotypes.append(child_gen)
    current_gen_phase_matrix = np.vstack(gen_genotypes)
    randompersonctr = 1
    while len(dummypairs) != 0:
        generations = 2
        curr_postrecomb_matrix = recombine(current_gen_phase_matrix, chrom_num, num_loci, num_founders)
        curr_postdenovo_matrix,sel_coeff,chrom_num= add_denovo_hms(curr_postrecomb_matrix, sel_coeff, chrom_num, num_loci)
        newGen = [] #store kids created here
        for children in currGen:
            for pairs in dummypairs:
                if children.get_name() in pairs:
                    randompersonctrprev = randompersonctr
                    partner, randompersonctr = generate_mate(founders_shuffled, numSNPs, randompersonctr)
                    if randompersonctr == randompersonctrprev:
                        partner = all_founders[partner]
                        partner.set_genotype_snps(list(recombine_founder_genotype_phase_matrix[partner.genotype,:]))
                    generationLists.append(partner)
                    print("Marriage between %s and %s" % (children.get_name(), partner.get_name()))
                    noChildren = dummypairs[pairs].get_num_children()
                    print("They will have %s children" % noChildren)
                    for i in range(noChildren):
                        childname = dummypairs[pairs].get_children()[i]
                        child, child_gen = mating(childname, children, partner, recombine_founder_genotype_phase_matrix, curr_postdenovo_matrix)
                        child.set_genotype(i)
                        child.set_genotype_snps(child_gen)
                        print("Child #%s created" % (i+1))
                        newGen.append(child)
                        gen_cont = True
                        generationLists.append(child)
                        gen_genotypes.append(child_gen)
                    dummypairs.pop(pairs)
                    break
        currGen = newGen
        current_gen_phase_matrix = np.vstack(gen_genotypes)
    return generationLists, sel_coeff
