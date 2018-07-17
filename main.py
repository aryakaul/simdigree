#!/net/home/akaul/miniconda3/bin/python3.6

import time
import argparse
import sys
import math
import os
from simdigree.io import read_vcf, write_fam, read_fam, write_effect_snps, write_genotype_matrix, read_vcf_founderliab
from simdigree.generate import *
from simdigree.pedigree import recreate_pedigree, return_subset_number
from simdigree.person import mating
from simdigree.dummy import DummyPair


def parser_args(args):
    parser = argparse.ArgumentParser(prog="simdigree")
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    parser_pedigree = subparsers.add_parser('pedigree', help="Use a sample pedigree file?")
    parser_pedigree.add_argument('-i', '--inputvcf', help="Path to the vcf file being used", type=str, required=True)
    parser_pedigree.add_argument('-T', '--noLoci', help="Number of loci being used", type=int, required=True)
    parser_pedigree.add_argument('-t', '--tau', help="Tau value for effect calculation. Default is 0.5", type=float, required=False, nargs='*')
    parser_pedigree.add_argument('-l', '--liabilityThreshold', help="Float value(s) representing the liability threshold percentage(s) being used. Default is 0.01", type=float, required=False, nargs = '*')
    parser_pedigree.add_argument('-p', '--fam', help="Path to the ped file to use", type=str, required=True)
    parser_pedigree.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is working directory under /simdigree_output", type=str, required=False, default="./simdigree_output")



    parser_generate = subparsers.add_parser('generate', help="Make a fake family history")
    parser_generate.add_argument('-i', '--inputvcf', help="Path to the vcf file being used", type=str, required=True)
    parser_generate.add_argument('-T', '--noLoci', help="Number of loci being used", type=int, required=True)
    parser_generate.add_argument('-f', '--founders', help="Number of founders matrix should be subsetted to. Default is 15", type=int, required=False, default=15)
    parser_generate.add_argument('-t', '--tau', help="Tau value to use for effect calculation. Default is 0.5", type=float, required=False, nargs='*')
    parser_generate.add_argument('-l', '--liabilityThreshold', help="Float value(s) representing the liability threshold percentage(s) being used. Default is 0.01", type=float, required=False, nargs = '*')
    parser_generate.add_argument('-g', '--generation', help="How many generations will be created total?", type=int, required=False, default=4)
    parser_generate.add_argument('-r', '--reproduction', help="Given a parental pair, what's the mean number of kids they'll have? (If desired # of gen. is not achieved, I force at least one child)", type=int, required=False, default=3)
    parser_generate.add_argument('-m', '--marriagerate', help="Given an individual, what's the the probability they will marry? (If desired # of gen. is not achieved, I force at least one marriage", type=float, required=False, default=0.5)
    parser_generate.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is in the working directory under '/simdigree_output'", type=str, required=False, default="./simdigree_output")
    return parser.parse_args(args)

def calculate_dosage_matrix(founder_genotype_phase_matrix):

    """
    Given a matrix describing phased genotypes (0: 0|0, 1: 0|1, 2: 1|0, 3: 1|1)
    Return a matrix describing good ole genotypes (0: Homozygous reference,
    1: Heterozygous, 2: Homozygous alternate)
    """

    doses = {0:0, 1:1, 2:2, 3:2}
    x = np.copy(founder_genotype_phase_matrix)
    for k, v in doses.items(): x[founder_genotype_phase_matrix==k] = v
    return x

def calculate_scaling_constant(s, maf, tau, h2=0.5):

    """
    Calculate scaling constant given a selection coefficient, a minor allele
    frequency, a tau value and ?h2?
    """

    C = np.sqrt(h2 / np.sum(2 * maf * (1 - maf) * ((s ** tau) ** 2)))
    if math.isnan(C):
        print("WARNING - Negative value found during square root for scaling constant...\nC = 1")
        C = 1
    return C

def calculate_liability(X, s, C, tau, h2=0.5):

    """
    Calculate 'liability' given a good ole genotype matrix, the selection
    coefficient, the scaling constant, tau and ?h2? Return the effect
    """

    b = (s ** tau) * C
    G = X @ b
    G = G.reshape(len(G))
    E = np.sqrt(1. - h2) * np.random.randn(len(G))
    return G + E, b

def main():
    args = parser_args(sys.argv[1:])
    startTHUGE = time.time()

    #find the number of founders we want to subset the vcf too --> this is key to efficiency
    if args.command == "pedigree":
        subset_number = return_subset_number(args.fam)
    elif args.command == "generate":
        subset_number = args.founders
    else:
        print("Unrecognized sub-command, valid options are {pedigree, generate}. Command given: %s" % args.command)
        sys.exit(2)
    print("VCF will be subsetted to the following number of founders: %s" % subset_number)

    #read in vcf - check io.py for more info
    start = time.time()
    founder_genotype_phase_matrix, selection_coeff, allele_freqs, chrom_num, all_founders = read_vcf(args.inputvcf, subset_number)
    end = time.time()
    print("Time took to read in vcf file... %s" % (end-start))

    #create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    basepath = args.output


    if args.command == "generate":

        #generate the novel pedigree, check generate.py for more details
        start = time.time()
        generations, selection_coeff_new = generate_pedigree(founder_genotype_phase_matrix, args.reproduction, args.generation, args.marriagerate, all_founders, chrom_num, args.noLoci, selection_coeff)
        end = time.time()
        print("Time took to generate pedigree... %s" % (end-start))

    elif args.command=="pedigree":

        # create a set of marriage pairs for all the people in the family
        individs = read_fam(args.fam)
        pairs = {}
        for i in individs:
            dummy = individs[i]
            if dummy.is_founder(): continue
            pair = DummyPair(dummy.get_parents())
            if pair.get_pair() not in pairs:
                pairs[pair.get_pair()] = pair
            else: pair = pairs[pair.get_pair()]
            p1Children = individs[dummy.get_parents()[0]].get_children()
            p2Children = individs[dummy.get_parents()[1]].get_children()
            # following line in case one person has multiple children with multiple people
            pair.children = list(set(p1Children) & set(p2Children))

        #check pedigree.py for more info.
        start = time.time()
        generations, selection_coeff_new = recreate_pedigree(individs, pairs, all_founders, founder_genotype_phase_matrix, args.noLoci, chrom_num, selection_coeff)
        end = time.time()
        print("Time took to model pedigree... %s" % (end-start))

    #set default tauValues and liabilityThreshold values to loop through
    if args.tau is None:
        tauValues = [0.5]
    else:
        tauValues = args.tau

    if args.liabilityThreshold is None:
        lT = [0.01]
    else:
        lT = args.liabilityThreshold



    #start looping over all tauValues
    for tauValue in tauValues:

        #create output path
        taupath = os.path.join(basepath, "tau-"+str(tauValue)+"/")
        if not os.path.exists(taupath):
            os.makedirs(taupath)
        startbig = time.time()
        print("Calculating all values for tau value of %s" % tauValue)

        #go through founders and determine the threshold for the given values 
        start = time.time()
        founder_biggenotype_mat = read_vcf_founderliab(args.inputvcf)
        founder_dosage = calculate_dosage_matrix(founder_biggenotype_mat)
        C = calculate_scaling_constant(selection_coeff, allele_freqs, tauValue)
        effects_people, effects_snps = calculate_liability(founder_dosage, selection_coeff, C, tauValue)
        end = time.time()
        print("Time took to calculate liability of founders... %s" % (end-start))


        # figure out the founders who are affected at each liability threshold
        listOfThresholds = []
        founder_outliers = []
        for i in lT:
            listOfThresholds.append(np.percentile(effects_people, 100*(1-i)))
            founder_outliers.append(np.argwhere(effects_people > i))
        print("Looping through all liability thresholds")
        for tidx in range(len(lT)):
            print("Number of outliers in the founder population... %s" % founder_outliers[tidx].shape[0])

            # create output path for diff. liab. thresholds
            outputpath = os.path.join(taupath, "lb-"+str(lT[tidx])+"/")
            if not os.path.exists(outputpath):
                os.makedirs(outputpath)
            print("On liability threshold of %s" % lT[tidx])
            start = time.time()
            founder_derived_threshold = listOfThresholds[tidx]
            print("Founder derived threshold is... %s" % founder_derived_threshold)

            # create genotype matrix from all individuals in final pedigree
            gt_matrix = []
            maxlen = 0
            nonfounder_ctr = 0
            for person in generations:

                # if they're a founder, then their phenotype is already known!
                if person.is_founder():
                    personName = person.get_name()
                    founderNo = int(personName.split('i')[1])
                    if founderNo in founder_outliers[tidx]:
                        person.set_affected(True)
                    else:
                        person.set_affected(False)

                #if they're not a founder...
                else:
                    # get the maxlen of the snps
                    snps = person.get_genotype_snps()
                    if len(snps) > maxlen: maxlen = len(snps)

                    # attach their snps to the matrix
                    gt_matrix.append(snps)
                    person.ctr_liab = nonfounder_ctr
                    nonfounder_ctr += 1

            # create new gt_matrix with all genotypes filled in for non-founders to account for denovo mutations
            gt_matrix_wdenovo = []
            for rows in gt_matrix:
                diff = maxlen - len(rows)
                if diff == 0: row = rows
                elif diff > 0:
                    zeros = np.zeros(diff)
                    row = np.append(rows,zeros)
                else:
                    print("ERROR")
                    sys.exit(2)
                gt_matrix_wdenovo.append(row)
            gt_matrix = np.vstack(gt_matrix_wdenovo)

            # convert this to a dosage matrix and calculate allele freq.
            nonfounder_dosage = calculate_dosage_matrix(gt_matrix)
            allele_freq_nonfounder = []
            for column in nonfounder_dosage.T:
                genotypes, freqs = (np.unique(column, return_counts = True))
                freqdict = {}
                for i in range(len(genotypes)): freqdict[genotypes[i]] = freqs[i]
                if 0 not in freqdict: freqdict[0] = 0
                if 1 not in freqdict: freqdict[1] = 0
                if 2 not in freqdict: freqdict[2] = 0
                allele_freq_persnp = freqdict[1]+2*freqdict[2]
                allele_freq_persnp /= (2*(freqdict[1]+freqdict[2]+freqdict[0]))
                allele_freq_nonfounder.append(allele_freq_persnp)
            allele_freq_nonfounder = np.vstack(allele_freq_nonfounder)
            # calculate the non-founders who are effected per founder_derived_threshold
            C = calculate_scaling_constant(selection_coeff, allele_freq_nonfounder, tauValue)
            effects_nonfounders, effects_snps_wdenovo = calculate_liability(nonfounder_dosage, selection_coeff_new, C, tauValue)
            nonfounder_outliers = np.argwhere(effects_nonfounders > founder_derived_threshold)

            #set affected status
            for people in generations:
                if people.is_founder(): continue
                pos = people.ctr_liab
                if pos not in nonfounder_outliers or len(nonfounder_outliers)==0:
                    people.set_affected(False)
                else:
                    people.set_affected(True)
            end = time.time()
            print("Time took to calculate who is affected... %s" % (end-start))

            # write three files:
            # (1) the fam file with information about the affected status of all the individuals
            # (2) the column vector describing the effects of each SNP
            # (3) the genotype matrix describing dosage of each individual
            write_fam(generations, os.path.join(outputpath, "simdigree_out-liabThreshold-"+str(lT[tidx])+"-tau-"+str(tauValue)+".fam"))
            write_effect_snps(effects_snps_wdenovo, os.path.join(outputpath, "simdigree_out-effects_snps"))
            write_genotype_matrix(generations, os.path.join(outputpath, "simdigree_out-dosage"))
        endbig = time.time()
        print("Time to calculate given tau value was %s" % (endbig-startbig))

    endTHUGE = time.time()
    print("TOTAL TIME TOOK... %s" % (endTHUGE-startTHUGE))
if __name__ == "__main__":
    main()
