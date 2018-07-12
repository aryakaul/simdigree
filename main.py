#!/net/home/akaul/miniconda3/bin/python3.6

import time
import argparse
import sys
import math
import os
from simdigree.io import parse_inputNames, read_vcf, write_fam, read_fam, write_effect_snps, write_genotype_matrix
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
    parser_generate.add_argument('-p1', '--ancestor1', help="Line number of individual who will be one ancestral parent. Default is a randomly selected founder.", type=str, required=False, default=None)
    parser_generate.add_argument('-p2', '--ancestor2', help="Line number of individual who will be the other ancestral parent. Default is a randomly selected founder.", type=str, required=False, default=None)
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



    # Read in the vcf file and extract a couple useful things:
    ## founder_genotype_phase_matrix
    #### A matrix where each row is an individual, and each column is a SNP
    #### Every value is either a 0, 1, 2, or 3. A 0 corresponds to a homozygous
    #### reference genotype, a 1 is a '0|1' genotype, a 2 is a '1|0' genotype
    #### and a 3 is a homozygous alternate genotype.
    ## selection_coeff
    #### A column vector with a row for each SNP. Every value denotes the
    #### selective effect of each SNP, as defined by SLiM
    ## allele_freqs
    #### Column vector with row for each SNP, and the value associated with
    #### each SNP's allele frequency.
    ## chrom_num
    #### Column vector with row for each SNP, and the value associated with
    #### the chromosomal location of each SNP
    ## all_founders
    #### a dictionary with a key, value pair where the key is 'i' followed
    #### by the number where it appears in the founder matrix and the value
    #### associated with the Person object created with that founder
    if args.command == "pedigree":
        subset_number = return_subset_number(args.fam)
        print("VCF will be subsetted to this number...%s" % subset_number)
    elif args.command == "generate":
        subset_number = args.founders
    else:
        print("Unrecognized sub-command, valid options are {pedigree, generate}. Command given: %s" % args.command)
        sys.exit(2)

    start = time.time()
    founder_genotype_phase_matrix, selection_coeff, allele_freqs, chrom_num, all_founders = read_vcf(args.inputvcf, subset_number)
    end = time.time()
    print("Time took to read in vcf file... %s" % (end-start))

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    basepath = args.output


    if args.command == "generate":
        #parse_inputNames allows users to put in 'indiv1' or 'i1'
        if args.ancestor1: founder1 = parse_inputNames(args.ancestor1)
        else: founder1 = None
        if args.ancestor2: founder2 = parse_inputNames(args.ancestor2)
        else: founder2 = None

        #TODO COMMENTS HERE
        start = time.time()
        generations, selection_coeff_new = generate_pedigree(founder1, founder2, founder_genotype_phase_matrix, args.reproduction, args.generation, args.marriagerate, all_founders, chrom_num, args.noLoci, selection_coeff)
        end = time.time()
        print("Time took to generate pedigree... %s" % (end-start))

    elif args.command=="pedigree":
        individs = read_fam(args.fam)
        pairs = {}
        for i in individs:
            dummy = individs[i]
            if dummy.is_founder(): continue
            pair = DummyPair(dummy.get_parents())
            if pair.get_pair() not in pairs:
                pairs[pair.get_pair()] = pair
            else: pair = pairs[pair.get_pair()]
            pair.children = individs[dummy.get_parents()[0]].get_children()
        start = time.time()
        generations, selection_coeff_new = recreate_pedigree(individs, pairs, all_founders, founder_genotype_phase_matrix, args.noLoci, chrom_num, selection_coeff)
        end = time.time()
        print("Time took to model pedigree... %s" % (end-start))


    if args.tau is None:
        tauValues = [0.5]
    else:
        tauValues = args.tau

    if args.liabilityThreshold is None:
        lT = [0.01]
    else:
        lT = args.liabilityThreshold

    for tauValue in tauValues:
        taupath = os.path.join(basepath, "tau-"+str(tauValue)+"/")
        if not os.path.exists(taupath):
            os.makedirs(taupath)
        startbig = time.time()
        print("Calculating all values for tau value of %s" % tauValue)
        start = time.time()
        founder_dosage = calculate_dosage_matrix(founder_genotype_phase_matrix)
        C = calculate_scaling_constant(selection_coeff, allele_freqs, tauValue)
        effects_people, effects_snps = calculate_liability(founder_dosage, selection_coeff, C, tauValue)
        end = time.time()
        print("Time took to calculate liability of founders... %s" % (end-start))


        listOfThresholds = []
        founder_outliers = []
        for i in lT:
            listOfThresholds.append(np.percentile(effects_people, 100*(1-i)))
            founder_outliers.append(np.argwhere(effects_people > i))

        print("Looping through all liability thresholds")
        for tidx in range(len(lT)):
            outputpath = os.path.join(taupath, "lb-"+str(lT[tidx])+"/")
            if not os.path.exists(outputpath):
                os.makedirs(outputpath)
            print("On liability threshold of %s" % lT[tidx])
            start = time.time()
            founder_derived_threshold = listOfThresholds[tidx]
            gt_matrix = []
            maxlen = 0
            nonfounder_ctr = 0
            for person in generations:
                if person.is_founder():
                    if person.genotype in founder_outliers[tidx]:
                        person.set_affected(True)
                    else:
                        person.set_affected(False)
                else:
                    snps = person.get_genotype_snps()
                    if len(snps) > maxlen: maxlen = len(snps)
                    gt_matrix.append(snps)
                    person.ctr_liab = nonfounder_ctr
                    nonfounder_ctr += 1
            for rows in gt_matrix:
                while len(rows) < maxlen:
                    np.append(rows, 0)
            gt_matrix = np.vstack(gt_matrix)
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
            C = calculate_scaling_constant(selection_coeff, allele_freq_nonfounder, tauValue)
            effects_nonfounders, effects_snps_wdenovo = calculate_liability(nonfounder_dosage, selection_coeff_new, C, tauValue)
            nonfounder_outliers = np.argwhere(effects_nonfounders > founder_derived_threshold)
            for people in generations:
                if people.is_founder(): continue
                pos = people.ctr_liab
                if pos not in nonfounder_outliers or len(nonfounder_outliers)==0:
                    people.set_affected(False)
                else:
                    people.set_affected(True)
            end = time.time()
            print("Time took to calculate who is affected... %s" % (end-start))
            write_fam(generations, os.path.join(outputpath, "simdigree_out-liabThreshold-"+str(lT[tidx])+"-tau-"+str(tauValue)+".fam"))
            write_effect_snps(effects_snps_wdenovo, os.path.join(outputpath, "simdigree_out-effects_snps"))
            write_genotype_matrix(generations, os.path.join(outputpath, "simdigree_out-dosage"))
        endbig = time.time()
        print("Time to calculate given tau value was %s" % (endbig-startbig))

    endTHUGE = time.time()
    print("TOTAL TIME TOOK... %s" % (endTHUGE-startTHUGE))
if __name__ == "__main__":
    main()
