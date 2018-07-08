#!/net/home/akaul/miniconda3/bin/python3.6

import time
import argparse
import sys
import math
import os
from simdigree.io import parse_inputNames, read_vcf, write_fam, read_fam
from simdigree.generate import *
from simdigree.pedigree import recreate_pedigree
from simdigree.person import mating
from simdigree.dummy import DummyPair


def parser_args(args):
    parser = argparse.ArgumentParser(prog="simdigree")
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    parser_pedigree = subparsers.add_parser('pedigree', help="Use a sample pedigree file?")
    parser_pedigree.add_argument('-i', '--inputvcf', help="Path to the vcf file being used", type=str, required=True)
    parser_pedigree.add_argument('-T', '--noLoci', help="Number of loci being used", type=int, required=True)
    parser_pedigree.add_argument('-t', '--tau', help="Tau value for effect calculation. Default is 0.5", type=float, required=False, default=0.5)
    #TODO change this to accept a list of values
    parser_pedigree.add_argument('-l', '--liabilityThreshold', help="Float value(s) representing the liability threshold percentage(s) being used. Default is 0.01", type=float, required=False, nargs = '*')
    parser_pedigree.add_argument('-p', '--fam', help="Path to the ped file to use", type=str, required=True)
    parser_pedigree.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is working directory under /simdigree_output", type=str, required=False, default="./simdigree_output")



    parser_generate = subparsers.add_parser('generate', help="Make a fake family history")
    parser_generate.add_argument('-i', '--inputvcf', help="Path to the vcf file being used", type=str, required=True)
    parser_generate.add_argument('-T', '--noLoci', help="Number of loci being used", type=int, required=True)
    parser_generate.add_argument('-t', '--tau', help="Tau value to use for effect calculation. Default is 0.5", type=float, required=False, default=0.5)
    parser_generate.add_argument('-l', '--liabilityThreshold', help="Float value(s) representing the liability threshold percentage(s) being used. Default is 0.01", type=float, required=False, nargs = '*')
    parser_generate.add_argument('-p1', '--ancestor1', help="Line number of individual who will be one ancestral parent. Default is a randomly selected founder.", type=str, required=False, default=None)
    parser_generate.add_argument('-p2', '--ancestor2', help="Line number of individual who will be the other ancestral parent. Default is a randomly selected founder.", type=str, required=False, default=None)
    parser_generate.add_argument('-g', '--generation', help="How many generations will be created total?", type=int, required=False, default=4)
    parser_generate.add_argument('-r', '--reproduction', help="Given a parental pair, what's the mean number of kids they'll have? (If desired # of gen. is not achieved, I force at least one child)", type=int, required=False, default=3)
    parser_generate.add_argument('-m', '--marriagerate', help="Given an individual, what's the the probability they will marry? (If desired # of gen. is not achieved, I force at least one marriage", type=float, required=False, default=0.5)
    parser_generate.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is in the working directory under '/simdigree_output'", type=str, required=False, default="./simdigree_output")
    return parser.parse_args(args)

def calculate_dosage_matrix(founder_genotype_phase_matrix):
    doses = {0:0, 1:1, 2:2, 3:2}
    x = np.copy(founder_genotype_phase_matrix)
    for k, v in doses.items(): x[founder_genotype_phase_matrix==k] = v
    return x

def calculate_scaling_constant(s, maf, tau, h2=0.5):
    C = np.sqrt(h2 / np.sum(2 * maf * (1 - maf) * ((s ** tau) ** 2)))
    if math.isnan(C):
        print("WARNING - Negative value found during square root for scaling constant...\nC = 1")
        C = 1
    return C

def calculate_liability(X, s, C, tau, h2=0.5):
    b = (s ** tau) * C
    G = X @ b
    G = G.reshape(len(G))
    E = np.sqrt(1. - h2) * np.random.randn(len(G))
    return G + E, b

def main():
    args = parser_args(sys.argv[1:])
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    start = time.time()
    founder_genotype_phase_matrix, selection_coeff, allele_freqs, chrom_num, all_founders = read_vcf(args.inputvcf)
    end = time.time()
    print("Time took to read in vcf file... %s" % (end-start))
    start = time.time()
    founder_dosage = calculate_dosage_matrix(founder_genotype_phase_matrix)
    C = calculate_scaling_constant(selection_coeff, allele_freqs, args.tau)
    effects_people, effects_snps = calculate_liability(founder_dosage, selection_coeff, C, args.tau)
    end = time.time()
    print("Time took to calculate liability of founders... %s" % (end-start))
    listOfThresholds = []
    founder_outliers = []
    if args.liabilityThreshold is None:
        lT = [0.01]
    else:
        lT = args.liabilityThreshold
    for i in lT:
        listOfThresholds.append(np.percentile(effects_people, 100*(1-i)))
        founder_outliers.append(np.argwhere(effects_people > i))


    if args.command == "generate":
        if args.ancestor1: founder1 = parse_inputNames(args.ancestor1)
        else: founder1 = None
        if args.ancestor2: founder2 = parse_inputNames(args.ancestor2)
        else: founder2 = None

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

    else:
        print("Unrecognized sub-command, valid options are {pedigree, generate}. Command given: %s" % args.command)
        sys.exit(2)

    for tidx in range(len(lT)):
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
        C = calculate_scaling_constant(selection_coeff, allele_freq_nonfounder)
        effects_nonfounders, effects_snps_wdenovo = calculate_liability(nonfounder_dosage, selection_coeff_new, C)
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
        write_fam(generations, os.path.join(args.output, "simdigree_out-liabThreshold-"+str(lT[tidx])+".fam"))


if __name__ == "__main__":
    main()
