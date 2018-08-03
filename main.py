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
    parser_pedigree.add_argument('-n', '--nosamples', help="Number of individuals in vcf being used", type=int, required=True)
    parser_pedigree.add_argument('-t', '--tau', help="Tau value for effect calculation. Default is 0.5", type=float, required=False, nargs='*')
    parser_pedigree.add_argument('-l', '--liabilityThreshold', help="Float value(s) representing the liability threshold percentage(s) being used. Default is 0.01", type=float, required=False, nargs = '*')
    parser_pedigree.add_argument('-p', '--fam', help="Path to the ped file to use", type=str, required=True)
    parser_pedigree.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is working directory under /simdigree_output", type=str, required=False, default="./simdigree_output")



    parser_generate = subparsers.add_parser('generate', help="Make a fake family history")
    parser_generate.add_argument('-i', '--inputvcf', help="Path to the vcf file being used", type=str, required=True)
    parser_generate.add_argument('-T', '--noLoci', help="Number of loci being used", type=int, required=True)
    parser_generate.add_argument('-n', '--nosamples', help="Number of individuals in vcf being used", type=int, required=True)
    parser_generate.add_argument('-f', '--founders', help="Number of founders matrix should be subsetted to. Default is 15", type=int, required=False, default=40)
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

    doses = {0:0, 1:1, 2:1, 3:2}
    x = np.copy(founder_genotype_phase_matrix)
    for k, v in doses.items(): x[founder_genotype_phase_matrix==k] = v
    return x

#@profile
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

#@profile
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
    founder_genotype_phase_matrix, selection_coeff, allele_freqs, chrom_num, all_founders = read_vcf(args.inputvcf, subset_number, args.nosamples)
    end = time.time()
    print("Time took to read in vcf file: %s" % (end-start))

    #create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    basepath = args.output


    if args.command == "generate":

        #generate the novel pedigree, check generate.py for more details
        start = time.time()
        generations, selection_coeff_new, GT_MATRIX = generate_pedigree(founder_genotype_phase_matrix, args.reproduction, args.generation, args.marriagerate, all_founders, chrom_num, args.noLoci, selection_coeff)
        end = time.time()
        print("Time took to generate pedigree: %s" % (end-start))

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
        generations, selection_coeff_new, GT_MATRIX = recreate_pedigree(individs, pairs, all_founders, founder_genotype_phase_matrix, args.noLoci, chrom_num, selection_coeff)
        numdenovo = GT_MATRIX.shape[1] - founder_genotype_phase_matrix.shape[1]
        print("Number of denovo mutations added is %s" % (numdenovo))
        end = time.time()
        print("Time took to model pedigree: %s" % (end-start))

    #set default tauValues and liabilityThreshold values to loop through
    if args.tau is None:
        tauValues = [0.5]
    else:
        tauValues = args.tau

    if args.liabilityThreshold is None:
        lT = [0.01]
    else:
        lT = args.liabilityThreshold

    # create founder dosage matrix
    start = time.time()
    founder_dosage = read_vcf_founderliab(args.inputvcf)
    end = time.time()
    print("Time took to get big dosage: %s" % (end-start))

    # create and write dosage matrix
    GEN_DOSAGE_MATRIX = calculate_dosage_matrix(GT_MATRIX)
    write_genotype_matrix(generations, GEN_DOSAGE_MATRIX, os.path.join(basepath, "simdigree_out-dosage"))

    start = time.time()
    #start looping over all tauValues
    for tauValue in tauValues:

        #create output path
        taupath = os.path.join(basepath, "tau-"+str(tauValue)+"/")
        if not os.path.exists(taupath):
            os.makedirs(taupath)
        startbig = time.time()
        print("Calculating all values for tau value of %s" % tauValue)

        #go through founders and determine the threshold and scaling constant for the given tau value
        start = time.time()
        C = calculate_scaling_constant(selection_coeff, allele_freqs, tauValue)
        print("C for founders = %s" % C)
        effects_people, effects_snps = calculate_liability(founder_dosage, selection_coeff, C, tauValue)
        end = time.time()
        print("Time took to calculate liability of founders at this tau: %s" % (end-start))

        # figure out the founders who are affected at each liability threshold
        listOfThresholds = []
        founder_outliers = []
        for i in range(len(lT)):
            listOfThresholds.append(np.percentile(effects_people, 100*(1-lT[i])))
            founder_outliers.append(np.argwhere(effects_people > listOfThresholds[i]))

        # calculate the non-founders who are effected per founder_derived_threshold
        start = time.time()
        EFFECTS_GEN, EFFECTS_SNPS_GEN = calculate_liability(GEN_DOSAGE_MATRIX, selection_coeff_new, C, tauValue)
        end = time.time()
        print("Time took to calculate liability of simulated pedigree: %s" % (end-start))

        # write the column vector describing the effects of each SNP
        write_effect_snps(EFFECTS_SNPS_GEN, os.path.join(taupath, "simdigree_out-effects_snps"))

        # loop through every given liability threshold
        print("Looping through all liability thresholds")
        for tidx in range(len(lT)):

            # create output path for diff. liab. thresholds
            outputpath = os.path.join(taupath, "lb-"+str(lT[tidx])+"/")
            if not os.path.exists(outputpath): os.makedirs(outputpath)

            print("On liability threshold of %s" % lT[tidx])
            print("Number of outliers in the founder population... %s" % founder_outliers[tidx].shape[0])
            founder_derived_threshold = listOfThresholds[tidx]
            print("Founder derived threshold is... %s" % founder_derived_threshold)

            # Determine the affected status of non-founders
            OUTLIERS = np.argwhere(EFFECTS_GEN > founder_derived_threshold)
            print("Number of affected in pedigree is %s out of %s" % (len(OUTLIERS), len(EFFECTS_GEN)))

            for person in generations:
                index = person.gt_matrix_ctr
                if index not in OUTLIERS: person.set_affected(False)
                else:
                    print(person.get_name())
                    person.set_affected(True)

            # write the fam file with information about the affected status of all the individuals
            write_fam(generations, os.path.join(outputpath, "simdigree_out-liabThreshold-"+str(lT[tidx])+"-tau-"+str(tauValue)+".fam"))

        endbig = time.time()
        print("Time to calculate given tau value was: %s" % (endbig-startbig))

    endTHUGE = time.time()
    print("Total time took: %s" % (endTHUGE-startTHUGE))
if __name__ == "__main__":
    main()
