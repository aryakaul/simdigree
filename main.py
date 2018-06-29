#!/net/home/akaul/miniconda3/bin/python3.6

import argparse
import sys
import os
from simdigree.io import readin_matrix
from simdigree.generate import *
from simdigree.person import mating


def parser_args(args):
    parser = argparse.ArgumentParser(prog="simdigree")
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    parser_pedigree = subparsers.add_parser('pedigree', help="Use a sample pedigree file?")
    parser_pedigree.add_argument('-i', '--input', help="Path to the snp-genotype matrix being used", type=str, required=True)
    parser_pedigree.add_argument('-e', '--effectMatrix', help="Path to the effect matrix file being used", type=str, required=True)
    parser_pedigree.add_argument('-p', '--ped', help="Path to the ped file to use", type=str, required=True)

    parser_generate = subparsers.add_parser('generate', help="Make a fake family history")
    parser_generate.add_argument('-i', '--input', help="Path to the snp-genotype matrix being read", type=str, required=True)
    parser_generate.add_argument('-e', '--effectMatrix', help="Path to the effect matrix file being used", type=str, required=True)
    parser_generate.add_argument('-p1', '--ancestor1', help="Line number of individual who will be one ancestral parent", type=str, required=False)
    parser_generate.add_argument('-p2', '--ancestor2', help="Line number of individual who will be the other ancestral parent", type=str, required=False)
    parser_generate.add_argument('-g', '--generation', help="How many generations will be created total?", type=int, required=False, default=4)
    parser_generate.add_argument('-r', '--reproduction', help="Given a parental pair, what's the mean number of kids they'll have?", type=int, required=False, default=3)
    parser_generate.add_argument('-m', '--marriagerate', help="Given an individual, what's the the probability they will marry?", type=float, required=False, default=0.5)

    return parser.parse_args(args)

def main():
    args = parser_args(sys.argv[1:])
    if args.command == "generate":
        """
        people,numSNPs = readin_matrix(args.input)
        ped_out = []
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
            currGen = newGen
        """

    elif args.command=="pedigree":
        people = readin_matrix(args.input)
        pedigree_people = readin_ped(args.ped)


    else:
        print("Unrecognized sub-command, valid options are {pedigree, generate}. Command given: %s" % args.command)
        sys.exit(2)

if __name__ == "__main__":
    main()
