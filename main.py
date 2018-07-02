#!/net/home/akaul/miniconda3/bin/python3.6

import argparse
import sys
import os
from simdigree.io import determine_healthy_affected, write_pedigree
from simdigree.generate import *
from simdigree.person import mating


def parser_args(args):
    parser = argparse.ArgumentParser(prog="simdigree")
    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    parser_pedigree = subparsers.add_parser('pedigree', help="Use a sample pedigree file?")
    parser_pedigree.add_argument('-i', '--inputvcf', help="Path to the vcf file being used", type=str, required=True)
    parser_pedigree.add_argument('-t', '--tau', help="Tau value to use to raise the S coeff", type=float, required=False, default=0.5)
    #TODO put actually right thing here ^
    parser_pedigree.add_argument('-l', '--liabilityThreshold', help="Float value representing the liability threshold percentage being used. Default is 0.01", type=float, required=False, default=0.01)
    parser_pedigree.add_argument('-p', '--ped', help="Path to the ped file to use", type=str, required=True)
    parser_pedigree.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to", type=str, required=False, default="./simdigree_output")

    parser_generate = subparsers.add_parser('generate', help="Make a fake family history")
    parser_generate.add_argument('-i', '--inputvcf', help="Path to the vcf file being used", type=str, required=True)
    parser_generate.add_argument('-t', '--tau', help="Tau value to use to raise the S coeff", type=float, required=False, default=0.5)
    #TODO put actually right thing here ^
    parser_generate.add_argument('-l', '--liabilityThreshold', help="Float value representing the liability threshold percentage being used. Default is 0.01", type=float, required=False, default=0.01)
    parser_generate.add_argument('-p1', '--ancestor1', help="Line number of individual who will be one ancestral parent. Default is a randomly selected healthy individual.", type=str, required=False, default=None)
    parser_generate.add_argument('-p2', '--ancestor2', help="Line number of individual who will be the other ancestral parent. Default is a randomly selected sick individual", type=str, required=False, default=None)
    parser_generate.add_argument('-g', '--generation', help="How many generations will be created total?", type=int, required=False, default=4)
    parser_generate.add_argument('-r', '--reproduction', help="Given a parental pair, what's the mean number of kids they'll have?", type=int, required=False, default=3)
    parser_generate.add_argument('-m', '--marriagerate', help="Given an individual, what's the the probability they will marry?", type=float, required=False, default=0.5)
    parser_generate.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to", type=str, required=False, default="./simdigree_output")

    return parser.parse_args(args)

def main():

    args = parser_args(sys.argv[1:])
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    healthy, affected = determine_healthy_affected(args.snpMatrix, args.effectVec, args.liabilityThreshold)
    healthy_shuffle = generate_shuffle(healthy)
    affected_shuffle = generate_shuffle(affected)


    if args.command == "generate":
        generations = generate_pedigree(healthy, affected, healthy_shuffle, affected_shuffle, args.ancestor1, args.ancestor2, args.reproduction, args.generation, args.marriagerate)
        write_pedigree(generations, os.path.join(args.output, "test-pedigree.ped"))
    elif args.command=="pedigree":
        people = readin_matrix(args.input)
        pedigree_people = readin_ped(args.ped)


    else:
        print("Unrecognized sub-command, valid options are {pedigree, generate}. Command given: %s" % args.command)
        sys.exit(2)

if __name__ == "__main__":
    main()
