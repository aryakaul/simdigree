#!/bin/env/python

import argparse
import sys
import os

def parser_args(args):
    parser = argparse.ArgumentParser(prog="simdigree")
    parser.add_argument('-i', '--input', help="Path to the matrix file being read", type=argparse.FileType('r'), required=True)
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_pedigree = subparsers.add_parser('pedigree', help="Use a sample pedigree file?")
    parser_pedigree.add_argument('-p', '--ped', help="Path to the ped file to use", type=argparse.FileType('r'), required=True)

    parser_generate = subparsers.add_parser('generate', help="Make a fake family history")
    parser_generate.add_argument('-p1', '--ancestor1', help="Line number of individual who will be one ancestral parent", type=str, required=False)
    parser_generate.add_argument('-p2', '--ancestor2', help="Line number of individual who will be the other ancestral parent", type=str, required=False)
    parser_generate.add_argument('-g', '--generation', help="How many generations will they create?", type=int, required=False)
    parser_generate.add_argument('-r', '--reproduction', help="Given a parental pair, what's the mean number of kids they'll have?", type=int, required=False)

    return parser.parse_args(args)

def main():
    args = parser_args(sys.argv[1:])

if __name__ == "__main__":
    main()
