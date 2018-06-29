#!/bin/env/python

import sys
import os

def readin_matrix(filein):
    matrix=[]
    with open(filein, 'r') as f:
        for lines in f:
            line = lines.rstrip().split()
            matrix.append(line)
    return matrix
