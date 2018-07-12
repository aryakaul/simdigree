# simdigree
A tool to apply SLiM population simulations onto pedigrees, and to generate artificial user-defined pedigrees.

## Installation

Install `simdigree` through `git clone`. Run it by using the `main.py` script found in the root directory and ensure you're using `Python3` and `numpy` is installed!

## Usage
`simdigree` is broken down into two main use cases:
    1. `simdigree pedigree` - Using a pedigree as a template
    2. `simdigree generate` - randomly generate a pedigree

For both of those use cases. `simdigree` has 2 required arguments.
1. A vcf file as outputted by SLiM
2. The number of loci concatenated


### pedigree

This function takes in a pedigree in plink's `fam` format. It then simulates variants on the structure provided. Usage is described below:
```
usage: simdigree pedigree [-h] -i INPUTVCF -T NOLOCI [-t TAU]
                          [-l [LIABILITYTHRESHOLD [LIABILITYTHRESHOLD ...]]]
                          -p FAM [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTVCF, --inputvcf INPUTVCF
                        Path to the vcf file being used
  -T NOLOCI, --noLoci NOLOCI
                        Number of loci being used
  -t TAU, --tau TAU     Tau value for effect calculation. Default is 0.5
  -l [LIABILITYTHRESHOLD [LIABILITYTHRESHOLD ...]], --liabilityThreshold [LIABILITYTHRESHOLD [LIABILITYTHRESHOLD ...]]
                        Float value(s) representing the liability threshold
                        percentage(s) being used. Default is 0.01
  -p FAM, --fam FAM     Path to the ped file to use
  -o OUTPUT, --output OUTPUT
                        Path to the output folder to dump simdigree's output
                        to. Default is working directory under
                        /simdigree_output
```

### generate

Here is the generate function's help message. It stochastically generates a pedigree per the specifications given.

```
usage: simdigree generate [-h] -i INPUTVCF -T NOLOCI [-t TAU]
                          [-l [LIABILITYTHRESHOLD [LIABILITYTHRESHOLD ...]]]
                          [-p1 ANCESTOR1] [-p2 ANCESTOR2] [-g GENERATION]
                          [-r REPRODUCTION] [-m MARRIAGERATE] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTVCF, --inputvcf INPUTVCF
                        Path to the vcf file being used
  -T NOLOCI, --noLoci NOLOCI
                        Number of loci being used
  -t TAU, --tau TAU     Tau value to use for effect calculation. Default is
                        0.5
  -l [LIABILITYTHRESHOLD [LIABILITYTHRESHOLD ...]], --liabilityThreshold [LIABILITYTHRESHOLD [LIABILITYTHRESHOLD ...]]
                        Float value(s) representing the liability threshold
                        percentage(s) being used. Default is 0.01
  -p1 ANCESTOR1, --ancestor1 ANCESTOR1
                        Line number of individual who will be one ancestral
                        parent. Default is a randomly selected founder.
  -p2 ANCESTOR2, --ancestor2 ANCESTOR2
                        Line number of individual who will be the other
                        ancestral parent. Default is a randomly selected
                        founder.
  -g GENERATION, --generation GENERATION
                        How many generations will be created total?
  -r REPRODUCTION, --reproduction REPRODUCTION
                        Given a parental pair, what's the mean number of kids
                        they'll have? (If desired # of gen. is not achieved, I
                        force at least one child)
  -m MARRIAGERATE, --marriagerate MARRIAGERATE
                        Given an individual, what's the the probability they
                        will marry? (If desired # of gen. is not achieved, I
                        force at least one marriage
  -o OUTPUT, --output OUTPUT
                        Path to the output folder to dump simdigree's output
                        to. Default is in the working directory under
                        '/simdigree_output'
```
