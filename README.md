# simdigree
A tool to apply SLiM population simulations onto pedigrees, and to generate artificial user-defined pedigrees.

## Usage
`simdigree` is broken down into two main use cases:
    1. `simdigree pedigree` - Using a pedigree as a template
    2. `simdigree generate` - randomly generate a pedigree

For both of those use cases. `simdigree` has 2 required arguments.

The first is a snp matrix file, an example of which is reproduced below:

|SNP1|SNP2|...|
|:-:|:-:|:-:|
|1|2|...|
|0|2|...|

This file is a matrix with `m` rows and `n` columns where `m` is the number of individuals randomly generated, and `n` is the number of SNPs genotyped. The value in [i][j] denotes the genotype of Person i, at SNP j.

The second is a effect column vector file, an example of which is reproduced below:

|SNP-Effect|
|:--------:|
|0.4523|
|0.232323|
|...|

This file has exactly `n` rows, each denoting the contribution of that rows SNP to the phenotype being studied.

### pedigree

This function takes in a pedigree in plink's `ped` format. It is currently not implemented.  

### generate

Here is the generate function's help message.

```
usage: simdigree generate [-h] -i INPUT -e EFFECTMATRIX [-p1 ANCESTOR1]
                          [-p2 ANCESTOR2] [-g GENERATION] [-r REPRODUCTION]
                          [-m MARRIAGERATE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the matrix file being read
  -e EFFECTMATRIX, --effectMatrix EFFECTMATRIX
                        Path to the snp-genotype matrix file being used
  -p1 ANCESTOR1, --ancestor1 ANCESTOR1
                        Line number of individual who will be one ancestral
                        parent
  -p2 ANCESTOR2, --ancestor2 ANCESTOR2
                        Line number of individual who will be the other
                        ancestral parent
  -g GENERATION, --generation GENERATION
                        How many generations will be created total?
  -r REPRODUCTION, --reproduction REPRODUCTION
                        Given a parental pair, what's the mean number of kids
                        they'll have?
  -m MARRIAGERATE, --marriagerate MARRIAGERATE
                        Given an individual, what's the the probability they
                        will marry?
```
