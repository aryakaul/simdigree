# Simdigree
#### Arya Kaul

#### 2018-06-28
The plan for this tool is to help Onuralp's analysis with regards to simulating theoretical population datasets. Towards this end I will be developing simdigree. Simdigree will make use of another package, pydigree; however, I will extend the functionality and background of pydigree to the specific use case of Onuralp's work.
**UPDATE** - pydigree no longer being used.

simdigree will have one required argument, a matrix with n rows and m columns. Each of the n rows will delineate an individual simulated by Onuralp through SLiM. Each of the m columns will denote a different SNP and that individual's genotype for that SNP (0,1,2). 

I've created two ways to use simdigree:
1. simdigree generate --> In this option, a user will not supply a sample pedigree to model. Instead, a 'random' pedigree will be generated per their specifications.
    * parents - user says who they want to be the ancestors
    * generations - how many generations will be created? Default will be 4?
    * reproduction - avg. number of kids generated per parental pair
2. simdigree pedigree --> In this option, user supplies a ped file to be modeled. 
    * PED file - contains information about familial relationships, and phenotype of individuals

I spent the rest of the day reading through the documentation for `pydigree`. I am now thoroughly convinced that it makes no sense. I think the simplest way forward would be to make my own set of pedigree functions/classes to handle what I would like to do. I learned some cool things that Mr. Shicks implemented; however, I feel the code is too much of a black box for what I would like it to do. I will start formulating the requisite structure tomorrow.

#### 2018-06-29
I've created a couple functions that will be useful. Currently I'm working on the `generate` functionality. The current workflow is as follows:
1. Read in SNP matrix. Generate `Person` objects for each row in the SNP matrix
2. Read in the effect column vector. 
3. Compute the phenotype of each individual by the matrix multiplication of SNP matrix vs effect column vector
4. Use liability threshold to determine which individuals are classified as 'Affected'

For the `generate` functionality:
1. Pick two random founders from entire population OR use user specified founders
2. Generate random number of children
3. For each generation desired:
    * Loop through current children in generation
        * Marry the current child?
            * Randomly generate some children

To this end I've generated the current structure of the simdigree. Read the `README` for more info.

#### 2018-07-02
The `generate` function is complete; however, I am in the process of developing a read in that takes in a VCF file, as opposed to, a effect vector and a genotype matrix. Additionally, I'll take in a tau coefficient, so that I may calculate the effect vector myself with s^tau. Currently on the #TODO block:

1. Add vcf input functionality
2. Add tau input functionality
3. Add in novel mutation simulation
4. Add in recombination events
5. Sample mates wrt to prevalence of affected vs. healthy (i.e. not 50/50 but based on the size of each population)
6. Calculate affected status after the fact using the founder population

Once this is done, I'll begin working on the `pedigree` functionality of `simdigree`
