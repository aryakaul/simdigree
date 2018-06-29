# Simdigree
#### Arya Kaul

#### 2018-06-28
The plan for this tool is to help Onuralp's analysis with regards to simulating theoretical population datasets. Towards this end I will be developing simdigree. Simdigree will make use of another package, pydigree; however, I will extend the functionality and background of pydigree to the specific use case of Onuralp's work.

simdigree will have one required argument, a matrix with n rows and m columns. Each of the n rows will delineate an individual simulated by Onuralp through SLiM. Each of the m columns will denote a different SNP and that individual's genotype for that SNP (0,1,2). 

I've created two ways to use simdigree:
1. simdigree generate --> In this option, a user will not supply a sample pedigree to model. Instead, a 'random' pedigree will be generated per their specifications.
    * parents - user says who they want to be the ancestors
    * generations - how many generations will be created? Default will be 4?
    * reproduction - avg. number of kids generated per parental pair
2. simdigree pedigree --> In this option, user supplies a ped file to be modeled. 
    * PED file - contains information about familial relationships, and phenotype of individuals

I spent the rest of the day reading through the documentation for `pydigree`. I am now thoroughly convinced that it makes no sense. I think the simplest way forward would be to make my own set of pedigree functions/classes to handle what I would like to do. I learned some cool things that Mr. Shicks implemented; however, I feel the code is too much of a black box for what I would like it to do. I will start formulating the requisite structure tomorrow.
