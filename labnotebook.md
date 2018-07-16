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

#### 2018-07-03
Started working on the **1** task identified yesterday, the vcf input functionality. I plan on adapting Onuralp's code for this purpose. One thing I need to also understand greater is Onuralp drawing a distinction between heterozygous individuals (0|1 vs 1|0).

I have successfully completed #1, #3, #4, and #5. I still need to develop the affected status calculation; however, I think I'm on the right track.

#### 2018-07-05
Completed the `generate` function! Output is a FAM file with information about status of each person generated. Next step is to code up the scenario where a hardcoded pedigree is given to me. To this end I'm considering using a new datatype? IDK not 100% sure. I'll start working on the code tomorrow and see how I feel. But looking good so far!

#### 2018-07-06
Finished the `pedigree` function! GUCCI GANG GUCCI GANG GUCCI GANG. I am currently running it on the large vcf provided by Onuralp (500 loci) and testing the generate/pedigree function.

The full runall command is produced below:

```
#!/bin/bash
source /net/lib/python3.3/bin/activate
cd /net/home/akaul/projs/Simdigree/jobs-simdigree
SIMDIGREE=/net/home/akaul/projs/Simdigree/git-simdigree/main.py
jobMaker=/net/home/akaul/projs/Simdigree/jobs-simdigree/jobMaker
vcf=/net/home/akaul/projs/Simdigree/jobs-simdigree/data-simdigree/vcf/A1_500_1.vcf

output=/net/home/akaul/projs/Simdigree/jobs-simdigree/output-simdigree/2018-07-06/PEDIGREE
for FAMFILES in /net/home/akaul/projs/Simdigree/git-simdigree/tests/pedigrees/struct/*; do
    fam=$(basename $FAMFILES)
    IFS='-' read -ra PEDSTRUCTS <<< "$fam"
    name="$PEDSTRUCTS"
    jobFile=/net/home/akaul/projs/Simdigree/jobs-simdigree/output-simdigree/2018-07-06/simdigree.job.$name
    cat $jobMaker > $jobFile
    echo "python $SIMDIGREE pedigree -T 500 -i $vcf -p $FAMFILES -l 0.0005 0.001 0.01 0.05 0.10 -o $output/$name" >> $jobFile
    qsub $jobFile
done

output=/net/home/akaul/projs/Simdigree/jobs-simdigree/output-simdigree/2018-07-06/GENERATE
jobFile=/net/home/akaul/projs/Simdigree/jobs-simdigree/output-simdigree/2018-07-06/simdigree.job.generate
cat $jobMaker > $jobFile
echo "python $SIMDIGREE generate -T 500 -i $vcf -p $FAMFILES -l 0.0005 0.001 0.01 0.05 0.10 -o $output" >> $jobFile
qsub $jobFile 
```
Future things to do:
1. Add Tau value list functionality (give me a shitton of tau values)
2. Make hardcoded pedigrees based on medical cases

#### 2018-07-07
Fixed bug where a denovo mutation in one child would cause a different genotype length with the married partner. I have also rerun the jobs.


#### 2018-07-09
So over the weekend all the jobs ran and subsequently crashed. All were due to memory issues. To rectify this issue I have created a new 'utils' folder and under that folder I have put one `subset-vcf.sh` script. This handy dandy little script reads in a given pedigree file, determines how many founders need to be sampled from it, and then randomly subsets the given (large) vcf into a smaller one. There are further improvements I can make to the `simdigree` code that are fairly straightforward, and I plan on tackling those tomorrow. 

#### 2018-07-10
I am rerunning the jobs from yesterday. Realized a small bug in the initial bash script. Today, I'd like to work on the Tau list functionality, and then tackle memory consumption. If I have time, I would also like to fully comment everything.

I am currently testing the tau value list functionality. I placed it at the bottom of my method call, right before the liability threshold calculation...

I have modified the code from yesterday. Instead of subsetting the vcf beforehand, I'll read it into memory, and then subset it while it's in memory only storing and operating on the subsetted matrix. So far, it appears to be working well; however, we'll see how long that lasts. In addition, I have implemented tau list functionality and it appears to be working properly. I've begun commenting things as the jobs are running in the background.

#### 2018-07-11
All the jobs from yesterday have finished running. With the exception of the 'deep' pedigree, they all completed within ~20 minutes and used ~11GB of memory. I don't understand why the deep pedigree did not finish running, and am gonna run some TESTS to find out why.

So oddly, the deep pedigree worked when I ran it on a smaller vcf swimmingly. I'm gonna try again and see if it was a mysterious 'cluster' issue. IT'S RUNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNING.

AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH.

IT WORKED. idk what the initial issue was, but it is working WELL. 

Commit #1 ---> changed name of saved files to be better, also changed print statements around.

I created the pedigree files from clinical data (not on github cause of obvious hippa violations) and am in the process of running simdigree on them. One thing I noticed flawed with my `pedigree` approach is that I do not allow pedigrees with 'multiple' founders. I'll start working on a fix for that next.

#### 2018-07-12
Added output functionality where I also output the genotype matrix and the snp effects column vector. Also patched some bugs with `pedigree` function and dealing with denovo mutations.

NOW next on list is the multiple founders stuff

#### 2018-07-13 to 2018-07-15

Multiple founders has been coded. I have run 100 trials on each of the pedigree clinical cases, as well as, the pedigree structures provided before. In addition, I generated 100 novel pedigrees.

#### 2018-07-16

Today I will comment EVERYTHING.

In addition to commenting everything, I also fixed a piece of logic. The founder liability threshold will now be determined by the COMPLETE vcf as opposed to the subset we extracted earlier.

I'm running tests to make sure this works.
