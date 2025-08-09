## Exploratory data analyses: unsupervised clustering methods

For the practical exercises on unsupervised clustering we will be analysing sequencing data from three unknown ancient canids. You get to pick one of them and throughout the day you'll characterise the ancestry of these canids and find their identity.

We will start with BAM files (like the one we created yesterday), that has sequencing data mapped to the wolf reference genome <sup>0</sup>, and compare it with a reference dataset of present-day and ancient canids using two approaches: a pseudohaploid approach and a genotype-likelihoods approach. 

### Outline 

* Familiarise with the data and file formats
* Pseudo haploid (random read) approach
    + Incorporate a new sample from a BAM file to a SNP dataset by randomly sampling a read per site using ```Frantk``` <sup>1</sup>
    + Use ADMIXTURE <sup>3</sup> to estimate ancestry proportions
    + Create a PCA using ```smartpca``` <sup>4</sup> with and without projection
* Genotype-likelihood (GL) approach 
    + Estimate GL and perform SNP calling using ```ANGSD``` <sup>2</sup>
    + Use ```NGSadmix``` <sup>5</sup> to estimate ancestry proportions
    + Create a PCA using ```pcangsd``` <sup>6</sup> and the GL


-----------------------------------

### Interactive node

We will start by getting an interactive node:

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# first request one CPU using salloc like this:
salloc --partition=cpuqueue --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000 --account=teaching --reservation=aDNA_PHD_course

# once the job has been allocated, you can login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```

### Data

Once you login to you interactive node, we will create a directory for today's exercises (make sure to change "username" for your own directory name):
 
```{bash, eval = FALSE}
# create directory
directoryExpA="/home/_username_/ExploratoryAnalyses/"
mkdir -p $directoryExpA
# go into the directory
cd $directoryExpA
```

Define some paths and file names. We have three unknow samples, pick only one (copy&paste only the section of the sample you selected): 

```{bash, eval = FALSE}
# Sample 1
BAM="/projects/course_1/people/clx746/Data/sample1.bam"
SAMPLENAME="AncientCanid"

# or sample 2 
BAM="/projects/course_1/people/clx746/Data/sample2.bam"
SAMPLENAME="AncientCanid"

# or sample 3
BAM="/projects/course_1/people/clx746/Data/sample3.bam"
SAMPLENAME="AncientCanid"
```

Also define the path to the SNP panel we will be using: 

```{bash, eval = FALSE}
# Dataset with SNP data for wolves and dogs:
SNPbasename="/projects/course_1/people/clx746/Data/wolves_rand"
```
 
### PLINK files refresher

We will be working with `plink` <sup>7</sup> files during very often during the course, so here is a quick refresher on this format. 

PLINK files come in sets of 2 or 3 files and contain the genotype data, information about the SNP's genomic coordinates, and information about the samples. Some common sets are TPED/TFAM, BED/BIM/FAM and MAP/PED files, all of which are interchangeable using ```plink```. We will be working mostly with BED/BIM/FAM files, more information about the PLINK files and ```plink``` can be found [here](https://www.cog-genomics.org/plink/2.0/input).

The `BIM` file has the information about the SNP's genomic coordinates and allelic variants (one line per SNP):
```{bash, eval = FALSE}
head $SNPbasename".bim" 
```
```
scaffold_0	scaffold_0_9493	0	9493	A	C
scaffold_0	scaffold_0_23098	0	23098	A	C
scaffold_0	scaffold_0_43358	0	43358	A	C
scaffold_0	scaffold_0_53358	0	53358	T	G
scaffold_0	scaffold_0_56736	0	56736	A	T
scaffold_0	scaffold_0_71165	0	71165	A	T
scaffold_0	scaffold_0_96026	0	96026	A	C
scaffold_0	scaffold_0_102434	0	102434	T	G
scaffold_0	scaffold_0_111381	0	111381	T	A
scaffold_0	scaffold_0_113242	0	113242	A	T
```

The `FAM` file has information about the samples (one line per sample):
```{bash, eval = FALSE}
head $SNPbasename".fam"
```
```
Wolf_Alaska Wolf_Alaska 0 0 0 1
Dog_AlaskanHusky Dog_AlaskanHusky 0 0 0 1
Dog_AlaskanM Dog_AlaskanM 0 0 0 1
Wolf_AtlanticCoast Wolf_AtlanticCoast 0 0 0 1
Wolf_Yana Wolf_Yana 0 0 0 1
Wolf_BungeToll Wolf_BungeToll 0 0 0 1
Wolf_Tirekhtyakh Wolf_Tirekhtyakh 0 0 0 1
Wolf_UlakhanSular Wolf_UlakhanSular 0 0 0 1
Wolf_Chinese2 Wolf_Chinese2 0 0 0 1
Wolf_Chinese Wolf_Chinese 0 0 0 1
```

The `BED` file is a compressed file that contains the genotype data.  
```{bash, eval = FALSE}
ls $SNPbasename".bed"  
```

For the pseudohaploid approach we will be working with these plink files:

```{bash, eval = FALSE}
ls $SNPbasename.*  
```
```
/projects/course_1/people/clx746/Data/wolves_rand.bed
/projects/course_1/people/clx746/Data/wolves_rand.bim
/projects/course_1/people/clx746/Data/wolves_rand.fam
```


<span style="color:purple"> **Question:** </span> How many SNPs and samples are present in this dataset? Consider the format of each of the files.

<details>
<summary> <b>Show answer</b> </summary>

```{bash, eval = FALSE}
# You can check the number of samples by counting the lines of the FAM file:
wc -l $SNPbasename.fam

# and the number of SNPs by counting the lines in the BIM file:
wc -l $SNPbasename.bim
```
These files have pseudohaploid SNP data for 20 dogs and 29 wolves that will be useful for estimating the ancestry of our unknow ancient canids. For each of these samples we have 98,241 SNPs (only transversion polymorphisms).

</details>

### Checking the depth of coverage of the samples

Before deciding which approach might be better for you sample of interest it might be a good idea to know at what depth was sequenced. We can use `samtools` to count the number of reads in the BAM file: 

```{bash, eval = FALSE}
samtools view -c $BAM

# -q 30 discard reads with mapping quality lower than 30
samtools view -q 30 -c $BAM
```

<span style="color:purple"> **Question:** </span> Is the BAM file already filtered by mapping quality?

We can also check the number of positions covered in our BAM file using `samtools`:

```{bash, eval = FALSE}
samtools depth -q 30 $BAM  |head
```
```
scaffold_0      9429  1
scaffold_0      9430  2
scaffold_0      9431  2
scaffold_0      9432  2
scaffold_0      9433  2
scaffold_0      9434  2
scaffold_0      9435  2
scaffold_0      9436  2
scaffold_0      9437  2
scaffold_0      9438  2
```

<span style="color:purple"> **Question:** </span> How does the output from `samtools depth` looks like? Can you infer what each column shows?

<details>
<summary> <b>Show answer</b> </summary>

```{bash, eval = FALSE}
Column 1: chromosome/scaffold name
Column 2: genomic position
Column 3: depth of coverage
```
</details>

Now let's count how many positions are covered by at least one read:

```{bash, eval = FALSE}
samtools depth -q 30 $BAM  |grep -vP '\t0$' |wc -l 
```

<span style="color:purple"> **Q:** </span> Can you use the output number of positions covered to estimate the **coverage** of your sample? (Assume the wolf genome is 2312944218 bp long).

<details>
<summary> <b>Show answer</b> </summary>

```{bash, eval = FALSE}
samtools depth -q 30 $BAM  |grep -vP '\t0$' |wc -l 
```
```
13406533
# (the number you get will depend on which sample you selected at the begining)
```
Coverage = 13406533/2312944218 = 0.005796306 ~ 0.57% 

In this case it is important to mention that the BAM file you are using only contains reads overlapping with the SNPs of interest, thus the coverage will be underestimated. 
</details>

... and estimate its **depth of coverage**:

```{bash, eval = FALSE}
samtools depth -q 30 $BAM  |grep -vP '\t0$' |cut -f  3 |sort |uniq -c > CoverageHistogram.txt
```

now plot it as a histogram using ```R``` <sup>8</sup>:

```{r, eval = FALSE}
R
a<-read.table("CoverageHistogram.txt", as.is=T)
a<-t(a[order(as.numeric(a[,2])),])
colnames(a)<-a[2,]
pdf("depth_of_cov.pdf")
barplot(a, las=2, xlim=c(1,35), xlab="Depth of coverage", ylab="Frequency", cex.axis=0.7, cex.names=0.7, col="mediumpurple", border=F)
dev.off()
q("no")
```

You can use **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

Example of the **scp** command:

```{bash, eval = FALSE}
scp clx746@mjolnirgate.unicph.domain:/projects/course_1/people/clx746/ExploratoryAnalyses/depth_of_cov.pdf .
```
(remember to change the directory name to your own)

The result should look like this:

![Figure 1. Histogram showing the per site depth of coverage.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/Exercises/depth_of_cov.png)

<span style="color:purple"> **Question:** </span> Now that you know the depth of coverage of your sample, which approach (random read sampling/GL/genotype-calling) do you think is best?

-----------------------------------

### Pseudohaploid approach

To start exploring the ancestry of our mystery sample, we will estimate a PCA using ```smartpca``` <sup>4</sup> and use the model-based clustering  approach implemented in ```ADMIXTURE``` <sup>3</sup> to estimate ancestry proportions. Both programs can be used with pseudo-haploid data and need BED/BIM/FAM files as input.  

<p>&nbsp;</p>
-----------------------------------

#### Adding an ancient sample (BAM) to a reference SNP dataset by randomly sampling a read

The first step will be to incorporate the mystery sample to the reference dataset. We will use `bam2plink.py` tool from `FrAnTK` to randomly sample one allele/read for each of the 98,241 SNPs and merge it with the reference dataset. More information about  `FrAnTK` can be found [here](https://github.com/morenomayar/FrAnTK).

In order to use `bam2plink.py`, we need to have an additional `regions` file that looks like this: 
```{bash, eval = FALSE}
head $SNPbasename"_regions"
```
```
# chr name	    pos-1   pos   allele1 allele2
scaffold_0      9492    9493    A       C
scaffold_0      23097   23098   A       C
scaffold_0      43357   43358   A       C
scaffold_0      53357   53358   T       G
scaffold_0      56735   56736   A       T
scaffold_0      71164   71165   A       T
scaffold_0      96025   96026   A       C
scaffold_0      102433  102434  T       G
scaffold_0      111380  111381  T       A
scaffold_0      113241  113242  A       T
```
We already have this file, so we don't need to create it. 

Now we can use `bam2plink.py` to randomly sample a read per site from the BAM file of our ancient sample:

```{bash, eval = FALSE}
# first load frAnTK module
module load frantk/20220523

# and then run bam2plink.py
bam2plink.py bamfile=${BAM} plinkpref=${SNPbasename} trim=0 MinMQ=30 MinBQ=20 indname=${SAMPLENAME} popname=${SAMPLENAME}
```

It will take a moment, meanwhile let's take a look at the parameters we are using:

```
bamfile       the BAM file for our ancient sample
plinkpref     the base name of our plink files
trim          how many bases we want to exclude from each read end (this is very useful for aDNA, as you might want to exclude the end of the reads where there is the most ancient DNA damage. 
MinBQ         minimum mapping quality for the reads (remember our BAM is not filterd yet)
MinMQ         minimum base quality for the bases to be considered
indname       name of our ancient sample
popname       name of the population of our ancient sample (if we have one)
```

Once it is done running,  we will have three new PLINK that contain pseudohaploid data for each of the SNPs in our dataset:
```{bash, eval = FALSE}
ls AncientCanid_AncientCanid_wolves_rand.*
```
```
AncientCanid_AncientCanid_wolves_rand.bed
AncientCanid_AncientCanid_wolves_rand.bim
AncientCanid_AncientCanid_wolves_rand.fam
```

Now we can merge the PLINK files with our ancient sample with the SNP panel using `plink`:

```{bash, eval = FALSE}
#first load plink's module
module load plink/1.9.0

# then run plink 
plink --bfile $SNPbasename --make-bed --bmerge ${SAMPLENAME}_${SAMPLENAME}_wolves_rand --allow-extra-chr --allow-no-sex --out  wolves_mergedTv
```

We now have a new SNP dataset that contains the reference data and our ancient sample:

```{bash, eval = FALSE}
ls wolves_mergedTv*
```
```
wolves_mergedTv.bed
wolves_mergedTv.bim
wolves_mergedTv.fam 
```

We will use these files to run ```ADMIXTURE``` clustering and PCA, and tomorrow for the *f*-statistics. 

 <span style="color:cornflowerblue">**BONUS**</span> If you have left over time, you can try creating a SNP dataset with all three unknown canids. 

-----------------------------------


#### Clustering using ADMIXTURE

Now that we have a SNP panel (our BED/BIM/FAM files) with both the reference data and our unknown ancient canid we can use perform ancestry clustering using `ADMIXTURE`.

**NOTE**: In this case our reference SNP dataset was already pseudo-haploid before we incorporated the ancient sample, however if you start with diploid genotypes and incorporate an ancient sample as pseudohaploid, you'll need to sample one allele for each of the reference samples. All samples need to have the same ploidy before running `ADMIXTURE`. 

Defining some paths and file names:

```{bash, eval = FALSE}
#remember to change here for your username
username="write_your_username"
# files
SNPsDS="/projects/course_1/people/${username}/ExploratoryAnalyses/wolves_mergedTv"
SNPsDSchrs="/projects/course_1/people/${username}/ExploratoryAnalyses/wolves_mergedTv_chrs"
```

**NOTE**: A common problem when working with ```plink``` files and fragmented reference genomes is that different programs require the chromosomes to be numeric or to be limited to a certain amount of chromosomes. For example, `plink` complains when the chromosomes are numbers and are more than 99 so we need to add the 'scaffold_' part to the chromosome name and the flag `--allow-extra-chr`. `ADMIXTURE`, on the other hand, prefers chromosomes to be numeric and but doesn't have a restriction for the number of chromosomes. So we will be changing the chromosome/scaffold names in our BIM files depending on which program we use. If you are working with an organisms that has a good genome assembly, you can skip steps marked with  <span style="color: orchid;"> † </span>

<span style="color: orchid;"> † </span> We will edit the BIM file and create a new set of BED/BIM/FAM files: 

```{bash, eval = FALSE}
cat ${SNPsDS}.bim  |perl -pe 's/scaffold_//g;' > ${SNPsDSchrs}.bim 
cp ${SNPsDS}.bed ${SNPsDSchrs}.bed 
cp ${SNPsDS}.fam ${SNPsDSchrs}.fam 
```

An example of an ```ADMIXTURE``` run looks like this (don't run it):
```
admixture -j1 --seed=$RANDOM wolves_mergedTv_c.bed 2 &> out.log
```
 
What parameters are we providing in this case?
```
-j1                     number of CPUs we want to use
-seed                   a starting seed, in this case $RANDOM is picking a random number
wolves_mergedTv_c.bed   our SNPs dataset
2                       the number of Ks or components into which we want to split our individuals
out.log                 a file to print the STDOUT 
```

A peculiarity of `ADMIXTURE` is that it doesn't allow you to specify the name of the output and instead the name is created from the name of the input file. Given we would like to run several replicates for each value of K we will need to create a directory for each ```ADMIXTURE``` run. 

For example, let's say we want to run 5 replicates for each K and we would like split our samples into 2-4 ancestry components (Ks):
```{bash, eval = FALSE}
# first load admixture's module
module load admixture/1.3.0

# then run admixture:
# the first _for_ loop will determine the number of replicates
for rep in 1 2 3 4 5
do
 # make a directory for the current replicate
 mkdir $rep
 cd $rep
 # the second _for_ loop will determine the K
 for k in 2 3 4
 do
   admixture -j1 --seed=$RANDOM ${SNPsDSchrs}.bed $k &> ${rep}_${k}_out.log
 done 
 cd ../
done
```
(This step will take a few minutes to run, so just wait until it is done...)

<span style="color: purple;"> **Q:** </span> Do you remember why running replicates is important? How many parameters are we asking ADMIXTURE to estimate when we ask for K=3?

<details>
<summary> <b>Show answer</b> </summary>

The more parameters our model has, the most likely it is that we can reach a local maximum. When assuming 3 ancestry components we are estimating allele frequencies for 3 populations for 98,241 SNPs, and admixture proportions for 50 samples for each of the 3 components:
```{bash, eval = FALSE}
K=3
nSNPS=98241
nInds=50
TotalParms=(K*nSNPS)+(nIds*(K-1))
```
</details>

Now let's check the results. 

For each `ADMIXURE` run you will get three files, for example for replicate=1 and K=2:

```{bash, eval = FALSE}
ls -l 1/*2*
```
```
1/1_2_out.log
1/wolves_mergedTv_chrs.2.P
1/wolves_mergedTv_chrs.2.Q
```

You'll have a $.P$ file that contains the estimated allele frequencies for each site and for each of the 2 populations:

```{bash, eval = FALSE}
head -n 5 1/wolves_mergedTv_chrs.2.P
```
```
0.696213 0.999990
0.999990 0.887048
0.709935 0.917736
0.870633 0.980417
0.999990 0.963370
```

a $.Q$ file that contains the estimated ancestry proportions of the 3 populations for each of the individuals:

```{bash, eval = FALSE}
head -n 5 1/wolves_mergedTv_chrs.2.Q
```
```
0.999990 0.000010
0.999990 0.000010
0.999990 0.000010
0.809885 0.190115
0.900579 0.099421
```

Check how many lines the .Q file has, is that expected?
```{bash, eval = FALSE}
wc -l 1/wolves_mergedTv_chrs.2.Q
```

How about the .P file?
```{bash, eval = FALSE}
wc -l 1/wolves_mergedTv_chrs.3.P
```

Finally, we have the log file which contains important information about the run:
```{bash, eval = FALSE}
cat 1/1_3_out.log
```
and at the end, there is a summary of the run:

```{bash, eval = FALSE}
tail -n 9 1/1_3_out.log
```
```
Summary: 
Converged in 20 iterations (15.813 sec)
Loglikelihood: -3286365.627536
Fst divergences between estimated populations: 
        Pop0    Pop1
Pop0
Pop1    0.156   
Pop2    0.228   0.253   
Writing output files.
```
In this case it says it converged after 20 iterations and the final Loglikelihood is -3286365.627536

Now we need to decide which replicate is more likely to have reached a global maximum. Let's check the *log* likelihood of all replicates for K=2:

```{bash, eval = FALSE}
grep "^Loglikelihood" */*_2_out.log
```
```
1/1_2_out.log:Loglikelihood: -3470719.048495
2/2_2_out.log:Loglikelihood: -3470719.048495
3/3_2_out.log:Loglikelihood: -3470719.048496
4/4_2_out.log:Loglikelihood: -3470719.048495
5/5_2_out.log:Loglikelihood: -3470719.048501
```

##### How do we chose the replicate with the best likelihood?

Remember that Likelihoods are the product of probabilities, therefore they become very small numbers very quickly. 

It also means that the best model is one that has a higher likelihood (higher probabilities will produce a higher likelihood).

Because programming languages (and computers) are not good at dealing with numbers with a large number of digits, we normally work with *log* likelihoods. 

Let's plot the likelihood value vs. its log and see how it compares:

```{R, eval = TRUE}
# create a list of likelihoods:
Likelihood<-seq(0, 1, by=0.001)
# calculate the log likelihood:
logLikelihood<-log(Likelihood)

# Let's visualize the numbers
plot(Likelihood, logLikelihood)
```

![Likelihoods.](https://github.com/JazminRM/aDNA_SummerCourse_UCPH/tree/main/Figures/loglikevslike.png)

As you can see, as the likelihood increases the *log* likelihood gets closer to 0. 

**Note** that some programs instead of estimating the log likelihood, they estimate the *-log* likelihood, in which case you will have positive numbers. 

<span style="color: purple;"> **Q:** </span> So, which replicate would you choose for each of the Ks?

Remember you can check all the likelihoods like this:
```{bash, eval = FALSE}
grep "^Loglikelihood" */*_3_out.log
```
```
1/1_3_out.log:Loglikelihood: -3259015.293804
2/2_3_out.log:Loglikelihood: -3282163.829765
3/3_3_out.log:Loglikelihood: -3259037.402389
4/4_3_out.log:Loglikelihood: -3259015.401782
5/5_3_out.log:Loglikelihood: -3259037.160064
```
```{bash, eval = FALSE}
grep "^Loglikelihood" */*_4_out.log
```
```
1/1_4_out.log:Loglikelihood: -3080225.071145
2/2_4_out.log:Loglikelihood: -3117689.851233
3/3_4_out.log:Loglikelihood: -3120583.207948
4/4_4_out.log:Loglikelihood: -3080246.524353
5/5_4_out.log:Loglikelihood: -3080903.640304
```
<p>&nbsp;</p>

<button class="btn btn-primary" button style="background-color:purple; border-color:purple; color:white" data-toggle="collapse" data-target="#BlockName77"> Show/hide solution </button>  
<div id="BlockName77" class="collapse">  

In the example above for admixture with K=3, the replicate with the best *log* likelihood is replicate 1 (-3259015.293804) and for K=4, the best one is also 1 (-3080225.071145)

Which one did you get?
</div>

<p>&nbsp;</p>

Now, let's try plotting all the replicates and see how they look:
(This is an option for ploting the results in R, but you can also check [$pong$](https://github.com/ramachandran-lab/pong), which is a very nice option to plot the results from $ADMIXTURE$)
```{r, eval = FALSE}
R

# Read the FAM file with information about the sample names and a file that contains a category for each sample:
fam<-read.table("wolves_mergedTv_chrs.fam", as.is=T)[,1]
info<-read.table("/home/ec2-user/Data/SNPs/wolves_rand_tv_info.txt", as.is=T, sep="\t")

# Find the category for each of the samples in the FAM file
category<-NULL
for(i in 1:length(fam)){
 if(sum(info[,1]==fam[i])>0){
  category<-c(category, info[info[,1]==fam[i],2])
 }else{
  category<-c(category, "X sample") 
 }
}

pdf("admixture_K3.pdf", width=7, height=9)
par(mar=c(7,4,1.5,4), mfrow=c(5,1))
for(i in 1:5){
 a<-read.table(paste0( i, "/wolves_mergedTv_chrs.3.Q"))
 barplot(t(a[order(category),]), names=fam[order(category)], las=2, col=c("lightsalmon1", "#F3E96B", "#6465A5"), cex.names=0.75, border=F)
}
dev.off()


pdf("admixture_K4.pdf", width=7, height=9)
par(mar=c(7,4,1.5,4), mfrow=c(5,1))
for(i in 1:5){
 a<-read.table(paste0( i, "/wolves_mergedTv_chrs.4.Q"))
 barplot(t(a[order(category),]), names=fam[order(category)], las=2, col=c("lightsalmon1", "purple4", "#6465A5", "#F3E96B"), cex.names=0.75, border=F)
}
dev.off()
q("no")
```


Remember you can download the plot we just created to your local computer using **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). In this case the **scp** command would look like this:

```{bash, eval = FALSE}
scp -i apgc-2021-key.pem.txt ec2-user@3.249.84.19:/home/$username/exploratoryA/admixture_K3.pdf .
scp -i apgc-2021-key.pem.txt ec2-user@3.249.84.19:/home/$username/exploratoryA/admixture_K4.pdf .
```

(use your own username)

<span style="color: purple;"> **Q:** </span> Take a look at the five replicates for K=3, does the results vary? 

Identify the replicate with the best log likelihood, is it different from the others? 

<center>

![Figure 2. Admixture results using $ADMIXTURE$ assuming 3 ancestry components (K=3). The results for the mystery sample are hidden. Different colors indicate the inferred ancestry components and the proportion of each corresponds to the estimated admixture proportions.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/IntroPalaeogenomics2023/Figures/admixture_K3.png){width=80%}

</center>

Now take a look at the five replicates for K=4 and repeat the same exercise. 

<center>

![Figure 3. Admixture results using $ADMIXTURE$ assuming 4 ancestry components (K=4). The results for the mystery sample are hidden. Different colors indicate the inferred ancestry components and the proportion of each corresponds to the estimated admixture proportions.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/IntroPalaeogenomics2023/Figures/admixture_K4.png){width=80%}

</center>

Finally, let's plot the best replicate for each of the K's:

```{bash, eval = FALSE}
# Make a directory to copy the best replicates:
mkdir BestRep

# select the best replicate:
R
options(digits=20)
for(k in 3:4){
	logs<-read.table(pipe(paste0("grep '^Loglikelihood' */*_", k, "_out.log")), as.is=T, sep=" ")
	ll<-paste0("cp ", unlist(lapply(strsplit(logs[which.max(logs[,2]),1], "/"), "[[", 1)), "/wolves_mergedTv_chrs.", k, ".* BestRep/")
	system(ll)
}
q("no")
```
And plot the results:

```{r, eval = FALSE}
R

# Read the FAM file with information about the sample names and a file that contains a category for each sample:
fam<-read.table("wolves_mergedTv_chrs.fam", as.is=T)[,1]
info<-read.table("/home/ec2-user/Data/SNPs/wolves_rand_tv_info.txt", as.is=T, sep="\t")

# Find the category for each of the samples in the FAM file
category<-NULL
for(i in 1:length(fam)){
 if(sum(info[,1]==fam[i])>0){
  category<-c(category, info[info[,1]==fam[i],2])
 }else{
  category<-c(category, "X sample") 
 }
}

pdf("admixture_bestLL.pdf", width=10, height=5.5)
par(mar=c(7,4,1.5,4), mfrow=c(2,1))
a<-read.table("BestRep/wolves_mergedTv_chrs.3.Q", as.is=T)
barplot(t(a[order(category),]), names=fam[order(category)], las=2, col=c("lightsalmon1", "#F3E96B", "#6465A5"), cex.names=0.75, border=F)

a<-read.table("BestRep/wolves_mergedTv_chrs.4.Q", as.is=T)
barplot(t(a[order(category),]), names=fam[order(category)], las=2, col=c("lightsalmon1", "purple4", "#6465A5", "#F3E96B"), cex.names=0.75, border=F)
dev.off()

q("no")
```

![Figure 4. Admixture results using $ADMIXTURE$ assuming 3 and 4 ancestry components showing the replicate with the best loglikelihood. Different colors indicate the inferred ancestry components and the proportion of each corresponds to the estimated admixture proportions.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/IntroPalaeogenomics2023/Figures/admixture_bestLL.png){width=80%}

<span style="color: purple;"> **Q:** </span> What can we say about the mystery sample? 

---------------------------------------------------

**Estimating cross-validation errors for the best ADMIXTURE replicate** <span style="color: cornflowerblue;"> BONUS </span>

One option to identify the K (number of clusters) that best explains your data is to perform a cross-validation (CV) for each K and identify the one with the smallest CV error. 

The best way to do that is once you identified the *seed* that took you to the best likelihood for each replicate. 

For example, we saw before that the replicate with the best likelihood for K=3 was replicate **1** (*log likelihood* = -3276680.135746):
```{bash, eval = FALSE}
grep "^Loglikelihood" */*_3_out.log
```
```
1/1_3_out.log:Loglikelihood: -3259015.293804
2/2_3_out.log:Loglikelihood: -3282163.829765
3/3_3_out.log:Loglikelihood: -3259037.402389
4/4_3_out.log:Loglikelihood: -3259015.401782
5/5_3_out.log:Loglikelihood: -3259037.160064
```
and the seed that took you to that loglikelihood is **2556**:
```{bash, eval = FALSE}
grep "Random seed:" 1/1_3_out.log
```
```
Random seed: 24952
```
We can do the same for K=4 and we would get the the replicate with the best *log likelihood* (-3080225.071145) is replicate **1** as well and the seed for that replicate was **14034**.

If we want to perform CV we can do that by re-running ADMIXTURE for each K starting at the seeds we identified:

```{bash, eval = FALSE}
username="write_your_username"

# Let's define our variables again:
SNPsDSchrs="/home/$username/exploratoryA/wolves_mergedTv_chrs"

# and go to your working directory:
cd /home/$username/exploratoryA

# Run cross-validation error for K=3
admixture -j2  --cv --seed=24952 $SNPsDSchrs".bed" 3 &> bestR_3_out.log

# and for K=4
admixture -j2 --cv  --seed=14034 $SNPsDSchrs".bed" 4 &> bestR_4_out.log

```

Once this is done, you will have the cross-validation error in the $.log$ file:
```{bash, eval = FALSE}
grep CV bestR_*_out.log
```
```
bestR_3_out.log:CV error (K=3): 1.13559
bestR_4_out.log:CV error (K=4): 1.24349
```

Theoretically, the K that best explains your data will be the one with the smallest CV-error. Ideally, your CV-error would look like this:

<center>
![Figure 5. Example of a cross-validation error plot for K=2..12](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/IntroPalaeogenomics2022/E1_2/CV_error_CGG6.png){width=50%}
</center>
In the example above, K=9 is the one that best fits the data. 

<p>&nbsp;</p>

-----------------------------------

##### Principal Component Analysis with $smartPCA$

Now we will estimate a PCA using the same BED/BIM/FAM files we used for $ADMIXTURE$

Once again, let's start by defining some directories, paths and file names:

```{bash, eval = FALSE}
username="write_your_username"

# files:
SNPDS="/home/$username/exploratoryA/wolves_mergedTv"
SNPDSnum="/home/$username/exploratoryA/wolves_mergedTv_numChrs"
```

**NOTE** ```smartpca``` takes different types of input files, it can directly take your BED/BIM/FAM files or you can transform them into EIGENSTRAT format and then use them to run ```smartpca```. However, a complication  with ```smartpca``` (given our fragmented reference genome) is that it requires that the chromosome names are numeric but only allows up to 99 chromosomes (at least up to end of last year). So in order for us to use ```smartpca``` we will edit the BIM file and change the chromosome name and positions so that:

1. it has a numeric chromosome name, and
2. the distance between SNPs deesn't change

<span style="color: orchid;"> † </span> We will use the following ```R-script``` to do that:
(You can just copy&paste this in your terminal, and if you get an error try to copy&paste some lines at a time.)
```{r, eval = FALSE}
Rscript /home/ec2-user/Software/renameChrs.R wolves_mergedTv.bim wolves_mergedTv_numChrs.bim
```

Take a look at the old and new BIM files and check that they look as we want to:

```{bash, eval = FALSE}
# old one:
head wolves_mergedTv.bim
```
```
scaffold_0	scaffold_0_9493	0	9493	A	C
scaffold_0	scaffold_0_23098	0	23098	A	C
scaffold_0	scaffold_0_43358	0	43358	A	C
scaffold_0	scaffold_0_53358	0	53358	T	G
scaffold_0	scaffold_0_56736	0	56736	A	T
scaffold_0	scaffold_0_71165	0	71165	A	T
scaffold_0	scaffold_0_96026	0	96026	A	C
scaffold_0	scaffold_0_102434	0	102434	T	G
scaffold_0	scaffold_0_111381	0	111381	T	A
scaffold_0	scaffold_0_113242	0	113242	A	T
```

```{bash, eval = FALSE}
# new one:
head wolves_mergedTv_numChrs.bim
```
```
1 scaffold_0_9493 0 1 A C
1 scaffold_0_23098 0 13606 A C
1 scaffold_0_43358 0 33866 A C
1 scaffold_0_53358 0 43866 T G
1 scaffold_0_56736 0 47244 A T
1 scaffold_0_71165 0 61673 A T
1 scaffold_0_96026 0 86534 A C
1 scaffold_0_102434 0 92942 T G
1 scaffold_0_111381 0 101889 T A
1 scaffold_0_113242 0 103750 A T
```

<span style="color: purple;"> **Q:** </span>Do the SNPs kept their original distance? Why do you think it is important to keep the information about the distance between SNPs?

<span style="color: orchid;"> † </span> Then we will copy the FAM and BED file without editing them: 

```{bash, eval = FALSE}
cp $SNPDS".bed" $SNPDSnum".bed"
cp $SNPDS".fam" $SNPDSnum".fam"
```

We will start by changing the PLINK files into EIGENSTRAT (which is the format that $smartpca$ needs). To do so,we need to create a parameters file that contains the following:

```
genotypename:       wolves_mergedTv_numChrs.bed
snpname:            wolves_mergedTv_numChrs.bim
indivname:          wolves_mergedTv_numChrs.fam
outputformat:       EIGENSTRAT
genotypeoutname:    wolvesTv.eigenstratgeno
snpoutname:	        wolvesTv.snp
indivoutname:       wolvesTv.ind
familynames:        YES
pordercheck:        NO
```

Take a look at the options. Some important options are:

The *genotypename*, *snpname* and *indivname* are input files (BED, BIM, FAM respectively). 

The *outputformat* indicates we want the output in EIGENSTRAT format. 

The *genotypeoutname*, *snpoutname*, and *indivoutname* are the names for the output files. 

And *familynames* allows you to keep the full names of you samples in the output. 

Create the file like this:

```{bash, eval = FALSE}
echo "genotypename: "$SNPDSnum".bed
snpname: "$SNPDSnum".bim
indivname: "$SNPDSnum".fam
outputformat: EIGENSTRAT
genotypeoutname: wolvesTv.eigenstratgeno
snpoutname: wolvesTv.snp
indivoutname: wolvesTv.ind
familynames: YES
pordercheck: NO" > wolves.par
```

Do you know of other ways to create a file in a server? 

Now let's used ```convertf``` to create our EIGENSTRAT files:

```{bash, eval = FALSE}
convertf -p wolves.par
```

**Note**: ```convertf```, ```smartpca``` and some programs of ```admixtools``` require that the samples names are shorter than 39/2 characters, so if you ever get an error saying: *idnames too long*, just change the names in the FAM file for a shorter version

After this we will have three new files:
```{bash, eval = FALSE}
ls wolvesTv.*
```
```
wolvesTv.eigenstratgeno
wolvesTv.ind
wolvesTv.snp
```

Check how the EIGENSTRAT format looks like. You'll have three files that contain information about the genotypes:

```{bash, eval = FALSE}
# SNPs
head wolvesTv.eigenstratgeno
```
```
22900900029000022009000099000090099000900009009000
00000900000000000000000002000000000000000009200299
00200902029000029029000900000002000002000009000090
00000900000000220000000000900020099000000009000099
00900000000000000000000020000000000000090000000090
00090000000000000000222000900900000000000029000000
00900090000000220000020000099000000220009209000009
20002222209020002202002000900000099000000020000290
00000000000000000000002020020202299220020209200002
00000000000000000000000022020202200220200209020209
```

, information about the SNPs:

```{bash, eval = FALSE}
# SNPs info
head  wolvesTv.snp
```
```
     scaffold_0_9493     1        0.000000               1 A C
    scaffold_0_23098     1        0.000136           13606 A C
    scaffold_0_43358     1        0.000339           33866 A C
    scaffold_0_53358     1        0.000439           43866 T G
    scaffold_0_56736     1        0.000472           47244 A T
    scaffold_0_71165     1        0.000617           61673 A T
    scaffold_0_96026     1        0.000865           86534 A C
   scaffold_0_102434     1        0.000929           92942 T G
   scaffold_0_111381     1        0.001019          101889 T A
   scaffold_0_113242     1        0.001038          103750 A T

```

and information about the individuals:

```{bash, eval = FALSE}
# individuals
head wolvesTv.ind 
```
```
Dog_AlaskanHusky:Dog_AlaskanHusky U    Control
Dog_AlaskanM:Dog_AlaskanM U    Control
 Dog_Dingo:Dog_Dingo U    Control
 Dog_Gansu:Dog_Gansu U    Control
 Dog_GMums:Dog_GMums U    Control
Dog_Greenland:Dog_Greenland U    Control
       Dog_GS:Dog_GS U    Control
Dog_Guizhou:Dog_Guizhou U    Control
 Dog_Hebei:Dog_Hebei U    Control
       Dog_ID:Dog_ID U    Control
```


Now we will use these files to run ```smartpca```. Again, we will first create a parameters file that should look like this:

```
genotypename:   wolvesTv.eigenstratgeno
snpname:        wolvesTv.snp
indivname:      wolvesTv.ind
evecoutname:    wolvesTv.evec
evaloutname:    wolvesTv.eval
familynames:    YES
numoutevec:     4
numthreads:     1
pordercheck:    NO
lsqproject:     NO
```
Take a look at the addtitional options:

The *genotypename*, *snpname* and *indivname* in this case are the EIGENSTRAT files we just created. 

The *evaloutname* and *evaloutname* are the names for the output files. 

*numoutevec* indicates the number of PCs that you want printed in the output.

*lsqproject* indicates whether we want to project the low coverage samples or not. For now we will not use this parameter, but we will try it later! **IMPORTANT**: It seems that in the newest version of ```smartpca``` *lsqproject: NO* will not longer be supported, so check which version you are using and its specific documentation. 

Create the file like this:
```{bash, eval = FALSE}
echo "genotypename: wolvesTv.eigenstratgeno
snpname: wolvesTv.snp
indivname: wolvesTv.ind
evecoutname: wolves_mergedTv.evec
evaloutname: wolves_mergedTv.eval
familynames: YES
numoutevec: 4
numthreads: 1
pordercheck: NO
lsqproject: NO" > wolves_pca.par
```

Once you have this file, you can run ```smartpca``` like this:

```{bash, eval = FALSE}
smartpca -p wolves_pca.par
```

You should have two new files now:
```{bash, eval = FALSE}
ls wolves_mergedTv.e*
```
```
wolves_mergedTv.evec
wolves_mergedTv.eval
```

The $.evec$ file contains sample names, followed by the coordinates for the samples in the first 4 PCs (or the number o PCs you indicatd in the parameters file):

```{bash, eval = FALSE}
head -n 5 wolves_mergedTv.evec |column -t
```
```
#eigvals:                          4.176    2.107    1.809    1.727    
Dog_AlaskanHusky:Dog_AlaskanHusky  -0.1875  0.0348   0.0063   -0.0280  Control
Dog_AlaskanM:Dog_AlaskanM          -0.1994  0.0635   -0.0598  0.0725   Control
Dog_Dingo:Dog_Dingo                -0.0896  -0.0217  0.0191   -0.1123  Control
Dog_Gansu:Dog_Gansu                -0.1242  -0.0653  0.0943   -0.1685  Control
```

And the $.eval$ file contains the variance explained in each of the PCs: 

```{bash, eval = FALSE}
head -n 5 wolves_mergedTv.eval 
```
```
    4.334354
    2.105549
    1.766116
    1.663662
    1.539779
```

Now let's plot the results, we will use the same file we used for ```ADMXITURE``` that contains the categories of the samples to assign colors:

```{r, eval = FALSE}
R

library(ggplot2)

evec<-read.table("wolves_mergedTv.evec", as.is=T)
info<-read.table("/home/ec2-user/Data/SNPs/wolves_rand_tv_info.txt", as.is=T, sep="\t")
eval<-as.numeric(readLines("wolves_mergedTv.eval"))

# Here we are taking the first part of the column 1. It looks like this: Alaska1:Alaska1, and we only want this: Alaska1
samps<-sapply(strsplit(evec[,1], ":"), "[[", 1)

# Find the category for each of the samples
category<-NULL
for(i in 1:length(samps)){
 if(sum(info[,1]==samps[i])>0){
  category<-c(category, info[info[,1]==samps[i],2])
 }else{
  category<-c(category, "X sample") 
 }
}

d<-data.frame(category=category, PC1=evec[,2], PC2=evec[,3], sampleid=samps)
# Create a vector of colors (we have 8 categories + our mystery sample), we will manually assign a color to each group:
colvalues <-c("#F3E96B", "coral4", "orange", "lightsalmon1", "#F05837","#6465A5", "darkmagenta", "darkred", "black")
names(colvalues)<-c("Grey wolf America", "Grey wolf Asia Highland", "Grey wolf Middle East", "Grey wolf Europe", "Grey wolf Asia", "Dog", "Dog Arctic", "Ancient grey wolf", "X sample")

p<-ggplot(d, aes(x=PC2, y=PC1))+
geom_vline(xintercept=0, colour="gray", linetype = "longdash", size=0.5)+
geom_hline(yintercept=0, colour="gray", linetype = "longdash", size=0.5)+
geom_point(aes(colour=category), shape=16, size=2.5, alpha=0.8, show.legend = T)+
geom_point(data=d[d$category=="X sample",],  aes(x=PC2, y=PC1), colour="black", shape=18, size=4, alpha=1, show.legend = F)+
scale_colour_manual(values=colvalues)+
ylab(paste0("PC1 (", round(eval[1], 2), "%)"))+
xlab(paste0("PC2 (", round(eval[2], 2), "%)"))+
theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=10), legend.key=element_rect(fill="white", colour="white"), legend.position = "right", axis.title=element_text(size=14), legend.text=element_text(size=10))

pdf("PCA1vsPCA2_smartpca.pdf", width=7, height=5)
print(p)
dev.off()

q("no")
```

Download and take a look at the results (use **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

<center>

![Figure 6. PCA created using smartpca. Note that the mystery sample was not included.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/IntroPalaeogenomics2023/Figures//PCA1vsPCA2_smartpca.png){width=60%}
</center>

<p>&nbsp;</p>


<span style="color: purple;"> **Q:** </span> Where do the mystery sample falls? Is this in agreement with the ADMIXTURE results? What could we infer about its ancestry so far?

<p>&nbsp;</p>

**Projecting samples ** 

A useful feature of ```smartpca``` is that it allows us to 'project' samples into a PCA space build with a subset of the individuals. For example, if your ancient sample is low coverage and high in aDNA error (which can potentially cause problems), you can build a PCA with the reference data and project you ancient sample on top. 

There are two options to project the samples:

1. changing the parameter *lsqproject* to *yes*, without specifying which individuals should be used to estimate the PCs, this way $smartpca$ will only project low coverage samples, 
2. changing the parameter *lsqproject* to *yes* and providing a list (*poplistname*) of samples we want to use to estimate the PCs, this way samples not included in the list will be projected

We will use the second option since, for the purpose of the exercise, we want to project samples that are not necessarily low coverage. 

First will editing the *.ind* file from our EIGENSTRAT files so that it has the sample name in the last column:

```{r, eval=FALSE}
R
a<-read.table("wolvesTv.ind")
p<-as.character(a[,1])
p<-unlist(lapply(strsplit(p, ":"), "[[", 1))
a[,3]<-p
write.table(a, file="wolvesTv.ind", sep=" ", quote=F, row.names=F, col.names=F)
q("no")
```

Then, we will create a list with the populations/samples that we want to use to estimate the PCs (the samples we DO NOT want to project), for this exercise we will project all the Arctic dogs.

```{bash, eval=FALSE}
cut -f 3 -d " " wolvesTv.ind | sort | uniq  |grep -v Dog_GMums |grep -v Dog_Greenland |grep -v Dog_Ilulissat |grep -v Dog_Qaanaaq |grep -v Dog_Sisimiut |grep -v Dog_Tasiilaq |grep -v Dog_AlaskanHusky |grep -v Dog_AlaskanM |grep -v Dog_SibHusk  > wolvesTv_estimatePCs.txt
```

Create the parameters file and include the list we just created in the option *poplistname*:

```{bash, eval = FALSE}
echo "genotypename: wolvesTv.eigenstratgeno
snpname: wolvesTv.snp
indivname: wolvesTv.ind
evecoutname: wolves_mergedTv_GWprojected.evec
evaloutname: wolves_mergedTv_GDprojected.eval
familynames: YES
numoutevec: 4
numthreads: 1
pordercheck: NO
poplistname: wolvesTv_estimatePCs.txt
lsqproject: YES" > wolves_pca_project.par
```

...run ```smartpca```:

```{bash, eval = FALSE}
smartpca -p wolves_pca_project.par
```

...and plot the results:

```{r, eval = FALSE}
R

library(ggplot2)

evec<-read.table("wolves_mergedTv_GWprojected.evec", as.is=T)
info<-read.table("/home/ec2-user/Data/SNPs/wolves_rand_tv_info.txt", as.is=T, sep="\t")
eval<-as.numeric(readLines("wolves_mergedTv_GDprojected.eval"))

# Here we are taking the first part of the column 1. It looks like this: Alaska1:Alaska1, and we only want this: Alaska1
samps<-sapply(strsplit(evec[,1], ":"), "[[", 1)

# Find the category for each of the samples
category<-NULL
for(i in 1:length(samps)){
 if(sum(info[,1]==samps[i])>0){
  category<-c(category, info[info[,1]==samps[i],2])
 }else{
  category<-c(category, "X sample") 
 }
}

d<-data.frame(category=category, PC1=evec[,2], PC2=evec[,3], sampleid=samps)

# Create a vector of colors (we have 8 categories + our mystery sample), we will manually assign a color to each group:
colvalues <-c("#F3E96B", "coral4", "orange", "lightsalmon1", "#F05837","#6465A5", "darkmagenta", "darkred", "black")
names(colvalues)<-c("Grey wolf America", "Grey wolf Asia Highland", "Grey wolf Middle East", "Grey wolf Europe", "Grey wolf Asia", "Dog", "Dog Arctic", "Ancient grey wolf", "X sample")

p<-ggplot(d, aes(x=PC2, y=PC1))+
geom_vline(xintercept=0, colour="gray", linetype = "longdash", size=0.5)+
geom_hline(yintercept=0, colour="gray", linetype = "longdash", size=0.5)+
geom_point(aes(colour=category), shape=16, size=2.5, alpha=1, show.legend = T)+
scale_colour_manual(values=colvalues)+
ylab(paste0("PC1 (", round(eval[1], 2), "%)"))+
xlab(paste0("PC2 (", round(eval[2], 2), "%)"))+
theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=10), legend.key=element_rect(fill="white", colour="white"), legend.position = "right", axis.title=element_text(size=14), legend.text=element_text(size=10))

pdf("PCA1vsPCA2_smartpca_projected.pdf", width=7, height=5)
print(p)
dev.off()
q("no")
```

Download and take a look at the results (use **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

<span style="color: purple;"> **Q:** </span> How do the results from the projection compare with the previous ones? 
<center>
<p>&nbsp;</p>

![Figure 7. PCA created using smartpca without projection (left) and with projection of the Arctic dogs (right).](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/IntroPalaeogenomics2023/Figures/PCA1vsPCA2_smartpca_projectedvsnonprojected.png)

</center>

<p>&nbsp;</p>

-----------------------------------

#### Genotype-likelihoods (GL)

An alternative to sampling a random allele from the BAM is to estimate genotype-likelihoods (GL). There are several methods that have been implemented with GL. We will try to of them today: PCA and Admixture clustering.

Create a new directory to run these exercises:

```{bash, eval = FALSE}
username="write_your_username"

directoryGL="/home/$username/exploratorya_gl"

mkdir -p $directoryGL

cd $directoryGL
```

<p>&nbsp;</p>

-----------------------------------

##### Estimating GL

In order to estimate GL, you'll need to:

1. have BAM files for both your ancient samples and the reference data or, 
1. if you only have SNPs for the reference data you can create 'fake' GL for these and then incorporate you ancient sample (not always recommended and rarely used)

Here we will assume we have BAM files for our reference data. We will be using a similar dataset to the one we used for the pseudo-haploid approach but with less samples (since the GL estimation is a slow process).

Define some paths and file names (again you can choose to work with mystery sample 1 or 2):

```{bash, eval = FALSE}
# list with the paths to the BAMS for the reference and mystery sample (again you can chose from mystery sample1 or sample2)
BAMLIST="/home/ec2-user/Data/BAMS/BamList1.txt"
# or
BAMLIST="/home/ec2-user/Data/BAMS/BamList2.txt"
```

We are starting with list of BAM files with sequencing data from 6 dogs, 5 Eurasian wolves, 5 American wolves, and 5 ancient Siberian wolves. 

```{bash, eval = FALSE}
cat $BAMLIST
```
```
/home/ec2-user/Data/BAMS/Dog_Dingo.bam
/home/ec2-user/Data/BAMS/Dog_GMums.bam
/home/ec2-user/Data/BAMS/Dog_Ilulissat.bam
/home/ec2-user/Data/BAMS/Dog_Shanxi.bam
/home/ec2-user/Data/BAMS/Dog_TM.bam
/home/ec2-user/Data/BAMS/Dog_Tasiilaq.bam
/home/ec2-user/Data/BAMS/Wolf_Alaska.bam
/home/ec2-user/Data/BAMS/Wolf_AtlanticCoast.bam
/home/ec2-user/Data/BAMS/Wolf_BungeToll.bam
/home/ec2-user/Data/BAMS/Wolf_Chinese.bam
/home/ec2-user/Data/BAMS/Wolf_Ellesmere.bam
/home/ec2-user/Data/BAMS/Wolf_GW.bam
/home/ec2-user/Data/BAMS/Wolf_Iberian.bam
/home/ec2-user/Data/BAMS/Wolf_Israeli.bam
/home/ec2-user/Data/BAMS/Wolf_Portuguese.bam
/home/ec2-user/Data/BAMS/Wolf_Saskatchewan.bam
/home/ec2-user/Data/BAMS/Wolf_Taimyr.bam
/home/ec2-user/Data/BAMS/Wolf_Tirekhtyakh.bam
/home/ec2-user/Data/BAMS/Wolf_UlakhanSular.bam
/home/ec2-user/Data/BAMS/Wolf_VictoriaIsland.bam
/home/ec2-user/Data/BAMS/Wolf_Yana.bam
/home/ec2-user/Data/BAMS/sample1.bam
```

In contrast to the previous approach, in this case we don't have a list of sites that we know are polymorphic in our populations, so we need to do SNP calling while calculating the GL. 

To do so we will run the following:

```{bash, eval = FALSE}
angsd -GL 1 -out wolves_gl -doGlf 2 -doMajorMinor 1 -bam $BAMLIST -minMapQ 30 -minQ 20 -rmTrans 1  -doMaf 2 -SNP_pval 2e-6 -minMaf 0.1
```

Take a look at the parameters:

```
-GL 1            this specifies the model to estimate the genotype likelihoods 1 is samtools model 
-minMapQ 30      use reads with a minimum mapping quality of 30
-minQ 20         use bases with a minimum base quality of 20
-doGlf 2         the output should be coded in nucleotides (A/C/T/G)
-doMajorMinor 1  select the minor allele from the data
-doMaf 2         estimate the minor allele frequency
-SNP_pval        minimum p-value to consider a site variable
-minMaf          minimum minor allele frequency to be considered variable
-rmTrans 1       skip transition sites to reduce aDNA-derived error
```

Some important filters are the *-minMaf* this one together with *-SNP_pval* allows us to do the SNP calling. The *-SNP_pval* allows us to keep sites that are variable, and with *-minMaf* we can chose the minimum allele frequency of the SNPs that we want to keep. What does the 0.1 (1%) means? Consider the number of samples we have in our list. 

After this, we'll have three new files:

```{bash, eval = FALSE}
ls wolves_gl*
```
```
wolves_gl.arg
wolves_gl.beagle.gz
wolves_gl.mafs.gz
```

The *.mafs.gz* file contains the allele frequencies for each of the SNPs, the *.arg*, contains the arguments that we used, and the *.beagle.gz* file contains the GL and looks like this:

```{bash, eval = FALSE}
zcat wolves_gl.beagle.gz |head -n 3
```
```
marker  allele1 allele2 Ind0    Ind0    Ind0    Ind1    Ind1    Ind1    Ind2    Ind2    Ind2    Ind3    Ind3    Ind3    Ind4    Ind4    Ind4    Ind5    Ind5    Ind5    Ind6    Ind6    Ind6    Ind7    Ind7      Ind7    Ind8    Ind8    Ind8    Ind9    Ind9    Ind9    Ind10   Ind10   Ind10   Ind11   Ind11   Ind11   Ind12   Ind12   Ind12   Ind13   Ind13   Ind13   Ind14   Ind14   Ind14   Ind15   Ind15   Ind15     Ind16   Ind16   Ind16   Ind17   Ind17   Ind17   Ind18   Ind18   Ind18   Ind19   Ind19   Ind19   Ind20   Ind20   Ind20
scaffold_0_9458 1       3       0.000211        0.333261        0.666528        0.999992        0.000008        0.000000        0.666580        0.333287        0.000133        0.999996        0.000004 0.000000 0.800002        0.199997        0.000001        0.000000        1.000000        0.000000        0.800003        0.199997        0.000000        0.969698        0.030302        0.000000        0.666580  0.333287        0.000133        0.333333        0.333333        0.333333        0.992249        0.007751        0.000000        0.941178        0.058822        0.000000        0.333333        0.333333  0.333333        0.666528        0.333261        0.000211        0.333333        0.333333        0.333333        0.800002        0.199997        0.000001        0.000000        0.199997        0.800003  0.000133        0.333287        0.666580        0.000265        0.333243        0.666492        0.888891        0.111109        0.000000        0.000265        0.333243        0.666492
scaffold_0_23166        2       0       0.666598        0.333296        0.000106        0.941178        0.058822        0.000000        0.800003        0.199997        0.000000        0.999024        0.000976  0.000000        0.800003        0.199997        0.000000        0.999512        0.000488        0.000000        0.800003        0.199997        0.000000        0.984616        0.015384        0.000000  0.333333        0.333333        0.333333        0.000000        0.001949        0.998051        0.992249        0.007751        0.000000        0.800003        0.199997        0.000000        0.992249  0.007751        0.000000        0.333333        0.333333        0.333333        0.000000        1.000000        0.000000        0.800003        0.199997        0.000000        0.333333        0.333333  0.333333        0.666580        0.333287        0.000133        0.666580        0.333287        0.000133        0.999024        0.000976        0.000000        0.800003        0.199997        0.000000
```

<span style="color: purple;"> **Q:** </span> Check the output, do you remember what each of the columns show?

<p>&nbsp;</p>

-----------------------------------
##### Clustering using NGSadmix and genotype-likelihoods

Now we will use ```NGSadmix``` to do admixture clustering using the GL that we estimated. 

This would be an example of how to run ```NGSadmix``` (don't run it yet):

```
NGSadmix -likes wolves_gl.beagle.gz -outfiles wolves_ngsadmx_k3 -P 1 -K 3  -seed $RANDOM
```

Take a look at the parameters

```
-K       how many population should be assumed, depending on the number of samples you can run it with different values of K.
-seed    the starting seed, we want it to be different for each replicate
-minMaf  NGSadmix allow us to apply a minor allele frequency filter in this step, but since we filtered before we can skip it
```

Similar to what we did for ```ADMIXTURE```, we want to run several replicates starting a different seeds for each K. We can do that with the following: 

```{bash, eval = FALSE}
for rep in 1 2 3 4 5
 do
 for k in 3 4
  do
  NGSadmix -likes wolves_gl.beagle.gz -outfiles wolves_gl_ngs_rep"$rep"_k"$k" -P 1 -K $k  -seed $RANDOM
 done
done
```

In this case, since we can specify the name for the output file, we don't need to write the results in different directories. 

It will take a few minutes to run...

```NGSadmix``` will create several output files, for example for K=3 and replicate 1, we have the following:

```{bash, eval = FALSE}
ls wolves_gl_ngs_rep1_k3*
```
```
wolves_gl_ngs_rep1_k3.fopt.gz  wolves_gl_ngs_rep1_k3.log  wolves_gl_ngs_rep1_k3.qopt
```

...*wolves_gl_ngs_rep1_k3.qopt* contains the estimated admixture proportions for each sample 

```{bash, eval = FALSE}
head wolves_gl_ngs_rep1_k3.qopt
```
```
0.65188264262369011171 0.34811735637630991658 0.00000000099999999998 
0.99999999799999994554 0.00000000099999999977 0.00000000099999999977 
0.99999999799999994554 0.00000000099999999975 0.00000000099999999975 
0.62483463326171062580 0.37516536573828929146 0.00000000099999999998 
0.62512946473543706816 0.37487053426456296013 0.00000000099999999999 
0.99999999800000005656 0.00000000099999999969 0.00000000099999999969 
0.00000000099999999991 0.00000000099999999991 0.99999999800000005656 
0.00000000099999999989 0.00000000099999999989 0.99999999800000005656 
0.00000000100000000000 0.66456811258405623644 0.33543188641594384736 
0.00000000099999999992 0.99999999799999994554 0.00000000099999999992 
```
...*wolves_gl_ngs_rep1_k3.fopt.gz* contains allele frequencies for each of the three populations estimated, 

```{bash, eval = FALSE}
zcat wolves_gl_ngs_rep1_k3.fopt.gz |head 
```
```
0.09509227932802177818 0.19913084223263760597 0.06851712602559092979 
0.13587473466096117169 0.21538007266346331781 0.06664937123814040909 
0.00000095328058355000 0.77657997906593045911 0.00002438732240459886 
0.16309322664041536743 0.00004235824199188685 0.41805533849187698214 
0.63530114719519847899 0.00000106426167362266 0.30263125187683920236 
0.00000002049886683294 0.24957110342551805382 0.51978447185614584658 
0.50506170853743614924 0.99999999900192004798 0.00001186413122165328 
0.00000003565541978002 0.30026330690802444101 0.53627074661784002174 
0.00003284496585082848 0.35475894617005626630 0.57093690854113454680 
0.68500035785955004641 0.23457915533710413380 0.31595806384967023428 
```

...and the *wolves_gl_ngs_rep1_k3.log* file contains information about the run including the log likelihood:

```{bash, eval = FALSE}
cat wolves_gl_ngs_rep1_k3.log 
```

```
Input: lname=wolves_gl.beagle.gz nPop=3, fname=(null) qname=(null) outfiles=wolves_gl_ngs_rep1_k3
Setup: seed=9458 nThreads=1 method=1
Convergence: maxIter=2000 tol=0.000010 tolLike50=0.100000 dymBound=0
Filters: misTol=0.050000 minMaf=0.050000 minLrt=0.000000 minInd=0
Input file has dim: nsites=126546 nind=22
Input file has dim (AFTER filtering): nsites=126546 nind=22
	[ALL done] cpu-time used =  62.88 sec
	[ALL done] walltime used =  62.00 sec
best like=-2205731.883820 after 250 iterations
```

In this case the log likelihood is -2205731.883820. 

Similar to what we did with ```ADMIXTURE```, we can identify the replicate with the best log likelihood for K=3 like this:

```{bash, eval = FALSE}
grep "best like" wolves_gl_ngs_rep*_k3*
```
```
wolves_gl_ngs_rep1_k3.log:best like=-2205731.883820 after 250 iterations
wolves_gl_ngs_rep2_k3.log:best like=-2223283.109840 after 200 iterations
wolves_gl_ngs_rep3_k3.log:best like=-2223281.683961 after 300 iterations
wolves_gl_ngs_rep4_k3.log:best like=-2223281.645070 after 300 iterations
wolves_gl_ngs_rep5_k3.log:best like=-2205727.403859 after 350 iterations
```

Which is the replicate with the best log likelihood in this case?


Now, let's plot the results using ```R```:

```{r, eval = FALSE}
R

# Read the file with the list of BAMS, the GL file has the samples in the same order as this list
# Remember to change the to BamList2.txt or BamList3.txt if necessary!
fam<-readLines("/home/ec2-user/Data/BAMS/BamList1.txt")
fam<-gsub(".bam", "", basename(fam))
info<-read.table("/home/ec2-user/Data/SNPs/wolves_rand_tv_info.txt", as.is=T, sep="\t")

# Find the category for each of the samples in the list
category<-NULL
for(i in 1:length(fam)){
 if(sum(info[,1]==fam[i])>0){
  category<-c(category, info[info[,1]==fam[i],2])
 }else{
  category<-c(category, "X sample") 
 }
}

pdf("ngsadmix_K3.pdf", width=7, height=10)
par(mar=c(7,4,1,4), mfrow=c(5,1))
for(i in 1:5){
 a<-read.table(paste0( "wolves_gl_ngs_rep", i, "_k3.qopt"))
 barplot(t(a[order(category),]), names=fam[order(category)], las=2, col=c("#d1aae6", "#6c5ba5", "#e66e58"), cex.names=0.75, border=F)
}
dev.off()


pdf("ngsadmix_K4.pdf", width=7, height=10)
par(mar=c(7,4,1,4), mfrow=c(5,1))
for(i in 1:5){
 a<-read.table(paste0( "wolves_gl_ngs_rep", i, "_k4.qopt"))
 barplot(t(a[order(category),]), names=fam[order(category)], las=2, col=c("#d1aae6", "#6c5ba5", "#e66e58", "#7f346f"), cex.names=0.75, border=F)
}
dev.off()

q("no")
```

Download and take a look at the results (use **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

Check how the replicates differ, select the one with the best likelihood. 

<span style="color: purple;"> **Q:** </span> Did we obtain similar results to the random read approach? 

<center>

![Admixture results using $NGSadmix$ assuming 3 ancestry components (K=3). The results for the mystery sample are hidden.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/Exercises/ngsadmix_K3.png){width=80%}


![Admixture results using $NGSadmix$ assuming 4 ancestry components (K=4). The results for the mystery sample are hidden.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/Exercises/ngsadmix_K4.png){width=80%}

</center>

<p>&nbsp;</p>

-----------------------------------

##### Principal Component Analysis using PCAngsd

Now we will estimate a PCA using our GL. ```PCAngsd``` allows us to estimate a covariance matrix from GL which we can then use to compute a PCA. 

Let's declare some paths:

```{bash, eval = FALSE}
PCANGSD="/home/ec2-user/Software/pcangsd2/pcangsd/pcangsd/pcangsd.py"
```

We can directly take the same GL that we used for ```NGSadmix``` and run ```pcangsd``` like this:

```{bash, eval = FALSE}
python3 $PCANGSD --beagle wolves_gl.beagle.gz -o wolves_gl_pcangsd
```

After this we will have a covariance matrix (wolves_gl_pcangsd.cov) that we can use for the PCA. 

```{bash, eval = FALSE}
ls wolves_gl_pcangsd*
```
```
wolves_gl_pcangsd.args
wolves_gl_pcangsd.cov
```

```pcangsd``` has other interesting options to estimate admixture proportions, inbreeding and perform a selection scan using PC. You can check all the options [here](http://www.popgen.dk/software/index.php/PCAngsdv2#Overview). 

Now we will read the covariance matrix in ```R``` and estimate PCA using ```R``` *prcomp* <sup>8</sup> function:

```{r, eval = FALSE}
R

library(ggplot2)
# read the matrix generated with PCAngsd
covmat<-read.table("wolves_gl_pcangsd.cov", as.is=T)

# read the file with the list of BAMS, the matrix will be in the same order as this list, make sure to change 1 for 2 if you are working with sample 2.  
samps<-readLines("/home/ec2-user/Data/BAMS/BamList1.txt")
samps<-gsub(".bam", "", basename(samps))

# read the information about the categories
info<-read.table("/home/ec2-user/Data/SNPs/wolves_rand_tv_info.txt", as.is=T, sep="\t")

# Find the category for each of the samples in the list
category<-NULL
for(i in 1:length(samps)){
 if(sum(info[,1]==samps[i])>0){
  category<-c(category, info[info[,1]==samps[i],2])
 }else{
  category<-c(category, "X sample") 
 }
}

# perfom a PCA
pca<-prcomp(covmat)

# Create the labels for the axes
var_exp<-round(pca$sdev^2/sum(pca$sdev^2)*100, 4)

d<-data.frame(category=category, PC1=pca$x[,1], PC2=pca$x[,2], sampleid=samps)

# define colors for the different categories
colvalues <-c("#F3E96B", "coral4", "orange", "lightsalmon1", "#F05837","#6465A5", "darkmagenta", "darkred", "black")
names(colvalues)<-c("Grey wolf America", "Grey wolf Asia Highland", "Grey wolf Middle East", "Grey wolf Europe", "Grey wolf Asia", "Dog", "Dog Arctic", "Ancient grey wolf", "X sample")

p<-ggplot(d, aes(x=PC2, y=PC1))+
geom_vline(xintercept=0, colour="gray", linetype = "longdash", size=0.5)+
geom_hline(yintercept=0, colour="gray", linetype = "longdash", size=0.5)+
geom_point(aes(colour=category), shape=16, size=2.5, alpha=0.8, show.legend = T)+
geom_point(data=d[d$category=="X sample",],  aes(x=PC2, y=PC1), colour="black", shape=18, size=4, alpha=1, show.legend = F)+
scale_colour_manual(values=colvalues)+
ylab(paste0("PC1 (", round(var_exp[1], 2), "%)"))+
xlab(paste0("PC2 (", round(var_exp[2], 2), "%)"))+
theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=10), legend.key=element_rect(fill="white", colour="white"), legend.position = "right", axis.title=element_text(size=14), legend.text=element_text(size=10))

# now we plot pc1 vs pc2
pdf("PCA_gl.pdf", width=6, height=4)
print(p)
dev.off()
q("no")

```

Download the result and take a look at them. 
<p>&nbsp;</p>

<center>

![PCA using $pcangsd$. Note that in this figure there is no mystery sample.](/Users/Jazmin/Dropbox/Desktop/Teaching/TransmittingScience/IntroPalaeogenomics2023/Figures/PCA_gl.png){width=60%}
</center>

<p>&nbsp;</p>

<span style="color: purple;"> **Q:** </span> Are the results similar to the ones obtained from ```smartpca``` (besides that there are less samples)?  

<p>&nbsp;</p>
-----------------------------------

**Congratulations! You made it to the end of the exercise.  **

<p>&nbsp;</p>

### References
<p>&nbsp;</p>

0. Gopalakrishnan, S *et al.* 20147 **The wolf reference genome sequence (*Canis lupus lupus*) and its implications for Canis spp. population genomics.** BMC Genomics 18, 495

1. FrantK

2. Korneliussen, T.S *et al.* 2014. **ANGSD: Analysis of Next Generation Sequencing Data.** BMC Bioinformatics 15, 356.

3. Alexander, D. H, *et al.* 2009. **Fast Model-Based Estimation of Ancestry in Unrelated Individuals.** Genome Res. 19, 1655-1664.

4. Patterson, N *et al.* 2016. **Population Structure and Eigenanalysis**. PLoS genetics 2.12 e190.

5. Skotte, L *et al.* 2013. **Estimating Individual Admixture Proportions from Next Generation Sequencing Data.** Genetics 195, 693–702.

6. Meisner J. and Albrechtsen A. 2018. **Inferring Population Structure and Admixture Proportions in Low-Depth NGS Data**. Genetics 210, 2:719:731 

7. Chang CC *et al.* 2015. **Second-Generation PLINK: Rising to the Challenge of Larger and Richer Datasets**. GigaScience, 4.

8. R Core Team. 2020. **R: A language and environment for statistical computing.** R Foundation for  Statistical Computing, Vienna, Austria. URL: https://www.R-project.org/.

9. RStudio Team. 2019. **RStudio: Integrated Development Environment for R.** RStudio, Inc., Boston, MA. URL: http://www.rstudio.com/.

<p>&nbsp;</p>
