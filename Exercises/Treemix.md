### Admixture graphs

### Outline

* Admixture graphs using Treemix  <sup>1</sup>
   + Create Treemix input using `FrAnTK` <sup>2</sup>
   + Build a `TreeMix` tree

### Interactive node

Start by getting an interactive node: 

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# request one CPU using salloc like this:
salloc --partition=cpuqueue --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000 --account=teaching --reservation=aDNA_PHD_course

# once the job has been allocated, you can login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```

Go to the ExploratoryAnalyses directory (we will use the same PLINK files we generated earlier):
```{bash, eval = FALSE}
# remember to change ku_username for your username:
cd /projects/course_1/people/ku_username/ExploratoryAnalyses
```

### TreeMix trees

Now we now how our ancient sample looks like in terms of its general genetic ancestry. We know what is the population it is closest to, we also have an idea of potential admixture with other populations, so the next thing we want to do is to build a tree. 

In this section we will build a tree with the reference data and our ancient samples using ```TreeMix``` <sup>1</sup> . 

```TreeMix``` takes as input a file with allele counts and since we want to discard sites with missing data we will select only some relevant and high coverage samples for the tree (two samples from each category) and we will use the golden jackal as outgroup:

(Feel free to include/exclude more samples if you want!)

```{bash, eval=FALSE}
cut -f 1 wolves_rand_ftk_pop |grep -e GoldenJackal -e Dog_GMums -e Dog_Qaanaaq -e Dog_Hebei -e Dog_TM -e Wolf_Israel -e Wolf_UlakhanSular -e Wolf_Qamanirjuaq -e Wolf_NorthBaffin -e Wolf_Portuguese -e AncientCanid -e Dog_SibHusky -e Dog_Greenland -e Dog_AlaskanM -e Dog_AlaskanHusky >  pop2kep
```

Our *pop2kep* file now contains a list of the populations we want to include in the tree. We will use ```FrAnTK```'s ```Freqs2Treemix``` tool to create the input for ```TreeMix```

```{bash, eval=FALSE}
frantk Freqs2Treemix freqpref=wolves_rand_ftk popsofint=pop2kep tmpref=wolves_rand_ftk_treemix
```

Take a look at the parameters we are using with ```Freqs2Treemix```:

```
freqpref    the basename or prefix of our precomputed allele frequencies
popsofint   a list of population that we want to include in the treemix tree
tmpref      prefix for the output file
```

After it is done running we will have four files:
```{bash, eval=FALSE}
ls wolves_rand_ftk_treemix*
```
```
wolves_rand_ftk_treemix_ALL_tm.gz
wolves_rand_ftk_treemix_ALL_tmpos
wolves_rand_ftk_treemix_NT_tm.gz
wolves_rand_ftk_treemix_NT_tmpos
```
The ones ending in *ALL_tm* contain all sites and the ones ending in *NT_tm* contain only the transversion sites. (Since our original files do not have transitions, these will be the same length.)

<span style="color: purple;"> **Optional (but highly encouraged):** </span> Take a look at the output files, check what they contain and how many sites/SNPs we will be using for the ```TreeMix``` analysis.

<details>
<summary> <b>Show answer</b> </summary>

The *wolves_rand_ftk_treemix_ALL_tm.gz* file contains the allele counts and it is the one we will be using for ```TreeMix```:
```{bash, eval=FALSE}
zcat wolves_rand_ftk_treemix_ALL_tm.gz |head
```
```
Dog_AlaskanHusky Dog_AlaskanM Dog_GMums Dog_Greenland Dog_Hebei Dog_Qaanaaq Dog_SibHusky Dog_TM GoldenJackal AncientCanid Wolf_Israel Wolf_NorthBaffin Wolf_Portuguese Wolf_Qamanirjuaq Wolf_UlakhanSular
0,1 0,1 0,1 0,1 0,1 0,1 0,2 0,1 2,0 1,0 0,2 0,1 0,1 0,1 0,1
0,1 0,1 0,1 0,1 0,1 0,1 1,1 0,1 0,2 0,1 0,2 1,0 0,1 0,1 0,1
0,1 0,1 0,1 0,1 0,1 0,1 0,2 0,1 2,0 0,1 0,2 1,0 1,0 0,1 0,1
0,1 1,0 0,1 0,1 0,1 1,0 0,2 1,0 0,1 0,1 0,2 0,1 1,0 0,1 0,1
0,1 0,1 0,1 1,0 0,1 0,1 0,2 0,1 0,1 0,1 0,2 0,1 0,1 0,1 0,1
0,1 0,1 0,1 0,1 0,1 0,1 0,2 0,1 0,3 0,1 0,2 0,1 0,1 0,1 1,0
1,0 0,1 0,1 0,1 0,1 0,1 1,0 0,1 0,2 0,1 0,2 0,1 0,1 0,1 1,0
0,1 0,1 0,1 0,1 0,1 1,0 2,0 1,0 0,2 0,1 0,2 0,1 0,1 0,1 0,1
0,1 0,1 0,1 0,1 0,1 1,0 0,2 0,1 0,2 0,1 0,2 0,1 0,1 0,1 0,1
```
As you can see the first line is the header with the population names. Each column corresponds to a population and each row to a SNP. For each SNP and population we have allele counts, for example for the *Dog_AlaskanHusky* we have *0,1* (0 alleles reference and 1 allele alternative) for the first SNP. 

We can check the number of SNPs by counting the lines:
```{bash, eval=FALSE}
zcat wolves_rand_ftk_treemix_ALL_tm.gz |wc -l
```
```
32761
```
so, that would be 32761-1 = 32760 SNPs that we will be using for the tree. 

The other file *wolves_rand_ftk_treemix_ALL_tmpos* contains information about the SNPs (chromosome, coordinate, snpname, ref and alternative alleles). 

```{bash, eval=FALSE}
head wolves_rand_ftk_treemix_NT_tmpos
```
```
scaffold_0	71165	scaffold_0_71165	A	T
scaffold_0	96026	scaffold_0_96026	A	C
scaffold_0	113242	scaffold_0_113242	A	T
scaffold_0	134733	scaffold_0_134733	G	T
scaffold_0	147209	scaffold_0_147209	T	A
scaffold_0	162559	scaffold_0_162559	G	C
scaffold_0	202531	scaffold_0_202531	G	T
scaffold_0	314745	scaffold_0_314745	A	T
scaffold_0	346152	scaffold_0_346152	T	A
scaffold_0	360221	scaffold_0_360221	T	A
```

</details>

```TreeMix``` allow us to model migration edges on top of our tree to include potential admixture events, so we will test from 0 to 3 possible migration edges (*-m* parameter). Since ```TreeMix``` uses a maximum-likelihood approach (same as ```ADMIXURE```), we want to run ```TreeMix``` several times starting at different seeds and keep the replicate with the best likelihood for each number of migrations (*-m*). 

We can use the following loop in ```bash``` to run 5 replicates for each value of *-m*:

```{bash, eval=FALSE}
mkdir treemixRes
cd treemixRes
for rep in 1 2 3 4 5
do
 for migration in 0 1 2 3
 do
  treemix -j1 -i ../wolves_rand_ftk_treemix_ALL_tm.gz -seed $RANDOM -o wolvesTM_m"$migration"_rep"$rep" -k 500 -noss -global -m $migration -root GoldenJackal
 done 
done
```

Take a look at the parameters:

```
-i      input file
-seed   starting seed, we want a different seed for each run, so we set it to random
-k      number of SNPs for each block, this will depend on the distribution of SNPs in your genome and the expected LD block size 
-m      number of migration edges to be modelled
-root   an individual to root the tree (GoldenJackal)
-o      base name for the output files
```

After TreeMix is done, we will have the following files for each migration and each replicate. For example, for migration 3 and replicate 5:

```{bash, eval=FALSE}
ls wolvesTM_m3_rep5*
```
```
wolvesTM_m3_rep5.cov.gz
wolvesTM_m3_rep5.covse.gz
wolvesTM_m3_rep5.edges.gz
wolvesTM_m3_rep5.llik
wolvesTM_m3_rep5.modelcov.gz
wolvesTM_m3_rep5.treeout.gz
wolvesTM_m3_rep5.vertices.gz
```

In this case the *.llik* that contains the final likelihood. 

For migration 3 replicate 5:
```{bash, eval=FALSE}
cat wolvesTM_m3_rep5.llik
```
```
Starting ln(likelihood) with 0 migration events: 385.131 
Exiting ln(likelihood) with 3 migration events: 436.014 
```

For all the replicates for migration 3:

```{bash, eval=FALSE}
 grep "^Exiting" wolvesTM_m3_*llik
 ```
 ```
wolvesTM_m3_rep1.llik:Exiting ln(likelihood) with 3 migration events: 436.014 
wolvesTM_m3_rep2.llik:Exiting ln(likelihood) with 3 migration events: 436.014 
wolvesTM_m3_rep3.llik:Exiting ln(likelihood) with 3 migration events: 435.144 
wolvesTM_m3_rep4.llik:Exiting ln(likelihood) with 3 migration events: 436.014 
wolvesTM_m3_rep5.llik:Exiting ln(likelihood) with 3 migration events: 436.014 
```

In this case there are four replicates with the same final likelihood: 436.014

Let's evaluate the likelihood of the replicates and choose the best one using ```R``` <sup>5</sup>, we will also plot the results using ```TreeMix``` ```R``` functions. 

We will need a file with a list of the populations, so let's create that one first:

```{bash, eval=FALSE}
zcat ../wolves_rand_ftk_treemix_ALL_tm.gz |head -n 1 |perl -pe 's/ /\n/g;' >poporderWolves
```

Plot the results using ```R```. 

```{r, eval=FALSE}
R

# load Treemix functions
source("/projects/course_1/people/clx746/Data/plotting_funcs.R")

pdf("wolves_treemix.pdf", width=14, height=7)
par(mfrow=c(1,2))

for(m in 0:3){
 llik<-read.table(pipe(paste0("grep 'Exiting' wolvesTM_m", m, "_rep*.llik")), as.is=T)
 reps<-as.numeric(gsub("rep", "", sapply(strsplit(sapply(strsplit(llik[,1], "_"), "[[", 3), "\\."), "[[", 1)))
 bestrep<-reps[which.max(as.numeric(llik[,7]))]
 basena<-paste0("wolvesTM_m", m, "_rep", bestrep[1])
 
 plot_tree(basena)
# par(mar=c(13.1,13.1,7.1,7.1))
 plot_resid(basena, "poporderWolves")
}

dev.off()

q("no")
```

Download the plot we just created to your local computer using **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

Example of the **scp** command:

```{bash, eval = FALSE}
scp clx746@mjolnirgate.unicph.domain:/projects/course_1/people/clx746/Mapping/ExploratoryAnalyses/wolves_treemix.pdf .
```
(remember to use your own username)

