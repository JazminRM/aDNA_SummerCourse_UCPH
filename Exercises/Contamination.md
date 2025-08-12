## Estimating modern human contamination in ancient human data

### Outline

* Requesting an interactive node
* Estimate contamination based on the mtDNA genome
    + Create a consensus FASTA sequence for your ancient genome using ANGSD <sup>1</sup>
    + Remap the mtDNA reads of your ancient genome to its own consensus using paleomix <sup>2</sup> 
    + Estimate contamination using contamMix <sup>3</sup>

### Interactive node

First, get an interactive node. 

```{bash, eval = FALSE}
# log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# request one CPU using salloc like this:
salloc --partition=cpuqueue --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000 --qos=teaching --reservation=aDNA_PHD_course --account=teaching

# login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```

Now that you are in a node of the server, we can create a directoy for the contamination exercises:
```{bash, eval = FALSE}
# remember to change ku_username for your username:
username="your_ku_username"
mkdir -p /projects/course_1/people/${username}/Contamination/

cd /projects/course_1/people/${username}/Contamination/
```

### Estimate mtDNA based contamination using contamMix

#### Get a consensus sequence for your ancient genome

The first step to estimate contamination based on the mtDNA genome will be to estimate a majority count consensus sequence for the mtDNA genome of our ancient individual. We will keep working on the same ~10,000 year-old individual from Lagoa Santa cave in Brazil <sup>4</sup> that you mapped earlier, but for this exercise we have a new BAM file that consist of only the mtDNA reads. 

To estimate the majority count consensus sequence we will use `ANGSD -doFasta`.
```{bash, eval=FALSE}
#first load angsd
module load angsd/0.921

#then run angsd like this:
angsd -minQ 20 -minMapQ 30 -doCounts 1 -setMinDepth 3 -doFasta 2 -basesPerLine 60 -i /projects/course_1/people/clx746/Data/Sumidouro6_mt.bam -out S6_mt
```

Take a look at the parameters:
```
-minQ 20           minimum base quality
-minMapQ 30        minimum mapping quality
-doCounts 1        this is required in order to build a consensus sequence
-setMinDepth 3     minimum depth of coverage for a site to be included
-doFasta 2         2 indicates you want a majority-count consensus sequence, if you want a secuence with a random base you can use 1 instead
-basesPerLine 60   Bases per line that you want printed in the output
-i                 input BAM file
-out               basename for the output
```

Once it is done running you'll have a fasta file with a majority count consensus sequence for the mtDNA genome: 
```{bash, eval=FALSE}
ls *fa.gz
```
```
S6_mt.fa.gz
```

#### Re-map the reads to the consensus sequence

The second step is to re-map the data to its own consensus sequence. To make the mapping process easier, we will use paleomix <sup>2</sup>. 

Start by preparing the consensus sequence (FASTA file) we just generated with `angsd` for mapping. 
```{bash, eval=FALSE}
# load the module we will need
module load bwa/0.7.15
module load samtools/1.21
module unload openjdk/20.0.0 
module unload jdk/1.8.0_291 
module unload picard/3.4.0
module load picard/3.4.0

# unzip and change the name of the FASTA file
zcat S6_mt.fa.gz > S6_mt.fasta

# create a samtools index
samtools faidx S6_mt.fasta

# create a bwa index
bwa index S6_mt.fasta

# create the picard-tools index
picard CreateSequenceDictionary R=S6_mt.fasta O=S6_mt.dict
```

Then we want to extract the reads from our BAM file with the mtDNA data and turn them into a FASTQ that we can use for mapping:
```{bash, eval=FALSE}
# use picard SamToFastq to create a FASTQ file with the ancient mtDNA reads:
picard SamToFastq VALIDATION_STRINGENCY=LENIENT I=/projects/course_1/people/clx746/Data/Sumidouro6_mt.bam F=S6_mt.fq
```

Next we will re-map the data to its own consensus sequence using `paleomix`:
```{bash, eval=FALSE}
# load paleomix module
module load paleomix/1.3.9

# create a YAML file by copying the template from my folder and changing my username for yours (remember to change your_username for your own user name): 
cat /projects/course_1/people/clx746/Data/mapS6MT.yaml |sed 's/clx746/your_username/' > mapS6MT.yaml 

# run paleomix 
paleomix bam run --jre-option=-Xmx2g --jar-root /projects/course_1/people/clx746/Data/picard mapS6MT.yaml
```
It will take a few minutes to run, so just wait...


#### Estimate contamination using contamMix

Now we have almost everything we need to run contamMix. We will now add our mtDNA consensus sequence to an alingment consisting of 311 mtDNA human genomes using `mafft`.

```{bash, eval=FALSE}
# load the modules you need:
module load gcc
module load R/3.6.1 
module load mafft/7.515

# run mafft:
mafft --auto --add S6_mt.fasta /projects/course_1/people/clx746/Data/311humans.fa > S6_alignment.fa
```

Finally, we will run contamMix using the alignment and our re-mapped data:
```{bash, eval=FALSE}
Rscript /projects/course_1/people/clx746/Data/estimate.R --samFn S6.MT.bam --malnFn S6_alignment.fa --figure S6_contamix.pdf --tabOutput --nChains 3 --save S6_chains > s6_contamination.txt
```

It will take a while to run, so once you get it running take a break. 

...

Once it is done running, you can take a look at the contamination estimate: 

```{bash, eval=FALSE}
cat s6_contamination.txt
```
```
#<inferred-error-rate> <MAP-authentic> <2.5% authentic> <97.5% authentic> <gelman diagnostic> <gelman diag upper bound>
0.0143    0.9920751       0.9705297       0.9990772       1.008536        1.030358
```
We are interested in the second column which is the posterior probability that the mtDNA genome is authentic (MAP-authentic), our contamination estimate is 1-(MAP-authentic). 

<span stype="color:red"> **Question** </span> Would you say our ancient genome is contaminated? How do you think this could affect our analyses?


### References

1. Korneliussen, T.S *et al.* 2014. **ANGSD: Analysis of Next Generation Sequencing Data.** BMC Bioinformatics 15, 356.

2. Schubert M, *et al.* (2014) **Characterization of ancient and modern genomes by SNP detection and phylogenomic and metagenomic analysis using PALEOMIX**. Nat Protoc. May;9(5):1056-82.

3. Fu Q, *et al.* (2013). **A revised timescale for human evolution based on ancient mitochondrial genomes**. Curr Biol, 23:553-9.

4. Moreno-Mayar J.V. *et al*. (2018). **Early human dispersals within the Americas.** Science. 362, 1128.












