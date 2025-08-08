## Basic steps for mapping ancient DNA sequencing data

### Outline: 

* Requesting an interactive node
* FASTQ files
* Data quality check with $fastqc$
* Adapter trimming
* Step-by-step mapping 
    + Preparing the reference genome
    + Mapping 
    + BAM files
    + Post-mapping filters
    + Identifying and removing PCR-duplicate
* Mapping statistics


### Interactive node

We will use an interactive node in mjolnir to run the exercises throughout the course, so you always start by opening an interactive session in the server. Follow the steps below.

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# first request one CPU using salloc like this:
salloc --partition=cpuqueue --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000 --account=teaching --reservation=aDNA_PHD_course

# once the job has been allocated, you can login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```

Now you are in a node of the server, where you can run the exercises. Each of you has a directory assigned where you can store and process data: 

```{bash, eval = FALSE}
# remember to change ku_username for your username:
cd /projects/course_1/people/ku_username/
```

**NOTE 1**: When you are done working, remember to logout and close the interactive session, you can do that by typing exit twice. It will look something like this: 

```
clx746@mjolnircomp04fl fst]$ exit
exit
clx746@mjolnirgate01fl fst]$ exit
exit
salloc: Relinquishing job allocation 31797328
```
**NOTE 2**: If you loose internet connection or disconnect from the server (for example during the coffee break or lunch), you'll need to start a new interactive session, or login again to the one you had. 

-----------------------------------


### Setting up your working directory and paths to data

Define some paths and file names. We will processes FASTQ files that correspond to sequencing data from three sequencing libraries from ~10,000 year-old individual from Lagoa Santa cave in Brazil <sup>0</sup>.

Create a directory for today's exercises:

```{bash, eval = FALSE}
# remember to change ku_username for your username:
mkdir -p /projects/course_1/people/ku_username/mapping/

# and go to your working directory:
cd /projects/course_1/people/ku_username/mapping/

FASTQ1="/projects/course_1/people/clx746/Data/PW13_E2_L1.fastq.gz"
FASTQ2="/projects/course_1/people/clx746/Data/PW13_E1_L1.fastq.gz"
FASTQ3="/projects/course_1/people/clx746/Data/PW13_E2_L2.fastq.gz"
```

#### FASTQ files refresher

Quick refreser on what are the FASTQ files and how they look.

```{bash, eval = FALSE}
# You can take a look at one of the FASTQ files like this:
zcat $FASTQ1 |head -n 8
```
```
@HISEQ:151:C8L9GANXX:1:1101:1419:2130 1:N:0:GTGTAT
TTCGGCGGCGGTGTCGCGGATGTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGTATATCTCGTATGCCGTCTT
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@HISEQ:151:C8L9GANXX:1:1101:1486:2132 1:N:0:GTGTAT
TTTGCATCTAACTGAAACTTTCTTCACTTTGAATACAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGTATATCT
+
BBCC0@CCFGCGFGGGGGGGGFGGGGGGGGGFGGGGGGGGG@CFEBGBGGGGGGEGGGGGGGGGEGGGGGGGGGGGGFF@F
```

<span stype="color:red"> **Question** </span> Do you remember what each line represents? How many reads are there in each FASTQ file? 

<details>
<summary> <b>Show answer</b> </summary>
Each read in the FASTQ file is represented by four lines: 
    1. read ID, 
    2. sequence, 
    3. a + sign and 
    4. the qualities. 
So you can count the number of read by dividing the number of lines in the file by 4: 

```{bash, eval = FALSE}
zcat $FASTQ1 |wc -l 
```
```
16000000
```
Total number of reads: 16000000/4=4000000
</details>

-----------------------------------



### Data quality check with fastqc

A common step before trimming and mapping sequencing data is to check its quality using `fastqc`<sup>1</sup>, which estimates some basic quality statistics. 

We will run $fastqc$ for each of our FASTQ files
```{bash, eval = FALSE}

# Load the necessary programs/modules:
module load perl/5.24.4 
module load openjdk
module load fastqc/0.12.1 

# Create a directory to save the results
mkdir RawQual
# run fastqc (run one line at a time):
fastqc --outdir ./RawQual -f fastq $FASTQ1
fastqc --outdir ./RawQual -f fastq $FASTQ2
fastqc --outdir ./RawQual -f fastq $FASTQ3
```

`fastqc` will create an HTML file for each of the FASTQs. Download the files to your local computer so you can take a look. 

You can use **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

Example of the **scp** command:

```{bash, eval = FALSE}
scp clx746@mjolnirgate.unicph.domain:/projects/course_1/people/clx746/Mapping/RawQual/PW13_E1_L1_fastqc.html .
scp clx746@mjolnirgate.unicph.domain:/projects/course_1/people/clx746/Mapping/RawQual/PW13_E2_L1_fastqc.html .
scp clx746@mjolnirgate.unicph.domain:/projects/course_1/people/clx746/Mapping/RawQual/PW13_E2_L2_fastqc.html .
```
(remember to change the directory name to your own)

-----------------------------------


### Adapter trimming

Before mapping the sequencing read to the reference genome, we will remove the adapter sequences and do some quality trimming using `AdapterRemoval`<sup>2</sup> 

```{bash, eval = FALSE}
# load AdapterRemoval module:
module load adapterremoval/2.3.3
 
# Assign a name for the output:
bn1=`basename $FASTQ1 |sed 's/.fastq.gz//'`
# Run adapterRemoval (it might take a few minutes, wait until it is done before running the next one): 
AdapterRemoval --file1 $FASTQ1 --basename $bn1 --trimns --trimqualities --minquality 2 --minlength 30  --mm 3 --gzip

# Run adapterRemoval for the other two files:
bn2=`basename $FASTQ2 |sed 's/.fastq.gz//'`
AdapterRemoval --file1 $FASTQ2 --basename $bn2 --trimns --trimqualities --minquality 2 --minlength 30  --mm 3 --gzip

bn3=`basename $FASTQ3 |sed 's/.fastq.gz//'`
AdapterRemoval --file1 $FASTQ3 --basename $bn3 --trimns --trimqualities --minquality 2 --minlength 30  --mm 3 --gzip
```

**NOTE:** We are assuming the sequences have the standard Illumina adapter sequence (default in AdapaterRemoval). If your sequences have a different adapter, make sure to provide it. 

Suggested filters:
```
--trimns            trim N's from the 5' and 3' ends
--trimqualities     trim low qualiy bases  from the 5' and 3' ends (min base quality defined by --minquality)
--minquality 2      default min quality for trimming low quality stretches from the ends
--minlength 30      minimum read length after removing adapters
--mm 3              max. error rate for the aligning the adapter sequence with the reads (1/3: default)
```
If you have paired-end data you can also use the ```--collapse``` option to collapse overlapping read pairs. 

Important output files to look at:

```{bash, eval = FALSE}
ls -l $bn1.*
```
```
PW13_E2_L1.discarded.gz
PW13_E2_L1.settings
PW13_E2_L1.truncated.gz
```
The `PW13_E2_L1.truncated.gz` and `PW13_E2_L1.discarded.gz` files contain the trimmed and discarded reads, respectively. The `PW13_E2_L1.settings` file has useful information about the number of sequences where the adapter was found, number of remaining reads, etc.


```{bash, eval = FALSE}
grep -A 10 'Trimming statistics' $bn1.settings
```
```
[Trimming statistics]
Total number of reads: 4000000
Number of unaligned reads: 134233
Number of well aligned reads: 3865767
Number of discarded mate 1 reads: 1598722
Number of singleton mate 1 reads: 2401278
Number of reads with adapters[1]: 3865767
Number of retained reads: 2401278
Number of retained nucleotides: 99392433
Average length of retained reads: 41.3915
```
Notice how the adapter sequence was found in the majority of the reads (3865767/4000000 reads), with ancient DNA data this is good indication that the data was trimmed correctly (and you had the correct adapter sequence).   

<details>
<summary> <b>Optional</b> </summary>
    
Now that we removed the adapter, we could re-run `fastqc` to check how the sequences look like:

```{bash, eval = FALSE}
# Create a new directory:
mkdir TrimQual
# run fastqc in the trimmed data:
fastqc --outdir TrimQual -f fastq ${bn1}.truncated.gz
fastqc --outdir TrimQual -f fastq ${bn2}.truncated.gz
fastqc --outdir TrimQual -f fastq ${bn3}.truncated.gz
```

Download the new HTML files and compare the results from before and after trimming the adapters. 

</details>

-----------------------------------


### Step-by-step mapping

#### Preparing the reference genome

Things to consider about the reference genome. Most likely you'll be using the same reference genome for every analysis, so make sure you have a good reference genome. Some things to consider are: chromosome names, no weird characters in the file, are you using the correct version? is it masked?, is it complete? are you including the mtDNA/cpDNA sequences too? 

For this exercise, we will be mapping our data to the Human reference genome (build37.1) <sup>3</sup>:

```{bash, eval = FALSE}
REF_GENOME="/projects/course_1/people/clx746/Data/Genome/hs.build37.1.fa"
```

First time you use a reference genome, you need to create three indexes: a ```samtools index```, a ```picard``` dictionary and a ```bwa index```. In this case we already created them since each can take some time to run. 

But you can see here how to create them here:

```
# bwa index:
bwa index hs.build37.1.fa

# samtools index:
samtools faidx hs.build37.1.fa

# picard dictionary/index:
picard CreateSequenceDictionary R=hs.build37.1.fa O=hs.build37.1.dict
```

#### Mapping

Here you can follow a step-by-step mapping of our ancient. As we mentioned in class, many times you'll use an automatic workflow that will do each of these steps for you, however, it is important to understand the main processing and filtering steps, so we are running each one manually here.

First step, map your data to the reference genome using ```bwa aln```<sup>4</sup>: 

```{bash, eval = FALSE}
# load bwa module
module load bwa/0.7.15
module load samtools/1.21

# map each of the files separately (each will take a few minutes, wait until it is done before running the next one)
bwa aln -l 1000 $REF_GENOME ${bn1}.truncated.gz >  ${bn1}.sai

bwa aln -l 1000 $REF_GENOME ${bn2}.truncated.gz >  ${bn2}.sai

bwa aln -l 1000 $REF_GENOME ${bn3}.truncated.gz >  ${bn3}.sai

```

One important parameter here is the `-l`. The `-l` specifies the length of the sequence that will be used as seed during the mapping (first 32bp by default), we want to use the complete read as seed, given that we expect most of the aDNA damage will located at the start of the read and only using the first bps will decrease the chances of the read mapping.

Second step, run ```bwa samse``` (or ```bwa samppe``` for paired-end data). 

Something important in this step is to assign a read group with the information about the sample name, library ID, read group ID and sequencing platform. This information is used by other programs to e.g. identify reads from the same library. 

You can define the read group like this:
```
@RG\tID:PW13_E2_L1\tLB:PW13_E2_L1\tPL:ILLUMINA\tSM:S6

ID = read group ID
LB = library ID
PL = sequencing platform
SM = sample name
```

```{bash, eval = FALSE}
# in each of these lines, we are first using bwa samse to create a SAM file with the read group information, and they we are passing this information to samtools to conver the SAM file into a BAM file: 
bwa samse -r '@RG\tID:PW13_E2_L1\tLB:PW13_E2_L1\tPL:ILLUMINA\tSM:S6' $REF_GENOME ${bn1}.sai ${bn1}.truncated.gz | samtools view -Sbho ${bn1}.bam

bwa samse -r '@RG\tID:PW13_E1_L1\tLB:PW13_E1_L1\tPL:ILLUMINA\tSM:S6' $REF_GENOME ${bn2}.sai ${bn2}.truncated.gz | samtools view -Sbho ${bn2}.bam

bwa samse -r '@RG\tID:PW13_E2_L2\tLB:PW13_E2_L2\tPL:ILLUMINA\tSM:S6' $REF_GENOME ${bn3}.sai ${bn3}.truncated.gz | samtools view -Sbho ${bn3}.bam

```

#### BAM files refresher

Let's look at how our BAM files. Given they are compressed files, we can only read them using `samtools`<sup>5</sup>. 

For example:

```{bash, eval = FALSE}
samtools view ${bn1}.bam |head -n 5
```
```
HISEQ:151:C8L9GANXX:1:1101:1486:2132    0       22      50581434        37      37M     *       0       0       TTTGCATCTAACTGAAACTTTCTTCACTTTGAATACA   BBCC0@CCFGCGFGGGGGGGGFGGGGGGGGGFGGGGG   RG:Z:FRC_1   XT:A:U   NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:37
HISEQ:151:C8L9GANXX:1:1101:1434:2178    16      8       111192000       37      41M     *       0       0       TGAGGATTCCCTTTCTTCCACATACTTGCCAGCATTTGTTA       >GF1GGGGGGGGGGFEGGGGGF1GGE1<@GGGGFFB=?3BA    RG:Z:FRC_1       XT:A:U  NM:i:1  X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  MD:Z:0C40
HISEQ:151:C8L9GANXX:1:1101:1392:2236    4       *       0       0       *       *       0       0       GAACTCCAGTCACCTATCAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAATCATTTTTCTG      3330AG;1EGGG1DGGGGF1=;==EFCGGCFGGGGB:@F>FFE:BFDBG///<//0?0=0:=:B@C    RG:Z:FRC_1
HISEQ:151:C8L9GANXX:1:1101:1539:2153    16      8       81727741        37      50M     *       0       0       ATAATGCTTAACTGTTGAATTATTTTAAAGGCTGGACAATAAGTAATAAT      GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGB@BB@    RG:Z:FRC_1      XT:A:U  NM:i:2  X0:i:1  X1:i:0  XM:i:2  XO:i:0  XG:i:0  MD:Z:1C25G22
HISEQ:151:C8L9GANXX:1:1101:1683:2171    4       *       0       0       *       *       0       0       CGTACAAGGGCAAGGGAATCCGGTACGCGGGCGA      ABBBACDFGGCDGGGGGGG==<EGB==EGGGCEG      RG:Z:FRC_1
```

Or if you want to also look at the header:
```{bash, eval = FALSE}
samtools view -H ${bn1}.bam 
```
```
@SQ     SN:1    LN:249250621
@SQ     SN:2    LN:243199373
@SQ     SN:3    LN:198022430
...
@SQ     SN:21   LN:48129895
@SQ     SN:22   LN:51304566
@SQ     SN:X    LN:155270560
@SQ     SN:Y    LN:59373566
@SQ     SN:MT   LN:16569
@RG     ID:FRC_1        LB:FRC  PL:ILLUMINA     SM:Tumat
@PG     ID:bwa  PN:bwa  VN:0.7.15-r1140 CL:bwa samse -r @RG\tID:FRC_1\tLB:FRC\tPL:ILLUMINA\tSM:Tumat /projects/course_1/people/clx746/Data/Genome/hs.build37.1.fa PW13_E2_L1.sai PW13_E2_L1.truncated.gz
```

We can also use $samtools$ to get some basic statistics about the number of reads mapped and unmapped:

```{bash, eval = FALSE}
samtools flagstat ${bn1}.bam
```

```
993038 + 0 in total (QC-passed reads + QC-failed reads)
993038 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
31739 + 0 mapped (3.20% : N/A)
31739 + 0 primary mapped (3.20% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
<span style="color: purple;"> **Question:** </span>  How many reads mapped to the reference genome? Is is the same percentage for the other two BAM files?  

<details>
<summary> <b>Show answer</b> </summary>
Check the stats on the other two files and compare the number mapped reads. 
```{bash, eval = FALSE}
samtools flagstat ${bn2}.bam
```

```{bash, eval = FALSE}
samtools flagstat ${bn3}.bam
```
</details>

**NOTE:** If you want to know more about what is coded in each of the BAM fields you can take a look at the the SAM/BAM format specification: https://samtools.github.io/hts-specs/SAMv1.pdf

##### Post-mapping filters

Next step is to remove the unmapped reads and filter by mapping quality using ```samtools```:

```{bash, eval = FALSE}
samtools view -bh -o ${bn1}.mapped.bam -F4 -q 30 ${bn1}.bam

samtools view -bh -o ${bn2}.mapped.bam -F4 -q 30 ${bn2}.bam

samtools view -bh -o ${bn3}.mapped.bam -F4 -q 30 ${bn3}.bam
```

Important parameters:
```
-F4        filter reads unmapped (with flag 4)
-q 30      filter reads with mapping quality lower than 30. 

```
You can check what different flags indicate here: https://broadinstitute.github.io/picard/explain-flags.html

Next, we need to sort the reads within the BAM files (```picard```, which we will use next, needs sorted BAM files). 

```{bash, eval = FALSE}
samtools sort -o ${bn1}.mapped.sorted.bam ${bn1}.mapped.bam 

samtools sort -o ${bn2}.mapped.sorted.bam ${bn2}.mapped.bam

samtools sort -o ${bn3}.mapped.sorted.bam ${bn3}.mapped.bam
```

##### Merge files by sample

Now that we have our mapped reads sorted we can merge all the BAMS that correspond to the same sample into a single BAM using ```picard MergeSamFiles```<sup>6</sup>: 

```{bash, eval = FALSE}
# Create a sample name for the merged BAM file: 
samplen="AncientWolf"
java -jar /home/ec2-user/Software/picard/picard.jar MergeSamFiles VALIDATION_STRINGENCY=LENIENT I=${bn1}.mapped.sorted.bam I=${bn2}.mapped.sorted.bam I=${bn3}.mapped.sorted.bam  O=${samplen}.merged.bam
```

Since we added a different read group ID to our BAM files before merging, we can always use this information to split the reads if needed.

##### PCR duplicates removal

The next step is to remove the PCR duplicates using $picard MarkDuplicates$:

```{bash, eval = FALSE}
java -jar /home/ec2-user/Software/picard/picard.jar MarkDuplicates I=${samplen}.merged.bam O=${samplen}.merged.dedup.bam METRICS_FILE=${samplen}.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
```

$picard MarkDuplicates$ will create a metrics file with the number of reads tha were PCR duplicates:

```{bash, eval = FALSE}
cat ${samplen}.metrics.txt |grep -v "#" |column -t
```
```
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
FRC	60113	0	0	0	430	0	0	0.007153	
Liver	18689	0	0	0	2265	0	0	0.121194	
```
<span style="color: purple;"> **Q:** </span> Do the libraries had the same proportion of pcr-duplicates? Which one is better? 

##### Index our BAM

Finaly, we will index our BAM file:
```{bash, eval = FALSE}
samtools index ${samplen}.merged.dedup.bam
```

Finally, we are ready to check how many reads remain after filtering: 

```{bash, eval = FALSE}
samtools view -c ${samplen}.merged.dedup.bam
```
```
76107
```
<span style="color: purple;"> **Q:** </span> How many useful reads remain? Try estimating the % of endogenous DNA/reads for each of the libraries. Hint: check which filtering options you can give to ```samtools``` by typing ```samtools view```

<button class="btn btn-primary" button style="background-color:purple; border-color:purple; color:white" data-toggle="collapse" data-target="#BlockName4"> Show/hide solution </button>  
<div id="BlockName4" class="collapse">  

You can use samtools to count the number of reads in each library
```{bash, eval = FALSE}
samtools view -l FRC -c  ${samplen}.merged.dedup.bam
```
```
59683
```
```{bash, eval = FALSE}
samtools view -l Liver -c  ${samplen}.merged.dedup.bam
```
```
16424
```
And you can get the number of initial reads from the *.settings* file from AdapterRemoval:
```{bash, eval = FALSE}
grep "Total" *settings
```
```
TOG_KCCS_FRC_1.settings:Total number of reads: 1000000
TOG_KCCS_FRC_2.settings:Total number of reads: 1549421
TOG_KCCS_L.settings:Total number of reads: 5827604
```
Use those numbers to estimate the endogenous content of each library. 
</div>

<p>&nbsp;</p>

---------------------------------

#### Authentication

Before start analysing the data, we would like to check whether the reads we mapped look authentic ancient DNA, so we will use $bamdamage$<sup>7</sup> to estimate substitution patterns in the reads for each of the two sequencing libraries.

##### MapDamage

Since we want to compare how the damage looks for each of the sequencing libraries, we will first split our BAM into two:
```{bash, eval = FALSE}
# separate the reads from the Liver library using samtools and index the BAM:
samtools view -r Liver -b -o  ${samplen}.merged.dedup.liver.bam  ${samplen}.merged.dedup.bam 
samtools index  ${samplen}.merged.dedup.liver.bam 

# Now separate the cartilage reads (since we have two different read groups, we need to create a file with the IDs: 
echo 'FRC_1
FRC_2
' > FRC_rg.txt
samtools view -R FRC_rg.txt -b -o  ${samplen}.merged.dedup.frc.bam ${samplen}.merged.dedup.bam
samtools index  ${samplen}.merged.dedup.frc.bam 

```

Now, let's run ```mapDamage```:

```{bash, eval = FALSE}
mapDamage -i  ${samplen}.merged.dedup.liver.bam -r $REF_GENOME -d ${samplen}_mapDamage_liver 
mapDamage -i  ${samplen}.merged.dedup.frc.bam -r $REF_GENOME -d ${samplen}_mapDamage_frc

```

```mapDamage``` creates several output files, including one PDF with the damage patterns and real length distribution: 
```{bash, eval = FALSE}
ls *_mapDamage_frc/*ploy.pdf
```
```
AncientWolf_mapDamage_frc/Fragmisincorporation_plot.pdf
AncientWolf_mapDamage_frc/Length_plot.pdf
```
Download both plots we just created to your local computer using **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

**NOTE:** PDFs for both libraries have the same names, so download each one at a time and change the names.  
Example of the **scp** command:

```{bash, eval = FALSE}
scp -i apgc-2021-key.pem.txt your_user_name@3.249.84.19:/home/your_user_name/mapping/*_mapDamage_liver/*.pdf .
scp -i apgc-2021-key.pem.txt your_user_name@3.249.84.19:/home/your_user_name/mapping/*_mapDamage_frc/*.pdf .
```
(remember to change the directory name to your own)

Do they look the same? Does it look authentic?
<p>&nbsp;</p>

---------------------------------

#### Mapping alternative using Paleomix

If you have many samples and FASTQ files that you want to map an alternative is to use ```paleomix```<sup>8</sup>, which run all of the steps above. 

First, create a directory so that we don't overwrite our previous file:

```{bash, eval = FALSE}
# First go into your own directory (remember to change for your own username):
username="your_own_username"
cd /home/$username/mapping/ 

# Then create a new directory:
directoryPaleomix="mapping_paleomix"
mkdir $directoryPaleomix
cd  $directoryPaleomix
```

In order to run $paleomix$, we need to fill in a template with the information about our FASTQs, the parameters we want to use and the samples information. 

Download the following pre-filled template:
```{bash, eval = FALSE}
wget https://sid.erda.dk/share_redirect/CaAx7aX67g/ancCanid.yaml
```

Take a look at the template and at the different options, do you recognise some of the parameters?

Once you have the template ready, $paleomix$ can be run like this:

```{bash, eval = FALSE}
paleomix bam run --jre-option=-Xmx2g --jar-root /home/ec2-user/Software/picard/ ancCanid.yaml
```

Once $paleomix$ is done running, you'll have your BAM file and two summary files:
```
AncCanid.summary
AncCanid.wolf.coverage
```

The first one (*AncCanid.summary*) is a summary of the number of reads mapped at different steps. The first lines contain information about the reference genome used:

```{bash, eval = FALSE}
 grep "^#" AncCanid.summary 
```

```
# Command:
#     /usr/local/bin/paleomix bam run --jre-option=-Xmx2g --jar-root /home/ec2-user/Software/picard/picard/build/libs/ /home/ec2-user/Day2/Paleomix/test.yaml
#
# Directory:
#     /home/ec2-user/Day2/Paleomix
#
# Genomes:
#     Name    Contigs    Size         Prefix
#     wolf    101        483573980    /home/ec2-user/TestData/RefGenome/Wolf.fasta
#
# Regions Of Interest:
#     Genome    ROI    Size    NFeatures    NIntervals    Path
```

Then, you have a summary table by sample followed by a summary table for each of the libraries.

This will get us the by-sample summary:
```{bash, eval = FALSE}
 grep -v "^#" AncCanid.summary |column -t | head -n 15
```
 
That looks like this:

 ```
Target    Sample Library  Measure                 Value                  #  Description                                                                                                     
AncCanid  *      *        lib_type                SE                     #  SE, PE, or * (for both)
AncCanid  *      *        seq_reads_se            8377025                #  Total number of single-ended reads
AncCanid  *      *        seq_trash_se            75028                  #  Total number of trashed reads
AncCanid  *      *        seq_trash_se_frac       0.008956401586482074   #  Fraction of SE reads trashed
AncCanid  *      *        seq_retained_reads      8301997                #  Total number of retained reads 

AncCanid  *      *        seq_retained_nts        581900328              #  Total number of NTs in retained reads                                          
AncCanid  *      *        seq_retained_length     70.09160904298086      #  Average number of NTs in retained reads                                          
AncCanid  *      *        hits_raw(wolf)          78928                  #  Total number of hits (prior to PCR duplicate filtering)            
AncCanid  *      *        hits_raw_frac(wolf)     0.009507110156749033   #  Total number of hits vs. total number of reads retained  
AncCanid  *      *        hits_clonality(wolf)    0.03343553618487738    #  Fraction of hits that were PCR duplicates                                     
AncCanid  *      *        hits_unique(wolf)       76289                  #  Total number of hits (excluding any PCR duplicates)                        
AncCanid  *      *        hits_unique_frac(wolf)  0.009189234831089436   #  Total number of unique hits vs. total number of reads retained
AncCanid  *      *        hits_coverage(wolf)     0.011778832682436718   #  Estimated coverage from unique hits
AncCanid  *      *        hits_length(wolf)       74.66262501802356      #  Average number of aligned bases per unique hit                                
 ```
 
 This will get us the by-library summary for the Liver library:
```{bash, eval = FALSE}
 grep -v "^#" AncCanid.summary |grep Liver |column -t
```

```
Target	  Sample	Library	Measure	              Value                # Description
AncCanid  AncCanid  Liver   lib_type              SE                   # SE, PE, or * (for both)
AncCanid  AncCanid  Liver   seq_reads_se          5827604              # Total number of single-ended reads
AncCanid  AncCanid  Liver   seq_trash_se          57676                # Total number of trashed reads
AncCanid  AncCanid  Liver   seq_trash_se_frac     0.009897034870591756 # Fraction of SE reads trashed
AncCanid  AncCanid  Liver   seq_retained_reads    5769928              # Total number of retained reads
AncCanid  AncCanid  Liver   seq_retained_nts      395275379            # Total number of NTs in retained reads
AncCanid  AncCanid  Liver   seq_retained_length   68.50611983373103    # Average number of NTs in retained reads 
AncCanid  AncCanid  Liver   hits_raw(wolf)        19369                # Total number of hits (prior to PCR duplicate filtering)
AncCanid  AncCanid  Liver   hits_raw_frac(wolf)   0.003356887642272139 # Total number of hits vs. total number of reads retained 
AncCanid  AncCanid  Liver   hits_clonality(wolf)  0.11905622386287362  # Fraction of hits that were PCR duplicates 
AncCanid  AncCanid  Liver   hits_unique(wolf)     17063                # Total number of hits (excluding any PCR duplicates) 
AncCanid  AncCanid  Liver   hits_unique_frac(wolf)0.002957229275651273 # Total number of unique hits vs. total number of reads retained
AncCanid  AncCanid  Liver   hits_coverage(wolf)   0.002547901770893462 # Estimated coverage from unique hits 
AncCanid  AncCanid  Liver   hits_length(wolf)     72.208814393717      # Average number of aligned bases per unique hit 
 ```
The second one (*AncCanid.wolf.coverage*) is a summary of the depth of coverage. The first lines correspond to a header describing what each column is:

```{bash, eval = FALSE}
grep "^#" AncCanid.wolf.coverage 
```

```
# Columns:
#   Contig:    Contig, chromosome, or feature for which a depth histogram was
#              created. Unnamed features are named after the chromosome or
#              contig on which they are located, with a star appended. For
#              example "chr1*".
#   Size:      The total size of the region. Multiple features with the same
#              name are combined into one row, with the size representing to
#              total of these. Note that overlapping bases are counted 2 (or
#              more) times.
#   Hits:      Sum of SE, PE_1, and PE_2 hits. Note that supplementary
#              alignments, duplicates, reads that failed QC, secondary
#              alignments, and unmapped reads are ignored.
#   SE, PE_*:  Number of Single Ended, and Pair Ended (mate 1 and 2) hits
#              overlapping the current contig or intervals. Note that a hit
#              may be counted multiple times if it overlaps multiple intervals
#   Collapsed: Number of hits for PE pair collapsed into a single read.
#   M, I, D:   Number of aligned (M), inserted (I) and deleted (D) bases
#              relative to references.
#   Coverage:  Average number of bases covering each position in the
#              contig(s)/intervals(s).
```

And the contents will look like this:
```{bash, eval = FALSE}
grep -v "^#" AncCanid.wolf.coverage  |column -t
```

```
Name      Sample    Library  Contig    Size       Hits   SE     PE_1  PE_2  Collapsed  M        I     D     Coverage
AncCanid  *         *        *         483573980  79054  79054  0     0     0          5888362  1058  1263  0.012176755250561661
AncCanid  *         *        <Genome>  483573980  79054  79054  0     0     0          5888362  1058  1263  0.012176755250561661
AncCanid  AncCanid  *        *         483573980  79054  79054  0     0     0          5888362  1058  1263  0.012176755250561661
AncCanid  AncCanid  *        <Genome>  483573980  79054  79054  0     0     0          5888362  1058  1263  0.012176755250561661
AncCanid  AncCanid  FRC_1    *         483573980  61991  61991  0     0     0          4656263  867   1095  0.009628853479668199
AncCanid  AncCanid  FRC_1    <Genome>  483573980  61991  61991  0     0     0          4656263  867   1095  0.009628853479668199
AncCanid  AncCanid  Liver    *         483573980  17063  17063  0     0     0          1232099  191   168   0.0025479017708934628
AncCanid  AncCanid  Liver    <Genome>  483573980  17063  17063  0     0     0          1232099  191   168   0.0025479017708934628
```
Do the number of final reads in the $paleomix$ BAM files match with our step-by-step BAM files?
<p>&nbsp;</p>

#### References
1. Andrews, S. (2010). **FastQC:  A Quality Control Tool for High Throughput Sequence Data.** Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
3. Schubert, M. *et al* (2016). **AdapterRemoval v2: rapid adapter trimming, identification, and read merging** BMC Research Notes. 
2. Gopalakrishnan, S *et al.* (2017) **The wolf reference genome sequence (*Canis lupus lupus*) and its implications for Canis spp. population genomics.** BMC Genomics 18, 495
4. Li, H. & Durbin, R. (2009) **Fast and accurate short read alignment with Burrows–Wheeler transform**. Bioinformatics 25, 1754–1760 .
5. Li, H. et al. (2009) **The Sequence Alignment/Map format and SAMtools**. Bioinformatics 25, 2078–2079.
6. picard-tools: http://picard.sourceforge.net
7. Malaspinas, A.-S. *et al.* (2014) **bammds: a tool for assessing the ancestry of low-depth whole-genome data using multidimensional scaling (MDS).** Bioinformatics 30, 2962–2964 
8. Schubert M, *et al.* (2014) **Characterization of ancient and modern genomes by SNP detection and phylogenomic and metagenomic analysis using PALEOMIX**. Nat Protoc. May;9(5):1056-82.


