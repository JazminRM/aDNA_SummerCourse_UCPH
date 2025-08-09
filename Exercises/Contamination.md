## Estimating modern human contamination in ancient human data

### Outline


### Interactive node

```{bash, eval=FALSE}

cd /projects/course_1/people/clx746/Contamination
```

### Estimate mtDNA based contamination using contamMix

Estimate a majority count consensus sequence using `ANGSD`

```{bash, eval=FALSE}
module load angsd/0.921

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

Prepare reference genome

```{bash, eval=FALSE}
module load bwa/0.7.15
module load samtools/1.21
module unload openjdk/20.0.0 
module unload jdk/1.8.0_291 
module unload picard/3.4.0
module load picard/3.4.0

gunzip S6_mt.fa.gz

mv  S6_mt.fa  S6_mt.fasta
samtools faidx S6_mt.fasta

bwa index S6_mt.fasta

picard CreateSequenceDictionary R=S6_mt.fasta O=S6_mt.dict

```

```{bash, eval=FALSE}

picard SamToFastq VALIDATION_STRINGENCY=LENIENT I=/projects/course_1/people/clx746/Data/Sumidouro6_mt.bam F=S6_mt.fq

module load paleomix/1.3.9

cat /projects/course_1/people/clx746/Data/mapS6MT.yaml |sed 's/cls746/your_username/' > mapS6MT.yaml 

# run paleomix 
paleomix bam run --jre-option=-Xmx2g --jar-root /projects/symbean/people/clx746/Scripts/picard/ mapS6MT.yaml



module load gcc
module load R/3.6.1 
module load mafft/7.515

mafft --auto --add S6_mt.fasta /projects/course_1/people/clx746/Data/311humans.fa > S6_alignment.fa

Rscript /projects/course_1/people/clx746/Data/estimate.R --samFn S6.MT.bam --malnFn S6_alignment.fa --figure S6_contamix.pdf --tabOutput --nChains 3 --save S6_chains


#  <inferred-error-rate> <MAP-authentic> <2.5% authentic> <97.5% authentic> <gelman diagnostic> <gelman diag upper bound>


module load angsd/0.921

angsd -i /projects/course_1/people/clx746/Data/Sumidouro6_final.bam -r X:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 30 -minQ 20 -out s6

contamination mapFile="/projects/course_1/people/clx746/Data/map100.chrX.gz" hapFile="/projects/course_1/people/clx746/Data/hapMapCeuXlift.map.gz" -a s6.icnts.gz -p 1

contamination -h /projects/course_1/people/clx746/Data/map100.chrX.gz -m /projects/course_1/people/clx746/Data/hapMapCeuXlift.map.gz -a s6.icnts.gz -p 1

```












