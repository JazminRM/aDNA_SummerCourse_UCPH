## Estimating heterozygocity with ancient sequencing data

### Outline of the exercises: 

* Estimating heterozygosity using ```ANGSD``` and GL
  + Estimate heterozygosity for all wolf and dog samples that we have been working with
  + Check if we have any bias in the heterozygosity estimates due to differences in depth of coverage


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
mkdir -p /projects/course_1/people/${username}/Diversity/

cd /projects/course_1/people/${username}/Diversity/
```

### Heterozygosity

We are going to use `ANGSD` and `realSFS` to estimate heterozygosity for the wolves and dogs we have been studying. 

Create a list of samples for which we want to estimate heterocigociy (this time we will all use the three ancient genomes): 

```{bash, eval = FALSE}
cat /projects/course_1/people/clx746/Bams/BamList1.txt > BamList.txt
ls /projects/course_1/people/clx746/Bams/sample2.bam >> BamList.txt
ls /projects/course_1/people/clx746/Bams/sample3.bam >> BamList.txt
```

The first step is to use ```NGSadmix``` to estimate a SAF file: 

```{bash, eval = FALSE}
# load angsd
module load angsd/0.921 

# we will do this for each sample independently (this will take a few minutes, so take a look at the parameters below while you wait):
cat BamList.txt | while read bam
do
samplename=$(basename $bam |cut -f 1 -d ".")
angsd -i $bam -ref /projects/course_1/people/clx746/Data/L.Dalen_14_wolf.scf.fasta -anc /projects/course_1/people/clx746/Data/GJ.fasta -out $samplename -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -noTrans 1 -C 50 -minMapQ 30 -minQ 20 -setMinDepth 3 -setMaxDepth 60 -doCounts 1 -GL 2 -doSaf 1
done
```

What parameters are we using:
```
-ref             Our reference genome that we used for mapping 
-anc             an outgroup or ancestral genome that can be used for polirizing the SFS
-uniqueOnly      remove reads mapping to more than one location
-C 50            reestimate the BAC, it is recommended by ANGSD 
-minMapQ         minimum mapping quality
-minQ            minimum base quality
-rmTrans 1       skip transition sites to reduce aDNA-derived error
-setMinDepth 3   minimun depth of coverage per site
-setMaxDepth 60  exclude sites with too high depth of coverage
-doCounts 1      necessary for doSaf
-doSaf 1         estimate GL 
```

Once it is done running, you will have 4 files for each of the samples:
```{bash, eval = FALSE}
ls Dog_Dingo*
```
```
Dog_Dingo.saf.gz
Dog_Dingo.saf.pos.gz
Dog_Dingo.saf.idx
Dog_Dingo.arg
```
The **.saf.gz** file corresponds to the GL (on this case they are written in a binary file).
The **.pos.gz** file contains the positions for which we were able to estimate GL for this particular sample. 
The **.saf.idx** and **.arg** files are the index and arguments from the run. 

Now we can use the GL and ```realSFS``` to estimate the SFS for each of the samples independently: 
```{bash, eval = FALSE}
cat BamList.txt | while read bam
do
samplename=$(basename $bam |cut -f 1 -d ".")
realSFS ${samplename}.saf.idx > ${samplename}.ml
done
```

Now let's look at the results for one of the samples:  
```{bash, eval = FALSE}
cat Dog_Dingo.ml
```
```
5782875.433847 7508.215427 22721.350728 
```
<span style="color: purple;"> **Q:** </span> Do you know what is represented in each column? How would you use that to estimate the heterozygosity?

<details>
<summary> <b>Show answer</b> </summary>

Each column represents the maximum likelihood estimate of the number of sites with 1 allele (ancestral), 2 alleles and 1 allele (derived). So we can use the second column as a proxy for the number of heterozygous sites in the genome and estimate the heterozygosity like this:
```
het= 7508.215427/(5782875.433847+7508.215427+22721.350728)
0.001291602
```
</details>

Now let's summarise and plot the results:
```{bash, eval = FALSE}
R
# get file names
f<-dir(pattern="ml$")
samplenames<-gsub(".ml", "", f)

# read files and estimate het
het<-NULL
for(i in 1:length(f)){
  a<-read.table(f[i], as.is=T)
  het<-c(het, as.numeric(a[2]/sum(a)))
}

# assign categories
info<-read.table("/projects/course_1/people/clx746/Data/wolves_rand_tv_info.txt", as.is=T, sep="\t")

# Find the category for each of the samples in the list
category<-NULL
for(i in 1:length(samplenames)){
 if(sum(info[,1]==samplenames[i])>0){
  category<-c(category, info[info[,1]==samplenames[i],2])
 }else{
  category<-c(category, samplenames[i]) 
 }
}

d<-data.frame(sample=samplenames, heterozygosity=het, category=category)
write.table(d, quote=F, sep="\t", col.names=T, row.names=F, file="Het.txt")

library(ggplot2)

colvalues <-c("#F3E96B", "coral4", "orange", "lightsalmon1", "#F05837","#6465A5", "darkmagenta", "darkred", "black", "black", "black")
names(colvalues)<-c("Grey wolf America", "Grey wolf Asia Highland", "Grey wolf Middle East", "Grey wolf Europe", "Grey wolf Asia", "Dog", "Dog Arctic", "Ancient grey wolf", "sample1", "sample2", "sample3")

p<-ggplot(d, aes(x=category, y=heterozygosity))+
geom_jitter(aes(colour=category), alpha=.5, size=2, shape=21, fill="black")+
geom_boxplot(aes(colour=category), outlier.shape=NA, fill=NA, notch=F, lwd=.7)+
scale_colour_manual(values=colvalues)+
theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=10), legend.key=element_rect(fill="white", colour="white"), legend.position = "right", axis.title=element_text(size=14), legend.text=element_text(size=10))

# now we the het
pdf("Het_GL.pdf", width=8, height=5)
print(p)
dev.off()

q("no")
```

Download and take a look at the results (use **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

What do you think about the results? 


#### Effect of depth of coverage

Now let's see if we can see an effect of the depth of coverage on our results. We will first estimate the depth of coverage for all of our samples and the plot it against the heterozygosity 

Estimate depth of coverage in the same way we did for our ancient sample a couple of days ago: 
```{bash, eval = FALSE}
cat BamList.txt | while read bam
do
samplename=$(basename $bam |cut -f 1 -d ".")
samtools depth -q 30 $bam  |grep -vP '\t0$' |cut -f  3 |sort |uniq -c > ${samplename}.doc.txt
done
```

Now let's summarise and plot the results:

```{bash, eval = FALSE}
R
# let's install ggrepel first (it might ask you to confirm that you want a local instalation):
install.package("ggrepel")

# get file names
f<-dir(pattern="doc.txt$")
samplenames<-gsub(".doc.txt", "", f)

# read files and estimate average doc
a_doc<-NULL
for(i in 1:length(f)){
  a<-read.table(f[i], as.is=T)
  a_doc<-c(a_doc, as.numeric(sum(a[,1]*a[,2])/sum(a[,1])))
}
d<-data.frame(sample=samplenames, a_doc=a_doc)
write.table(d, quote=F, sep="\t", col.names=T, row.names=F, file="Doc.txt")
het<-read.table("Het.txt", as.is=T, h=T, sep="\t")

if(sum(het$sample==samplenames)==length(samplenames)){
  d<-data.frame(sample=d$sample, a_doc=d$a_doc, heterozygosity=het$heterozygosity, category=het$category)
}else{
  message("Error")
}

library(ggplot2)
library(ggrepel)

colvalues <-c("#F3E96B", "coral4", "orange", "lightsalmon1", "#F05837","#6465A5", "darkmagenta", "darkred", "black", "black", "black")
names(colvalues)<-c("Grey wolf America", "Grey wolf Asia Highland", "Grey wolf Middle East", "Grey wolf Europe", "Grey wolf Asia", "Dog", "Dog Arctic", "Ancient grey wolf", "sample1", "sample2", "sample3")

p<-ggplot(d, aes(x=a_doc, y=heterozygosity))+
geom_point(aes(colour=category), shape=16, size=2.5, alpha=0.8, show.legend = T)+
scale_colour_manual(values=colvalues)+
geom_text_repel(data=d[d$sample %in% c("sample1", "sample3", "sample2"),], aes(x=a_doc, y=heterozygosity, label=sample), segment.size=0.2, segment.color="black", size=2.5, box.padding=unit(0.025, 'npc'))+
theme(panel.background = element_rect(fill='white', colour='black'), axis.text.x = element_text(angle = 45, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'), strip.background = element_rect(size=.2, colour="black", fill ="#FFA5004D"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=10), legend.key=element_rect(fill="white", colour="white"), legend.position = "right", axis.title=element_text(size=14), legend.text=element_text(size=10))

# now we plot het vs doc
pdf("Het_vs_doc.pdf", width=8, height=5)
print(p)
dev.off()

q("no")

```

Download and take a look at the results (use **WinSCP** (for Windows users) or **scp** (for Mac or Linux users). 

Take a look at the results, do they look as expected? Do you think we have a bias due to the coverage?
