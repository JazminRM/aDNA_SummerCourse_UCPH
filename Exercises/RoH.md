## Estimating Runs of Homozygocity with ancient genomes

We are going to estimate runs of homozygosity in a panel of 49 individuals (ancient and present-day) from different populations around the world. Africa: Yoruba; Oceania: Papuan; Europe: French, Icelandic; Asia: Punjabi, Han; Americas (ancient): Anzick1, Spirit Cave, Lovelock Cave, Taino; Americas (present-day): Mixe, Surui, Karitiana. 

For all individuals in the panel, we inferred their diploid genotypes using imputation (more on this tomorrow) and considered a subset of ~1.2M SNPs across the genome. 

### Interactive node 
We start by starting an interactive job:

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# request one CPU using salloc like this:
salloc --partition=cpuqueue --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000 --qos=teaching --reservation=aDNA_PHD_course --account=teaching

# once the job has been allocated, you can login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```

Create a directory for the exerise:

```{bash, eval = FALSE}
# change here for your username
username="your_username"

# create directory
directoryroh="/projects/course_1/people/${username}/RoHs/"
mkdir -p $directoryroh
# go into the directory
cd $directoryroh
```

### RoHs

Then, we load some modules:

```{bash, eval = FALSE}
module load gcc/13.2.0
module load openjdk/20.0.0
module load R/4.4.2
module load plink/1.9.0
```

Use plink to estimate runs of homozygosity in each individual

```{bash, eval = FALSE}
plink --bfile /projects/course_1/people/bkl835/RoH/Data/Ame_imp_1240ksites --homozyg --homozyg-kb 500 --homozyg-gap 100 --homozyg-density 50 --homozyg-snp 50 --homozyg-window-het 1 --homozyg-window-snp 50 --homozyg-window-threshold 0.05 --out Ame_imp_1240ksites
```

You can check the online manual to see what each parameter means: https://www.cog-genomics.org/plink/1.9/ibd

Now we plot the results in R

```{bash, eval = FALSE}
R

#read plink output
a<-read.table("Ame_imp_1240ksites.hom", as.is=T, h=T)

#compute run lengths
a$l<-a[,8]-a[,7]

#stratify runs per individual
bounds<-c(12000000, 20000000, 3000000000)

inds<-unique(a$IID)

indv<-NULL
lv<-NULL

for(i in inds){
	stratl<-NULL
	stratl<-c(stratl, sum(a[a$IID==i & a$l>=8000000 & a$l<bounds[1], ]$l))
	for(j in 2:length(bounds)){
		stratl<-c(stratl, sum(a[a$IID==i & a$l>=bounds[j-1] & a$l<bounds[j], ]$l))
	}
	indv<-c(indv, rep(i, length(stratl)))
	lv<-c(lv, cumsum(stratl))
}

stratroh<-data.frame(i=indv, l=lv, stringsAsFactors=F)

#Reorder our data before plotting
stratroh$rohl<-c("8-12Mb", "12-20Mb", ">20Mb")
stratroh$p<-sapply(strsplit(stratroh[,1], "_"), "[[", 1)

stratroh$rohl<-factor(stratroh$rohl, levels=rev(c("8-12Mb", "12-20Mb", ">20Mb")))

orderedids<-NULL

for(i in c("Yoruba", "Papuan", "French", "Icelandic", "Punjabi", "Han", "Anzick1", "SpCave", "Lovelock", "Lovelock3", "Taino", "Mixe", "SuruiHGDP", "Surui", "KaritianaHGDP")){
	orderedids<-c(orderedids, unique(stratroh$i[stratroh$p==i]))
}

stratroh$i<-factor(stratroh$i, levels=orderedids)

#Plot

library(ggplot2)

pdf("RoH.pdf", height=5, width=10)
ggplot(stratroh, aes(fill=rohl, y=l, x=i))+
geom_bar(position="stack", stat="identity", color="black", size=.15)+
ylab("Mean cumulative length\nof RoH (b)")+
theme(axis.text.x = element_text(angle = 90, hjust = 1, colour='black'), axis.text.y = element_text(colour='black'))
dev.off()

q("no")
```


Discuss with your partner:

1. What features of the data could affect our inference and in which direction? For example, how would low depth of coverage change our inferred runs? 
2. What can you say about the populations in the example in terms of their recent and deep evolutionary history?
3. Is there anything in the method/filters that makes you trust the results more/less?




















