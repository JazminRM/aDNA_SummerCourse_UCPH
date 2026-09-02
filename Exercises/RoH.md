## Estimating Runs of Homozygocity with ancient genomes

### Interactive node 
We start by starting an interactive job:

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# request one CPU using salloc like this:
salloc --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000  --reservation=3685-26-00-00 --account=teaching

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
plink --bfile /projects/course_1/people/clx746/DataRelatedness/2_plink/koszyce_1240k  --homozyg --out koszyce_1240k

```

You can check the online manual to see what each parameter means: https://www.cog-genomics.org/plink/1.9/ibd

Now we plot the results in R

```{bash, eval = FALSE}
R

#read plink output
a<-read.table("koszyce_1240k.hom", as.is=T, h=T)

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

library(ggplot2)

pdf("RoH.pdf", height=7, width=10)
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




















