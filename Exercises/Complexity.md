## Estimating library complexity 

### Outline: 

* Interactive node
* Estimating complexity with Preseq
* Designing sequencing strategy


### Interactive node

Start by opening an interactive session in the server (if you are not in one already).

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# request one CPU using salloc like this:
salloc --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000  --reservation=3685-26-00-00 --account=teaching

# once the job has been allocated, you can login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```
Go to your mapping directory:

```{bash, eval = FALSE}
# remember to change ku_username for your username:
cd /projects/course_1/people/ku_username/mapping/
```
### Estimating complexity with Preseq

Preseq can estimate a complexity curve to estimate how many unique sequences you can obtain if you sequence your library deeper. 

We will use Preseq to estimate a complexity curve for each of our libraries in our ancient genomes: 

```{bash, eval = FALSE}
# load preseq
module load preseq/3.1.2 
```

As input for preseq, we need to provide the BAM files with the reads mapped before we removed the PCR duplicates, so we can run Preseq like this:
```{bash, eval = FALSE}
# Extract 1, library 1
preseq lc_extrap -B PW13_E1_L1.mapped.sorted.bam > PW13_E1_L1_comp
# Extract 2, library 1
preseq lc_extrap -B PW13_E2_L1.mapped.sorted.bam > PW13_E2_L1_comp
# Extract 2, library 2
preseq lc_extrap -B PW13_E2_L2.mapped.sorted.bam > PW13_E2_L2_comp
```

You can plot the results like this (this is an example for E2_L2, but try to use this code to plot the results for all libraries):
```{R, eval = FALSE}
R
# output from preseq:
filename<-"PW13_E2_L2_comp"

# Total reads currently mapped before removing duplicates 
seqreads<-358736

# Total reads currently mapped after removing duplicates:
distinctreads<-357775 

a<-read.table(filename, h=T)

pdf("PW13_E2_L2_comp.pdf", height=7, width=7)
# Here modify the xlim and ylim accordingly:
plot(a[,1], a[,2], t="l", xlab="Sequenced reads", ylab="Expected distinct reads", col="slateblue4", xlim=c(0, 500000000), ylim=c(0, 150000000))
polygon(c(a[,1], rev(a[,1])), c(a[,3], rev(a[,4])), col="#473C8B40", border=NA)

abline(v=seqreads, col="red", lty=2)
abline(h=distinctreads, col="red", lty=2)

legend("topright", c("Current"), pch="+", col=c("red"), bty="n")

dev.off()
```

Question: Given complexity curves, the average read length observed from the mapDamage results and the human genome size (~3200000000 bp), how many more reads do we need to sequence to get to 1X of depth of coverage for each library?


