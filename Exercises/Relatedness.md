## Estimating relatedness with ancient sequencing data

### Outline:

* Estimating relatedness in the Koszyce invididuals
* Estimating relatedness is in ancient wolves
* Estimating relatedness in ancient grapes


### Interactive node

Get an interactive node: 

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# first request one CPU using salloc like this:
salloc --partition=cpuqueue --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000 --qos=teaching --reservation=aDNA_PHD_course --account=teaching

# once the job has been allocated, you can login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```

Now, let's create a directory for these exercises:
```{bash, eval = FALSE}
# remember you change for you username
username="write_your_username"
directory="/projects/course_1/people/${username}/relatedness/"
mkdir -p $directory

# go to the directory
cd $directory
```

### Estimating relatedness in the Koszyce invididuals

#### READ

`READ` <sup>1</sup> is a method that was developed specifically for low coverage data, so it requires pseudo-haploid calls like the ones we have been using in previous days. 

**Dataset**

For this session we will be working on the data from Schroeder et al. 2019 PNAS <sup>2</sup>.

We'll be working with these PLINK files:
```
/projects/course_1/people/clx746/DataRelatedness/2_plink/koszyce_1240k_subset.bed
/projects/course_1/people/clx746/DataRelatedness/2_plink/koszyce_1240k_subset.fam
/projects/course_1/people/clx746/DataRelatedness/2_plink/koszyce_1240k_subset.bim
```
These contain SNP data for 5 related individuals.  

```{bash, eval=FALSE}
cat /projects/course_1/people/clx746/DataRelatedness/2_plink/koszyce_1240k_subset.fam
```
```
RISE1160_rnd RISE1160_rnd 0 0 0 -9
RISE1161_rnd RISE1161_rnd 0 0 0 -9
RISE1164_rnd RISE1164_rnd 0 0 0 -9
RISE1168_rnd RISE1168_rnd 0 0 0 -9
RISE1169_rnd RISE1169_rnd 0 0 0 -9
```
Now let's see how we can run ```READ```. 

```READ``` will automatically assign a default name for the output, so we have to change the output name once it is done running.  

You can run it like this (we won't run it this time since our AWS servers don't have enouth RAM for it to work, so just take a look at how it is run and download the output after): 
```
# load python
module unload python
module load python/2.7.17

# run read:
python /projects/course_1/people/clx746/Data/read/READ.py /projects/course_1/people/clx746/DataRelatedness/2_plink/koszyce_1240k_subset
```

#### KIN


```{bash, eval=FALSE}
module load kin/3.1.4

KINgaroo

```

 
