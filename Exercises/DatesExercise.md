## Admixture dating with ancient samples

For this exercise we are going to try to date the admixture event that formed the genomic ancestry of Greenlandic Inuit. 

We have a diploid genotype panel that includes the genomes of 2 present-day Inuit and data from some populations that we are going to use as potential sources of admixture. 

- Present-day Mixe from Mexico. We know that around one third of the ancestry in Inuit is related to Indigenous Americans. 
- Present-day Yoruba genomes
- Present-day French genomes
- 3 genomes of 7.6ka old individuals from the Devil's Gate Cave in East Asia
- The genome of a 9ka old individual from the Kolyma River in Siberia
- The genome of a 4ka old individual from the PalaeoInuit Saqqaq culture in Greenland

For all individuals in the panel, we inferred their diploid genotypes using imputation and considered a subset of ~1.2M SNPs across the genome. 

We are going to use DATES to build admixture linkage disequilibrium curves (ALD) and estimate admixture dates for the present-day Inuit. 

We will model the ancestry of present-day Inuit as a mixture of Indigenous American (as represented by the Mixe) and one of the other five potential sources. 

Since the program takes a long time to run, we will divide the group in five groups and each of you will run one combination. 


### Interactive node

We will start by getting an interactive node:

```{bash, eval = FALSE}
# First log in to the server (remember to change ku_username for your username)
ssh ku_username@mjolnirgate.unicph.domain

# first request one CPU using salloc like this:
salloc --partition=cpuqueue --nodes=1 -D `pwd` --mem-per-cpu 5250 --ntasks-per-node=1 -t 1000 --qos=teaching --reservation=aDNA_PHD_course --account=teaching

# once the job has been allocated, you can login to the node with srun like this:
srun --pty -n 1 -c 1 bash -i
```

Once you login to you interactive node, create a directory for today's exercises (make sure to change "username" for your own directory name):
 
```{bash, eval = FALSE}
# create directory
directoryDates="/projects/course_1/people/${username}/Dates/"
mkdir -p $directoryDates
# go into the directory
cd $directoryDates
```

Then, we load dates

```{bash, eval = FALSE}
module dates
```

Data have already been converted into an eigenstrat file. 

We start by deciding which source population we will test. For this example we will use French. The options are: Saqqaq Kolyma DevilsCaveN French Yoruba

```{bash, eval = FALSE}
sourceB="French"
```

Remember that for all models, sourceA will be Mixe and we are only looping over different sourceB populations. 

Then we create a population list file that describes our model

```{bash, eval = FALSE}
echo "Mixe "$sourceB" Inuit Mixe_"$sourceB"_Inuit" >  Mixe_"$sourceB"_Inuit_poplist
```

Now we create an input file for dates:

```{bash, eval = FALSE}
echo "genotypename: /projects/course_1/people/bkl835/Dates/Data/Ame_imp_1240ksites_gm.eigenstratgeno" > Mixe_"$sourceB"_Inuit_par
echo "snpname: /projects/course_1/people/bkl835/Dates/Data/Ame_imp_1240ksites_gm.snp" >> Mixe_"$sourceB"_Inuit_par
echo "indivname: /projects/course_1/people/bkl835/Dates/Data/Ame_imp_1240ksites_gm.ind" >> Mixe_"$sourceB"_Inuit_par
echo "admixlist: Mixe_"$sourceB"_Inuit_poplist" >> Mixe_"$sourceB"_Inuit_par
echo "minparentcount: 1" >> Mixe_"$sourceB"_Inuit_par
echo "maxdis: 0.5" >> Mixe_"$sourceB"_Inuit_par
echo "jackknife: YES" >> Mixe_"$sourceB"_Inuit_par
echo "runfit: YES" >> Mixe_"$sourceB"_Inuit_par
echo "afffit: YES" >> Mixe_"$sourceB"_Inuit_par
echo "lovalfit: 0.45" >> Mixe_"$sourceB"_Inuit_par
```

Finally we run dates

```{bash, eval = FALSE}
dates -p Mixe_"$sourceB"_Inuit_par > Mixe_"$sourceB"_Inuit_log
```

Dates created new directory. Inside, you can find multiple files. 

The jout file has the admixture date estimates and standard errors (in generations).


The other important file is the pdf file which contains the ALD curve. 


1. Sit together with your team and compare and contrast your different ALD curves. How are they different? Do all of them show a signal of ALD? (Important: notice that y-axis scale could be different between plots)
2. Which one is the sourceB population that represents the non-Mixe-related ancestry in Inuit best? 
3. Look at the admixture date estimates for the best combination. Does this make sense in terms of when Inuit peoples move into the Americas? (you can check online)






