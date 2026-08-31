# aDNA_course_2026
Genotype imputation exercise for the 2026 aDNA course "Decoding the past: An introduction to ancient genomics"

### Introduction

This tutorial will show you how to run genotype imputation using GLIMPSE2 ([Rubinacci et al.,2021](https://doi.org/10.1038/s41588-020-00756-0), [Rubinacci et al.,2023](https://www.nature.com/articles/s41588-023-01438-3)), an imputation and phasing method tailored for low coverage sequencing data. We will test imputation accuracy in an ancient dog from Siberia, dated to 2,000 years before present (ybp). This sample is relatively high coverage for ancient standards (11.8x), so in order to test the imputation accuracy, we will downsample it to lower coverages (0.1x and 1x) and see how the imputed genotypes compare to the "true" high coverage ones.

We will run GLIMPSE from sequencing reads data (BAM files) to obtain imputed and phased genotypes.

For this tutorial, we will focus only on chromosome 22.
