# LDcor
LDcor is a linkage diseuilibrium statitic, which can be used to assessing the evilutionary potential of species. 
This uses individual-level genotype of populations to compute and compare global linkage disequilibrium of different mutation sets.
Almost all of the steps is command-line based.
The relative singularity container is on the way.

## Prerequisites
The software is develoged and tested in Linux environment.

To infer and add ancestral allele of specified VCF file, you need relative chain files (which can be obtained by [Last 961]) and [liftOver].

To prepare before exact calculation of LDcor from VCF files (with ancestral allele added in the INFO column), you need [VCFtools 0.1.16], [bcftools 1.18], [vep 110], and [Python 3.10.12], the localization of [vep 110] can be replaced by [singularity 3.8.4].

To calculte LDcor, you need [R 4.3.1], with the following libraries:
* [boot 1.3_28.1]
* [matrixStats 1.2.0 0.61.0] 
* [BSDA 1.2.2]
* [data.table 1.14.8]

To reproduce the pictures presented in paper, you need the following libraries in addition:
* [ggplot2 3.4.3]
* [grid 4.3.1]
* [gridExtra 2.3]
* [dplyr 1.1.3]
* [ggsci 3.0.1]
* [RColorBrewer]
* [cowplot 1.1.3]
* [ggalt 0.4.0]
* [ggpubr 0.6.0]
* [magrittr 2.0.3]

To reproduce the simulations presented in paper, you need [SLiM 3.6].

## Useful Input Data

The calculation based on population-level VCF file with ancestral allele added in the INFO column.

## Project Layout
When you prepared an environment with "Prerequisites" met and input data, you can follow the command lines listed in "the_command_list.sh".

## Source data of this study
The data and scripts for figures of the related paper are collected in "paper_fig_code/".

