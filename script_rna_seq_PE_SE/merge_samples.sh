#!/bin/bash

# SLURM directives
#SBATCH --job-name=merging
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH -p g100_usr_prod
#SBATCH --mem=30GB
#SBATCH --time 1:00:00
#SBATCH --account PHD_Dhori
#SBATCH --error sorting.err
#SBATCH --output sorting.out

# Module(s) loading
#add here the modules to load:
module load profile/bioinf
module load autoload samtools


# Add here the command line:
PATH_TO_SAMN=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}/${SAMN}
cd $PATH_TO_SAMN
ls SRR*/*.bam | xargs samtools merge -@ 8 -o $SAMN.bam
