#!/bin/bash

# SLURM directives
#SBATCH --job-name=mir_mapper
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p g100_usr_prod
#SBATCH --account PHD_Dhori 
#SBATCH -p g100_usr_prod
#SBATCH --mem=20GB
#SBATCH --time 00:30:00
#SBATCH --error mapping.err
#SBATCH --output mapping.out

# Commands
PATH_TO_RUN=$PATH_TO_PRJ/$SAMPLE/$RUN
cd $PATH_TO_RUN
mkdir -p $PATH_TO_PRJ/collapsed_read

mapper.pl ${RUN}.fastq -e -h -i -j -k TGGAATTCTCGGGTGCCAAGG \
-l 18 -m -s $PATH_TO_PRJ/collapsed_read/${RUN}.fa
