#!/bin/bash

# SLURM directives
#SBATCH --job-name=sorting
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

PATH_TO_RUN=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}/$SAMN/${RUN}
cd $PATH_TO_RUN
# Add here the command line:
samtools sort -@ 8 $RUN.sam -o $RUN.bam
