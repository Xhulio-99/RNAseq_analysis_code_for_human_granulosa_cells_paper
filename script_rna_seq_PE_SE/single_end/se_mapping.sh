#!/bin/bash

# SLURM directives
#SBATCH --job-name=alignment
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p g100_usr_prod
#SBATCH --account PHD_Dhori 
#SBATCH -p g100_usr_prod
#SBATCH --mem=50GB
#SBATCH --time 01:00:00
#SBATCH --error mapping.err
#SBATCH --output mapping.out

# Module(s) loading
module load profile/bioinf
module load autoload hisat2
module load illumina_genome_Homo_sapiens/


PATH_TO_RUN=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}/${SAMN}/$RUN

# Commands
cd $PATH_TO_RUN

hisat2 -p 8 -q \
-x $HISAT2_INDEX/genome \
-U trimmed_${RUN}.fq.gz \
-S ${RUN}.sam
