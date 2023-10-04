#!/bin/bash

# SLURM directives
#SBATCH --job-name=mir_mapper
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p g100_usr_prod
#SBATCH --account PHD_Dhori 
#SBATCH -p g100_usr_prod
#SBATCH --mem=40GB
#SBATCH --time 01:00:00
#SBATCH --error mapping.err
#SBATCH --output mapping.out

# Module(s) loading
module load autoload profile/bioinf
module load illumina_genome_Homo_sapiens/hg38
#PATH_TO_SAMPLE=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}/$SAMN/${SAMPLE}

# Commands

cd $PATH_TO_PRJ/collapsed_read

mapper.pl $PATH_TO_PRJ/*conf.txt -d -c -m -o ${SLURM_CPUS_PER_TASK} -p $BOWTIE_INDEX/genome \
-s reads.fa -t reads_vs_genome.arf

