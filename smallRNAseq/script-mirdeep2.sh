#!/bin/bash

# SLURM directives
#SBATCH --job-name=mirdeep2
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p g100_usr_prod
#SBATCH --account PHD_Dhori 
#SBATCH -p g100_usr_prod
#SBATCH --mem=40GB
#SBATCH --time 10:00:00
#SBATCH --error mapping.err
#SBATCH --output mapping.out

# Module(s) loading
module load autoload profile/bioinf
module load illumina_genome_Homo_sapiens/hg38

FASTAS=/g100/home/userexternal/xdhori00/data/fasta
cd $PATH_TO_PRJ/collapsed_read

miRDeep2.pl reads.fa $SEQUENCE/WholeGenomeFasta/genome.fa reads_vs_genome.arf \
$FASTAS/mmirBase-hsa.fa none $FASTAS/pmirBase-hsa.fa -P -t Human 2>report.log