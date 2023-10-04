#!/bin/bash

# SLURM directives
#SBATCH --job-name=mir_mapper
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p g100_usr_prod
#SBATCH --account PHD_Dhori 
#SBATCH -p g100_usr_prod
#SBATCH --mem=100GB
#SBATCH --time 04:00:00
#SBATCH --error mapping.err
#SBATCH --output mapping.out

# Module(s) loading
module load autoload profile/bioinf
module load illumina_genome_Homo_sapiens/hg38

# Commands
mkdir -p collapsed_read
cd $PATH_TO_PRJ/collapsed_read

Ampliseq_adapter=CTGTCTCTTATACACATCT
Truseq_adapter=TGGAATTCTCGGGTGCCAAGG

mapper.pl $PATH_TO_PRJ/*conf.txt -d -e -i -j -h -l 18 -m -o ${SLURM_CPUS_ON_NODE} \
-k TCGTATGCCGTCTTCTGCTTGT -p $BOWTIE_INDEX/genome \
-s reads.fa -t reads_vs_genome.arf