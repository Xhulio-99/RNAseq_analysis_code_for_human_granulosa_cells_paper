#!/bin/bash

# SLURM directives
#SBATCH --job-name=mir_quantifier
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p g100_usr_prod
#SBATCH --account PHD_Dhori 
#SBATCH -p g100_usr_prod
#SBATCH --mem=70GB
#SBATCH --time 10:00:00
#SBATCH --error mapping.err
#SBATCH --output mapping.out


cd $PATH_TO_PRJ/collapsed_read
FASTAS=/g100/home/userexternal/xdhori00/data/fasta

quantifier.pl -T ${SLURM_CPUS_ON_NODE} -p $FASTAS/pmirBase-hsa.fa -m $FASTAS/mmirBase-hsa.fa \
-r reads.fa -t Human