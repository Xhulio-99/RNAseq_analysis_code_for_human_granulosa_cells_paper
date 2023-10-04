#!/bin/bash

# SLURM directives
#SBATCH --job-name=stringTie
#SBATCH -N 1 
#SBATCH -n 8
#SBATCH -p g100_usr_prod
#SBATCH --mem=20GB
#SBATCH --time 1:00:00
#SBATCH --account PHD_Dhori
#SBATCH --error stringTie.err
#SBATCH --output stringTie.out

# Module(s) loading
#add here the modules to load:
module load profile/bioinf
module load autoload stringtie
module load illumina_genome_Homo_sapiens/


PATH_TO_RUN=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}/$SAMN/${RUN}
cd $PATH_TO_RUN
# Add here the command line:
stringtie $RUN.bam \
-B -e -p 8 \
-G ${GENES}/genes.gtf \
-o $RUN.gtf

