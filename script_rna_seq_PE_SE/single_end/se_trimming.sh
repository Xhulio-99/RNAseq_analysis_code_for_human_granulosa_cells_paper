#!/bin/bash
## SLURM directives
#SBATCH --job-name=Trim
#SBATCH -N 1
#SBATCH -n 18
#SBATCH -p g100_usr_prod
#SBATCH --mem=40GB
#SBATCH --time 01:00:00
#SBATCH --account PHD_Dhori
#SBATCH --error trimm.err
#SBATCH --output trimm.out

# Module(s) loading
module load profile/bioinf
module load autoload trimmomatic

PATH_TO_RUN=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}/$SAMN/${RUN}
# Command line
cd $PATH_TO_RUN
java -jar $TRIMMOMATIC_HOME/bin/trimmomatic-0.39.jar SE \
        -threads 18 \
        -phred33 \
        ${RUN}.fastq \
        trimmed_${RUN}.fq.gz \
        ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
