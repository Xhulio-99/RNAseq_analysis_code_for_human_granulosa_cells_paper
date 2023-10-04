#!/bin/bash
## SLURM directives
#SBATCH --job-name=Trim_$RUN
#SBATCH -N 1
#SBATCH -n 18
#SBATCH -p g100_usr_prod
#SBATCH --mem=40GB
#SBATCH --time 00:30:00
#SBATCH --account PHD_Dhori
#SBATCH --error trimm.err
#SBATCH --output trimm.out

# Module(s) loading
module load profile/bioinf
module load autoload trimmomatic

PATH_TO_RUN=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}/$SAMN/${RUN}
# Command line
cd $PATH_TO_RUN
java -jar $TRIMMOMATIC_HOME/bin/trimmomatic-0.39.jar PE \
        -threads 18 \
        -phred33 \
        ${RUN}_1.fastq ${RUN}_2.fastq \
        trimmed_${RUN}_R1_paired.fq.gz \
        trimmed_${RUN}_R1_unpaired.fq.gz \
        trimmed_${RUN}_R2_paired.fq.gz \
        trimmed_${RUN}_R2_unpaired.fq.gz \
        ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
