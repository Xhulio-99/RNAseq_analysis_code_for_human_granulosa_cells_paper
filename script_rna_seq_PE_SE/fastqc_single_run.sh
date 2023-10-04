#!/bin/bash
#SBATCH --job-name=Quality_PRJNA417973
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p g100_usr_prod
#SBATCH --mem=20GB
#SBATCH --time 00:30:00
#SBATCH --error fastqc_PRJNA417973.err
#SBATCH --output fastqc_PRJNA417973.out
#SBATCH --account PHD_Dhori

module load autoload profile/bioinf
module load autoload fastqc/0.11.9

PRJ=$1
PATH_TO_PRJ=/g100_scratch/userexternal/xdhori00/original_fastq/${PRJ}
cd $PATH_TO_PRJ
mkdir -p fastqc_out

while getopts ":t" option; do
   case $option in
        t) # fastqc on trimmed fastq
            mkdir -p fastqc_out/trimmed
            fastqc -t 10 --nogroup --extract -o fastqc_out/trimmed $PATH_TO_PRJ/SAMN*/SRR*/trim*.fq.gz
            exit;;
   esac
done

fastqc -t 10 --nogroup --extract -o fastqc_out $PATH_TO_PRJ/SAMN*/SRR*/*.fastq