#!/bin/bash
#SBATCH --job-name=Fastqc
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p g100_usr_prod
#SBATCH --mem=20GB
#SBATCH --time 01:00:00
#SBATCH --error fastqc_err.err
#SBATCH --output fastqc_out.out
#SBATCH --account PHD_Dhori

module load autoload profile/bioinf
module load autoload fastqc/0.11.9

cd $PATH_TO_PRJ
mkdir -p fastqc_out
start=`date +%H%s`
fastqc -t 10 --nogroup --extract -o fastqc_out $SAMN*/SRR*/trim*.fq
end=`date +%H%s`
echo Execution time for R1 was `expr $end - $start` seconds