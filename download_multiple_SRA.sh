#!/bin/bash
#SBATCH --job-name=Download_mul
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p g100_all_serial
#SBATCH --mem=30GB
#SBATCH --time 04:00:00
#SBATCH --error Download_all_SE.err
#SBATCH --output Download_all_SE.out
#SBATCH --account PHD_Dhori

##################################################
############### Module(s) loading ################
##################################################
module load autoload profile/bioinf
module load autoload sra/3.0.0
    
#Download Single End Sample from SRA
cd $SCRATCH/original_fastq/$PRJ

start=`date +%H%s`

prefetch -O . --option-file SRR_Acc_List.txt

end=`date +%H%s`
echo Execution time was `expr $end - $start` seconds
