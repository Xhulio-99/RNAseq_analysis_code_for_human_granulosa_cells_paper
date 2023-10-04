#!/bin/bash
#SBATCH --job-name=Download_mul
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p g100_usr_prod
#SBATCH --mem=100GB
#SBATCH --time 05:00:00
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

cat SRR_Acc_List.txt | xargs -I {} fasterq-dump -e ${SLURM_CPUS_ON_NODE} -O {} {}

end=`date +%H%s`
echo Execution time was `expr $end - $start` seconds
