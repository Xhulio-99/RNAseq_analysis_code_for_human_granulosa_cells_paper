cd $SCRATCH/original_fastq

for PRJ in `ls $SCRATCH/original_fastq/ | grep "PRJ"`;do

    esearch -db sra -query $PRJ | efetch -format runinfo | \
    cut -d ',' -f 1 | grep SRR > $PRJ/SRR_Acc_List.txt

    esearch -db sra -query $PRJ | efetch -format runinfo | \
    cut -d ',' -f 1,26 | grep SRR > $PRJ/SRR_SAMN_List.txt
done