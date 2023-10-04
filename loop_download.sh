for PRJ in `ls $SCRATCH/original_fastq/ | grep "PRJ"`;do

    job_id1=$(sbatch -J download_$PRJ --export=PRJ=$PRJ --error=$SCRATCH/original_fastq/$PRJ/download.err \
    --output=$SCRATCH/original_fastq/$PRJ/download.out --parsable \
    download_multiple_SRA.sh)
    echo Prefetching Job ID:$job_id1

    job_id2=$(sbatch -J dump_$PRJ --export=PRJ=$PRJ --error=$SCRATCH/original_fastq/$PRJ/dump.err \
    --output=$SCRATCH/original_fastq/$PRJ/dump.out --parsable --depend=afterok:$job_id1\
    fasterq-dump.sh)
    echo Dumping Job ID:$job_id2

done
