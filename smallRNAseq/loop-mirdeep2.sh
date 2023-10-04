PRJ=$1
#PATH_TO_PRJ=${SCRATCH}/original_fastq/${PRJ}
PATH_TO_PRJ=${SCRATCH}/mirSeq/${PRJ}
for SAMPLE in `ls $PATH_TO_PRJ | grep "SAMN"`;do
        break
        PATH_TO_SAMPLE=$PATH_TO_PRJ/$SAMPLE
        echo "Running SAMPLE: " $SAMPLE
        for RUN in `ls $PATH_TO_SAMPLE | grep "SRR"`;do
            job_id1=$(sbatch -J mapper_$PRJ --export=RUN=$RUN,PATH_TO_PRJ=$PATH_TO_PRJ,SAMPLE=$SAMPLE --error=$PATH_TO_PRJ/$SAMPLE/$RUN/mapper.err \
                        --output=$PATH_TO_PRJ/$SAMPLE/$RUN/mapper.out --parsable script-mapper1.sh)
                echo "Mapper1 Job ID is:" $job_id1
        done
done
job_id2=$(sbatch -J mapper_$PRJ --export=RUN=$RUN,PATH_TO_PRJ=$PATH_TO_PRJ,SAMPLE=$SAMPLE --error=$PATH_TO_PRJ/mapper2.err \
    --output=$PATH_TO_PRJ/mapper2.out --parsable $SCRATCH/mirSeq/script-mapper2.sh)
echo "Mapper2 Job ID is:" $job_id2

job_id3=$(sbatch -J mirdeep2_$PRJ --export=PATH_TO_PRJ=$PATH_TO_PRJ --error=$PATH_TO_PRJ/mirdeep2.err \
    --output=$PATH_TO_PRJ/mirdeep2.out --depend=afterok:$job_id2 --parsable $SCRATCH/mirSeq/script-mirdeep2.sh)
echo "miRdeep2 Job ID is:" $job_id3