PRJ=$1
PATH_TO_PRJ=${SCRATCH}/mirSeq/${PRJ}

job_id1=$(sbatch -J mapper_$PRJ --export=PATH_TO_PRJ=$PATH_TO_PRJ --error=$PATH_TO_PRJ/mapper.err \
    --output=$PATH_TO_PRJ/mapper.out --parsable $SCRATCH/mirSeq/script-mapper.sh)
echo "Mapper Job ID is:" $job_id1

job_id2=$(sbatch -J mirdeep2_$PRJ --export=PATH_TO_PRJ=$PATH_TO_PRJ --error=$PATH_TO_PRJ/mirdeep2.err \
    --output=$PATH_TO_PRJ/mirdeep2.out --depend=afterok:$job_id1 --parsable $SCRATCH/mirSeq/script-mirdeep2.sh)
echo "miRdeep2 Job ID is:" $job_id2

job_id3=$(sbatch -J quantifier_$PRJ --export=PATH_TO_PRJ=$PATH_TO_PRJ --error=$PATH_TO_PRJ/quantifier.err \
    --output=$PATH_TO_PRJ/quantifier.out --depend=afterok:$job_id1 --parsable $SCRATCH/mirSeq/script-quantifier.sh)
echo "Quantifier Job ID is:" $job_id3