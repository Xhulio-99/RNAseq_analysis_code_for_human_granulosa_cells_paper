job_id1=$(sbatch --export=SAMPLE=$SAMPLE --error=$SAMPLE/$SAMPLE.fastqc.err \
        --output=$SAMPLE/$SAMPLE.fastqc.out fastqc.sh)
echo "FastQC Job ID is:" $job_id1

job_id2=$(sbatch -J trim_$SAMPLE --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.trim.err \
        --output=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.trim.out --parsable se_trimming.sh)
echo "Trimming Job ID is:" $job_id2

job_id3=$(sbatch -J map_$SAMPLE --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.map.err \
        --output=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.map.out --depend=afterok:$job_id2 --parsable se_mapping.sh)
echo "Alignment Job ID is:" $job_id3

#Add here the submission of sorting.sh. Remember that this job can start if mapping.sh finishes without errors:
job_id4=$(sbatch -J sort_$SAMPLE --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.sort.err \
        --output=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.sort.out --parsable --depend=afterok:$job_id3 sorting.sh)
echo "Sorting Job ID is:" $job_id4
        

job_id5=$(sbatch -J merge_$SAMN --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMN.merge.err \
        --output=$PATH_TO_SAMN/$SAMN.merge.out --parsable merge_samples.sh)
echo "Merge Job ID is:" $job_id5

#Add here the submission of stringtie.sh. Remember that this job can start if sorting.sh finishes without errors:
job_id6=$(sbatch -J stringtie_$SAMN --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMN.stringtie.err \
        --output=$PATH_TO_SAMN/$SAMN.stringtie.out --parsable --depend=afterok:$job_id4 no_merge_stringtie.sh)
echo "StringTie Job ID is:" $job_id6