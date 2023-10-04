#!/bin/bash

PRJ=$1
PATH_TO_PRJ=${SCRATCH}/original_fastq/${PRJ}
for SAMN in `ls ${SCRATCH}/original_fastq/${PRJ} | grep "^SAMN"`;do
        PATH_TO_SAMN=$PATH_TO_PRJ/$SAMN
        echo "Running Samn: " $SAMN
        for SAMPLE in `ls $PATH_TO_SAMN | grep "SRR"`;do
        echo "Running Sample: " $SAMPLE
        
                #sbatch --export=SAMPLE=$SAMPLE --error=$SAMPLE/$SAMPLE.fastqc.err --output=$SAMPLE/$SAMPLE.fastqc.out fastqc.sh

                job_id2=$(sbatch -J trim_$SAMPLE --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.trim.err \
                        --output=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.trim.out --parsable se_trimming.sh)
                echo "Trimming ID is:" $job_id2

                job_id3=$(sbatch -J map_$SAMPLE --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.map.err \
                        --output=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.map.out --depend=afterok:$job_id2 --parsable se_mapping.sh)
                echo "alignment ID is:" $job_id3

                #Add here the submission of sorting.sh. Remember that this job can start if mapping.sh finishes without errors:
                job_id4=$(sbatch -J sort_$SAMPLE --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.sort.err \
                        --output=$PATH_TO_SAMN/$SAMPLE/$SAMPLE.sort.out --parsable --depend=afterok:$job_id3 sorting.sh)
                echo "sorting ID is:" $job_id4
        
        done
        : '
        job_id5=$(sbatch -J merge_$SAMN --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMN.merge.err \
                --output=$PATH_TO_SAMN/$SAMN.merge.out --parsable merge_samples.sh)
        echo "merge ID is:" $job_id5
        '

        #Add here the submission of stringtie.sh. Remember that this job can start if sorting.sh finishes without errors:
        job_id6=$(sbatch -J stringtie_$SAMN --export=SAMPLE=$SAMPLE,PRJ=$PRJ,SAMN=$SAMN --error=$PATH_TO_SAMN/$SAMN.stringtie.err \
                --output=$PATH_TO_SAMN/$SAMN.stringtie.out --parsable --depend=afterok:$job_id4 no_merge_stringtie.sh)
        echo "stringtie ID is:" $job_id6
done




