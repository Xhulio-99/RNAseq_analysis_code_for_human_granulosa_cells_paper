for i in `ls -d PRJ*`;do
    cd $SCRATCH/original_fastq/$i
    prepDE.py -g gene_count_matrix_$i.csv -t transcript_count_matrix_$i.csv
done