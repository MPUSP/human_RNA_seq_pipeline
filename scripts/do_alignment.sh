
sample_name=$1
fastq_file=$2
star_index=$3

num_cores=$4
#replace PATH with the path to the STAR executable, i.e. /local/tools/star/2.6.1d/source/STAR-2.6.1d/bin/Linux_x86_64/STAR
PATH --runMode alignReads --runThreadN ${num_cores} --twopassMode Basic --genomeDir ${star_index} --readFilesCommand zcat --readFilesIn ${fastq_file} --outFileNamePrefix ${sample_name}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 alignMatesGapMax 1000000 --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted --outSAMattrRGline ID:rg1 SM:sm1 --outSAMattributes NH HI AS nM NM --alignSoftClipAtReferenceEnds Yes --outSAMstrandField intronMotif
#replace PATH with the path to the samtools executable, i.e. /local/tools/samtools/1.3.1/bin/samtools
PATH sort -@ ${num_cores} -o ${sample_name}_Aligned.out.srt.bam  ${sample_name}_Aligned.out.bam
PATH index ${sample_name}_Aligned.out.srt.bam
