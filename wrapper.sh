#create sample fastq list
#create metadata file


echo `date` "-- Starting Human RNA-Seq Pipeline"
echo "See ./logs for details."
echo `date` "-- Setting envirnonment variables..."
DATA_DIR=`realpath ./example`
#path to star aligner based indexed reference sequence for homo sapiens, i.e. /local/genomes/h_sapiens/star_76bp_index
STAR_INDEX=PATH
#path to annotation gtf for reference, i.e. /local/genomes/h_sapiens/gencode.v26.annotation.genes.gtf
ANNOTATION=PATH

#adjust accordingly
NUM_JOBS=6
NUM_CORES=32

#adjust accordingly
CONTROL_GROUP=""
TREATMENT_GROUP=""


mkdir tmp
TMPDIR=`realpath tmp`
TMP_DIR=`realpath tmp`

echo `date` "-- Done."
mkdir logs
LOGDIR=`realpath logs`
echo "will cite" | parallel --bibtex 1>${LOGDIR}/parallel.out 2>${LOGDIR}parallel.err

echo `date` "-- Building Python virtual environment..."
#replace PATH with your path to your python3 executable, i.e. /opt/rh/rh-python36/root/usr/bin/python3
virtualenv -p PATH rnaseq_env 1>${LOGDIR}/venv.out 2>${LOGDIR}/venv.err
. rnaseq_env/bin/activate
pip install -r venv_requirements.txt  1>${LOGDIR}/venv.out 2>${LOGDIR}/venv.err
echo `date` "-- Done."

echo `date` "-- Building R virtual environment..."
mkdir deseq
cd deseq
RENV_PATHS_ROOT=`realpath .`
Rscript ../scripts/build_renv.R 1>${LOGDIR}/renv.out 2>${LOGDIR}/renv.err
cd renv/library/
tar -xzvf ../../../human_rnaseq_r_libraries.tar.gz 1>${LOGDIR}/renv.out 2>${LOGDIR}/renv.err
cd ../../..
echo `date` "-- Done."
#create alignment commands to run
echo `date` "-- Starting STAR alignment..."
mkdir alignment
cd alignment
awk '{print "../scripts/do_alignment.sh "$1" "$2" "starindex" "cores}' starindex="${STAR_INDEX}" cores="${NUM_CORES}" ${DATA_DIR}/sample_fastq_list.txt > commands_to_run.sh
parallel -j ${NUM_JOBS}  --nice 18 --tmpdir $TMPDIR < commands_to_run.sh 1>${LOGDIR}/star.out 2>${LOGDIR}/star.err
cd ..
echo `date` "-- Done."


#create rnaseqc commands to run
# todo -- make awk FS explicit
echo `date` "-- Starting RNA-SeQC quantification..."
mkdir rnaseqc
cd rnaseqc

cut -f1 ${DATA_DIR}/sample_fastq_list.txt | while read data; do realpath ../alignment/${data}_Aligned.out.srt.bam; done > bamlist.txt
paste ${DATA_DIR}/sample_list.txt bamlist.txt  > sample_bam_list.txt

#cat sample_bam_list.txt | while read data; do sample=`echo $data| awk '{print $1}'`; bam=`echo $data| awk '{print $2}'`; echo $sample; echo $bam; (../scripts/rnaseqc.v2.1.0.linux ../references/gencode.v26.annotation.genes.sorted.gtf ${bam} ${sample} -v -v -u -s ${sample} &) done

awk '{print "../scripts/rnaseqc.v2.1.0.linux " annotation" "$2" "$1" -v -v -u -s "$1}' annotation="${ANNOTATION}"  sample_bam_list.txt > commands_to_run.sh
parallel -j ${NUM_JOBS}  --nice 18 --tmpdir $TMPDIR < commands_to_run.sh 1>${LOGDIR}/rnaseqc.out 2>${LOGDIR}/rnaseqc.err

python ../scripts/join_gct_reads.py */*gene_tpm.gct all_readtpms.tsv 1>>${LOGDIR}/rnaseqc.out 2>>${LOGDIR}/rnaseqc.err
python ../scripts/join_gct_reads.py */*gene_reads.gct all_readcounts.tsv 1>>${LOGDIR}/rnaseqc.out 2>>${LOGDIR}/rnaseqc.err

zcat all_readcounts.tsv.gct.gz| tail -n +3 > all_readcounts.tsv
zcat all_readtpms.tsv.gct.gz | tail -n +3 > all_readtpms.tsv

cd ..
echo `date` "-- Done."

#create deseq commands to run
cd deseq
mkdir svadata

echo `date` "-- Starting DESeq..."
Rscript ../scripts/deseq_preprocessing.R ../rnaseqc/all_readcounts.tsv ${DATA_DIR}/sample_metadata.tsv ${CONTROL_GROUP} ${TREATMENT_GROUP} 1>${LOGDIR}/deseq.out 2>${LOGDIR}/deseq.err
Rscript ../scripts/deseq_postprocessing.R ${CONTROL_GROUP} ${TREATMENT_GROUP} Treatment 1>>${LOGDIR}/deseq.out 2>>${LOGDIR}/deseq.err

echo `date` "-- Done."

echo `date`" -- Starting bokeh plotting..."
python ../scripts/bokeh_volcano_ma.py svadata/${TREATMENT_GROUP}vs.${CONTROL_GROUP}_Treatment.csv ${TREATMENT_GROUP}vs.${CONTROL_GROUP} 1>${LOGDIR}/bokeh.out 2>${LOGDIR}/bokeh.err
echo `date` "-- Done."

cd ..

mkdir output
cd output
ln -s ../deseq/interactive_volcano_ma_${TREATMENT_GROUP}vs.${CONTROL_GROUP}.html
ln -s ../deseq/svadata/rlog_transformed_reads.tsv
ln -s ../deseq/svadata/${TREATMENT_GROUP}vs.${CONTROL_GROUP}_Treatment.csv
ln -s ../rnaseqc/all_readcounts.tsv
ln -s ../rnaseqc/all_readtpms.tsv
ln -s ../rnaseqc/sample_bam_list.txt
cd ..

echo `date` "-- Pipeline complete.  See ./output for results and ./logs for logs"
