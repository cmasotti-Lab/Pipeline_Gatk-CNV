#!/usr/bin/bash

#Pipeline atualizado em 30-10-2023
#1. Ajustes das variaves de ambiente para facilitar na reexecução quando os diretorios forem atualizado

# PARAMETROS OBRIGATORIOS
SCRATCH60="/home/scratch60/vlira_20nov2023/"

#DATA=$(date "+%F") # EDITE SE QUISER USAR UMA PASTA DE UMA DATA ESPECIFICA 
DATA="2023-11-07"
OUTPUT_DIR=$SCRATCH60"/Result_Gatk-CNV.XL."$DATA

INPUT_DIR="/home/scratch60/rtorreglosa_12jan2024/preprocessing_READ_result/"
BAM_FILES=$(find "$INPUT_DIR" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam')
JOBS=5
mem=200
#MAXmem=$((mem / JOBS))
MAXmem=200

COPY_RATIO=5
CUTOFF_AMP=0.5
CUTOFF_DEL=-1

#TOOLS e DATABASES
REF_FASTA="/home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/"
ANNOVAR="$SCRATCH60/tools/annovar/table_annovar.pl"
ANNOVAR_DB="$SCRATCH60/humandb/"
GATK="$SCRATCH60/tools/gatk-4.3.0.0/./gatk"
TARGET="$SCRATCH60/references/xgen-exome-research-panel-v2-targets-hg38.autossome.bed"
BLACKLIST="$SCRATCH60/references/CNV_and_centromere_blacklist.hg38liftover.list"
PON="/home/users/vlira/PanelOfNornal/PON.100COVID.100-eigensamples.XL.hdf5"

mkdir $OUTPUT_DIR

find "$INPUT_DIR" -maxdepth 1 -mindepth 1  -name '*.dedup.tags.bqsr.bam' | grep -Pv "ROP-25-|ROP-26-|ROP-27-|ROP-29-" > $OUTPUT_DIR/samples.list
#head -3 $OUTPUT_DIR/samples.list >  $OUTPUT_DIR/TOY.samples.list
OUTPUT_LOG="$OUTPUT_DIR.log"

export OUTPUT_DIR
export OUTPUT_LOG
export REF_FASTA
export ANNOVAR
export ANNOVAR_DB
export GATK
export PON
export TARGET
export BLACKLIST
export MAXmem
export COPY_RATIO
export CUTOFF_AMP
export CUTOFF_DEL

step1_PreprocessIntervals (){
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step1_PreprocessIntervals <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G"  PreprocessIntervals \
    -L ${TARGET} \
    -XL ${BLACKLIST} \
    -R ${REF_FASTA}/Homo_sapiens_assembly38.fasta \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $OUTPUT_DIR/step1_PreprocessIntervals/targets.preprocessed.interval_list 2> $OUTPUT_DIR/step1_PreprocessIntervals/step1_PreprocessIntervals.log
}
export -f step1_PreprocessIntervals

step2_AnnotateIntervals (){
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step2_AnnotateIntervals <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" AnnotateIntervals \
    -R ${REF_FASTA}/Homo_sapiens_assembly38.fasta \
    -L $OUTPUT_DIR/step1_PreprocessIntervals/targets.preprocessed.interval_list \
    -XL $BLACKLIST \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $OUTPUT_DIR/step2_AnnotateIntervals/annotated_intervals.tsv  2> $OUTPUT_DIR/step2_AnnotateIntervals/step2_AnnotateIntervals.log
}
export -f step2_AnnotateIntervals

step3_CollectReadCounts (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step3_CollectReadCounts para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" CollectReadCounts \
    -I $SAMPLE  \
    -L $OUTPUT_DIR/step1_PreprocessIntervals/targets.preprocessed.interval_list \
    -XL $BLACKLIST \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O  $OUTPUT_DIR/step3_CollectReadCounts/${NAME}.counts.hdf5  2> $OUTPUT_DIR/step3_CollectReadCounts/$NAME.log
}
export -f step3_CollectReadCounts

step4_DenoiseReadCounts (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step4_DenoiseReadCounts para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" DenoiseReadCounts \
    -I $OUTPUT_DIR/step3_CollectReadCounts/${NAME}.counts.hdf5 \
    --annotated-intervals $OUTPUT_DIR/step2_AnnotateIntervals/annotated_intervals.tsv \
    --count-panel-of-normals ${PON} \
    --standardized-copy-ratios $OUTPUT_DIR/step4_DenoiseReadCounts/${NAME}.standardizedCR.tsv \
    --denoised-copy-ratios $OUTPUT_DIR/step4_DenoiseReadCounts/${NAME}.denoisedCR.tsv  2> $OUTPUT_DIR/step4_DenoiseReadCounts/$NAME.log
}
export -f step4_DenoiseReadCounts

step5_PlotDenoisedCopyRatios (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step5_PlotDenoisedCopyRatios para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" PlotDenoisedCopyRatios \
    --standardized-copy-ratios $OUTPUT_DIR/step4_DenoiseReadCounts/${NAME}.standardizedCR.tsv \
    --denoised-copy-ratios $OUTPUT_DIR/step4_DenoiseReadCounts/${NAME}.denoisedCR.tsv \
    --sequence-dictionary ${REF_FASTA}/Homo_sapiens_assembly38.dict \
    --minimum-contig-length 46709983 \
    --output $OUTPUT_DIR/step5_PlotDenoisedCopyRatios/ \
    --output-prefix ${NAME}  2> $OUTPUT_DIR/step5_PlotDenoisedCopyRatios/$NAME.log
}
export -f step5_PlotDenoisedCopyRatios

step5_PlotDenoisedCopyRatios_DOCKER (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step5_PlotDenoisedCopyRatios_DOCKER para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  docker run \
    --rm \
    -v $OUTPUT_DIR/:/gatk/my_data \
    -v ${REF_FASTA}:/gatk/hg38/ \
    -u $(id -u):$(id -g) broadinstitute/gatk \
    gatk --java-options "-Xmx${MAXmem}G" PlotDenoisedCopyRatios \
    --standardized-copy-ratios /gatk/my_data/step4_DenoiseReadCounts/${NAME}.standardizedCR.tsv \
    --denoised-copy-ratios /gatk/my_data/step4_DenoiseReadCounts/${NAME}.denoisedCR.tsv \
    --sequence-dictionary /gatk/hg38//Homo_sapiens_assembly38.dict \
    --minimum-contig-length 46709983 \
    --output /gatk/my_data/step5_PlotDenoisedCopyRatios/ \
    --output-prefix ${NAME}  2> $OUTPUT_DIR/step5_PlotDenoisedCopyRatios/$NAME.log
}
export -f step5_PlotDenoisedCopyRatios_DOCKER


step6_CollectAllelicCounts (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step6_CollectAllelicCounts para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" CollectAllelicCounts \
    -L $OUTPUT_DIR/step1_PreprocessIntervals/targets.preprocessed.interval_list \
    -I $SAMPLE \
    -R ${REF_FASTA}/Homo_sapiens_assembly38.fasta \
    -O $OUTPUT_DIR/step6_CollectAllelicCounts/${NAME}.allelicCounts.tsv  2> $OUTPUT_DIR/step6_CollectAllelicCounts/$NAME.log
}
export -f step6_CollectAllelicCounts

step7_ModelSegments (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step7_ModelSegments para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" ModelSegments \
     --denoised-copy-ratios $OUTPUT_DIR/step4_DenoiseReadCounts/${NAME}.denoisedCR.tsv  \
     --allelic-counts  $OUTPUT_DIR/step6_CollectAllelicCounts/${NAME}.allelicCounts.tsv \
     --output-prefix ${NAME} \
     -O $OUTPUT_DIR/step7_ModelSegments/  2> $OUTPUT_DIR/step7_ModelSegments/$NAME.log
}
export -f step7_ModelSegments

step8_CallCopyRatioSegments (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step8_CallCopyRatioSegments para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" CallCopyRatioSegments \
    --input $OUTPUT_DIR/step7_ModelSegments/${NAME}.cr.seg \
    --output $OUTPUT_DIR/step8_CallCopyRatioSegments/${NAME}.called.seg  2> $OUTPUT_DIR/step8_CallCopyRatioSegments/$NAME.log
}
export -f step8_CallCopyRatioSegments

step9_PlotModeledSegments (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step9_PlotModeledSegments para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  ${GATK} --java-options "-Xmx${MAXmem}G" PlotModeledSegments \
    --denoised-copy-ratios  $OUTPUT_DIR/step4_DenoiseReadCounts/${NAME}.denoisedCR.tsv  \
    --allelic-counts $OUTPUT_DIR/step7_ModelSegments/${NAME}.hets.tsv \
    --segments $OUTPUT_DIR/step7_ModelSegments/${NAME}.modelFinal.seg \
    --sequence-dictionary ${REF_FASTA}/Homo_sapiens_assembly38.dict \
    --minimum-contig-length 46709983 \
    --output $OUTPUT_DIR/step9_PlotModeledSegments/ \
    --output-prefix ${NAME}  2> $OUTPUT_DIR/step9_PlotModeledSegments/$NAME.log
}
export -f step9_PlotModeledSegments


step9_PlotModeledSegments_DOCKER (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step9_PlotModeledSegments_DOCKER para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

  docker run \
    --rm \
    -v $OUTPUT_DIR/:/gatk/my_data \
    -v ${REF_FASTA}:/gatk/hg38/ \
    -u $(id -u):$(id -g) broadinstitute/gatk \
    gatk --java-options "-Xmx${MAXmem}G" PlotModeledSegments \
    --denoised-copy-ratios  /gatk/my_data/step4_DenoiseReadCounts/${NAME}.denoisedCR.tsv  \
    --allelic-counts /gatk/my_data/step7_ModelSegments/${NAME}.hets.tsv \
    --segments /gatk/my_data/step7_ModelSegments/${NAME}.modelFinal.seg \
    --sequence-dictionary /gatk/hg38/Homo_sapiens_assembly38.dict \
    --minimum-contig-length 46709983 \
    --output /gatk/my_data/step9_PlotModeledSegments/ \
    --output-prefix ${NAME}  2>  $OUTPUT_DIR/step9_PlotModeledSegments/$NAME.log

  # docker run \
  # -v /home/scratch60/vlira_20nov2023/Result_Gatk-CNV.2023-11-07:/gatk/my_data \
  # -v /home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/:/gatk/hg38/ \
  # -u $(id -u):$(id -g) broadinstitute/gatk \
  # gatk PlotModeledSegments \
  # --denoised-copy-ratios /gatk/my_data/step4_DenoiseReadCounts/ROP-94-ExC85-xgenV2_S62.dedup.tags.bqsr.bam.denoisedCR.tsv \
  # --allelic-counts /gatk/my_data/step7_ModelSegments/ROP-94-ExC85-xgenV2_S62.dedup.tags.bqsr.bam.hets.tsv \
  # --segments /gatk/my_data/step7_ModelSegments/ROP-94-ExC85-xgenV2_S62.dedup.tags.bqsr.bam.modelFinal.seg \
  # --sequence-dictionary /gatk/hg38/Homo_sapiens_assembly38.dict \
  # --minimum-contig-length 46709983 \
  # --output /gatk/my_data/step9_PlotModeledSegments/ \
  # --output-prefix ROP-94-ExC85-xgenV2_S62.dedup.tags.bqsr.bam

  #docker run -v /home/scratch60/vlira_20nov2023/Result_Gatk-CNV.2023-11-07:/gatk/my_data -v /home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/:/gatk/hg38/ -it -u $(id -u):$(id -g) broadinstitute/gatk
}
export -f step9_PlotModeledSegments_DOCKER

step10_FilterCallCopyRatioSegments (){
  local SAMPLE=$1
  local NAME="${SAMPLE##*/}"
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Executando step10_FilterCallCopyRatioSegments para amostra: $NAME <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG
  
  COPY_RATIO=5
  CUTOFF_AMP=0.5
  CUTOFF_DEL=-1

  grep -v "^chr" $OUTPUT_DIR/step8_CallCopyRatioSegments/${NAME}.called.seg > $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/${NAME}.called.filt.seg
  # grep "^chr" $OUTPUT_DIR/step8_CallCopyRatioSegments/${NAME}.called.seg | awk -F "\t" '{ if(($5 <= ${CUTOFF_DEL} || $5 >= ${CUTOFF_AMP} ) && ($4 > ${COPY_RATIO})) print $_}' >> $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/${NAME}.called.filt.seg
  grep "^chr" $OUTPUT_DIR/step8_CallCopyRatioSegments/${NAME}.called.seg | awk -F "\t" '{ if(($5 <= -1 || $5 >= 0.5 ) && ($4 > 5)) print $0}' >> $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/${NAME}.called.filt.seg

  grep -v "^ROP" $OUTPUT_DIR/step8_CallCopyRatioSegments/${NAME}.called.igv.seg > $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/${NAME}.called.filt.igv.seg
  grep "^ROP" $OUTPUT_DIR/step8_CallCopyRatioSegments/${NAME}.called.igv.seg | awk -F "\t" '{ if(($7 <= -1 || $7 >= 0.5 ) && ($5 > 5)) print $0}' >> $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/${NAME}.called.filt.igv.seg

}
export -f step10_FilterCallCopyRatioSegments





echo "                           >>>>>> Starting Pipeline to Run Pipeline_GATK-CNV.XL.sh <<<<<<" >> $OUTPUT_LOG
date >> $OUTPUT_LOG

mkdir $OUTPUT_DIR/step1_PreprocessIntervals/
#step1_PreprocessIntervals

mkdir $OUTPUT_DIR/step2_AnnotateIntervals/
#step2_AnnotateIntervals

mkdir $OUTPUT_DIR/step3_CollectReadCounts/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step3_CollectReadCounts  "$@"' 'step3_CollectReadCounts'

mkdir $OUTPUT_DIR/step4_DenoiseReadCounts/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step4_DenoiseReadCounts  "$@"' 'step4_DenoiseReadCounts'

mkdir $OUTPUT_DIR/step5_PlotDenoisedCopyRatios/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step5_PlotDenoisedCopyRatios  "$@"' 'step5_PlotDenoisedCopyRatios'
xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step5_PlotDenoisedCopyRatios_DOCKER  "$@"' 'step5_PlotDenoisedCopyRatios_DOCKER'

mkdir $OUTPUT_DIR/step6_CollectAllelicCounts/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step6_CollectAllelicCounts  "$@"' 'step6_CollectAllelicCounts'

mkdir $OUTPUT_DIR/step7_ModelSegments/
#xargs -a $OUTPUT_DIR/TOY.samples.list -t -n1 -P${JOBS} bash -c 'step7_ModelSegments  "$@"' 'step7_ModelSegments'

mkdir $OUTPUT_DIR/step8_CallCopyRatioSegments/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step8_CallCopyRatioSegments  "$@"' 'step8_CallCopyRatioSegments'

mkdir $OUTPUT_DIR/step9_PlotModeledSegments/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step9_PlotModeledSegments  "$@"' 'step9_PlotModeledSegments'
xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step9_PlotModeledSegments_DOCKER  "$@"' 'step9_PlotModeledSegments_DOCKER'

mkdir $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/
#xargs -a $OUTPUT_DIR/samples.list -t -n1 -P${JOBS} bash -c 'step10_FilterCallCopyRatioSegments  "$@"' 'step10_FilterCallCopyRatioSegments'

#cat  $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/ROP-*.called.filt.igv.seg | grep "Sample"| head -1 > $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/ALL-ROP.called.filt.igv.seg
#cat  $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/ROP-*.called.filt.igv.seg | grep -v "Sample" >> $OUTPUT_DIR/step10_FilterCallCopyRatioSegments/ALL-ROP.called.filt.igv.seg

echo "" >> $OUTPUT_LOG
echo "                           >>>>>> End Pipeline <<< " >> $OUTPUT_LOG
date >> $OUTPUT_LOG
echo "" >> $OUTPUT_LOG
