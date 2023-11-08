WD="/home/venus/mar/vlira"

CRAM_DIR="$WD/samples/"
CRAM_FILES=$(find "$CRAM_DIR" -maxdepth 1 -mindepth 1  -name '*.cram')
REF_FASTA="$WD/reference/"
GATK="$WD/gatk-4.3.0.0/./gatk"
TARGET="$WD/reference/xgen-exome-research-panel-v2-targets-hg38.autossome.bed"

OUTPUT_DIR="RESULT_PON-GATK_CNV-gatk4.3.0.0-2023-11-06"
mkdir $OUTPUT_DIR
mkdir $OUTPUT_DIR/PreprocessAnnotateIntervals/
mkdir $OUTPUT_DIR/CollectReadCounts/
mkdir $OUTPUT_DIR/PanelOfNormals/

TIME_FILE="$OUTPUT_DIR/$OUTPUT_DIR.log"




echo "                                                     >>>>>> Starting Pipeline  to create PON CNV <<<<<<" >> $TIME_FILE
date >> $TIME_FILE

STAGE_PreprocessIntervals(){
  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_PreprocessIntervals <<<<<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK PreprocessIntervals \
        -L $TARGET \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $OUTPUT_DIR/PreprocessAnnotateIntervals/preprocessed_intervals.interval_list  2>> $TIME_FILE
}

STAGE_AnnotateIntervals() {
  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_AnnotateIntervals <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK AnnotateIntervals \
      -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
      -L $OUTPUT_DIR/PreprocessAnnotateIntervals/preprocessed_intervals.interval_list  \
      --interval-merging-rule OVERLAPPING_ONLY \
      -O $OUTPUT_DIR/PreprocessAnnotateIntervals/annotated_intervals.tsv 2>> $TIME_FILE
}

STAGE_CollectReadCounts(){
  local SAMPLE=$1
  NAME="${SAMPLE##*/}"
  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_CollectReadCounts <<<" >> $TIME_FILE
  echo ">>>>>> Executando para Amostra: "$NAME" <<<" >> $TIME_FILE
  date >> $TIME_FILE


  $GATK CollectReadCounts \
        -I $SAMPLE \
        -L $OUTPUT_DIR/PreprocessAnnotateIntervals/preprocessed_intervals.interval_list  \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta  \
        --format TSV \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O $OUTPUT_DIR/CollectReadCounts/$NAME.counts.tsv 2> $OUTPUT_DIR/CollectReadCounts/$NAME.log
}


STAGE_CreateReadCountPanelOfNormals(){
 local SAMPLE_READCOUNT=$(find "$OUTPUT_DIR"/CollectReadCounts/ -maxdepth 1 -mindepth 1  -name '*.tsv')
 local SAMPLE_HDF5=$(echo $SAMPLE_READCOUNT| sed 's/\s/ -I  /g')

  echo "" >> $TIME_FILE
  echo ">>>>>> STAGE_CreateReadCountPanelOfNormals 100-eigensamples.hdf5 <<<" >> $TIME_FILE
  date >> $TIME_FILE

  $GATK  CreateReadCountPanelOfNormals \
    -I $SAMPLE_HDF5 \
    --annotated-intervals $OUTPUT_DIR/PreprocessAnnotateIntervals/annotated_intervals.tsv \
    --number-of-eigensamples 100 \
    -O $OUTPUT_DIR/PanelOfNormals/PON.100COVID.100-eigensamples.2023-11-06.hdf5 2> $OUTPUT_DIR/PanelOfNormals/PON.100COVID.100-eigensamples.2023-11-06.log
}

## STAGE_0
   STAGE_PreprocessIntervals
   STAGE_AnnotateIntervals

# # #STAGE_1
 for run in $CRAM_FILES; do 
  STAGE_CollectReadCounts "$run"  >> $TIME_FILE
 done

STAGE_CreateReadCountPanelOfNormals



echo "" >> $TIME_FILE
echo ">>>>>> End Pipeline <<< " >> $TIME_FILE
date >> $TIME_FILE
echo "" >> $TIME_FILE
