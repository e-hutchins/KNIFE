#!/bin/sh

# ./findCircularRNA.sh /srv/gsfs0/projects/salzman/Linda/DEBUG_MISEQ/trimmed complete /srv/gsfs0/projects/salzman/Linda/alignments NUGENpipeline 10
#Read directory
#readStyle (complete or appended)
#Alignment parent directory
#dataset name
#overlap
#Do analysis only, Sam or bam or sorted bam (sam, bam, bam_sort, analysis, sam_large, bam_large, bam_sort_large). Default sam
#report directory name. Default circReads
#ntrim
#bt1 directory
#bt2 directory
#script directory
#Rd1 threshold. 
#Rd2 threshold. 
#text appended to junction id output directory. Default none
#junction midpoint. Default 150

# some error checking
if [ $# -eq 14 ]
then
  echo "Read 1 threshold was specified without Read 2 threshold. Please include Read 2 threshold or check your arguments if you did not intend to specify Read 1 threshold ${14}"
  exit 2
fi

READ_DIR=${1}
READ_STYLE=${2}
ALIGN_PARDIR=${3}
DATASET_NAME=${4}
OVERLAP=${5}
if [ $# -ge 6 ]
then
  MODE=${6}
else
  MODE=sam
fi

if [ $# -ge 7 ]
then
  REPORTDIR_NAME=${7}
else
  REPORTDIR_NAME=circReads
fi

# ntrim for denovo (this isn't actually used here, but to be consistent with how findCircularRNA.sh is called on sherlock I am adding it in)
if [ $# -ge 8 ]
then
  NTRIM=${8}
else
  NTRIM=50
fi

# add directory for downloaded bt1 indices. this is only used in completeRun.sh script, just keeping same order of positional arguments.
BT1_INDEX=${9}

# add directory for downloaded bt2 indices
BT2_INDEX=${10}

#script directory
SCRIPT_DIR=${11}

# should denovo contains only circles (1 means only circles)
if [ $# -ge 12 ]
then
  DENOVOCIRC=${12}
else
  DENOVOCIRC=1
fi

JUNCTION_DIR_SUFFIX=${13}
RD1_THRESH=${14}
RD2_THRESH=${15}


if [ $# -ge 16 ]
then
  JUNCTION_MIDPOINT=${16}
else
  JUNCTION_MIDPOINT=150
fi

# a few special flags for the unaligned mode
if [[ $MODE = *unaligned* ]]
then
  echo "...Starting findCircularRNA.sh in unaligned mode..."
  TASK_FILE_READ_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/orig/unaligned
  UFLAG="-u"
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}_unaligned.txt
else
  TASK_FILE_READ_DIR=${READ_DIR}
  TASK_DATA_FILE=${ALIGN_PARDIR}/taskIdFiles/${DATASET_NAME}.txt
fi

# set up info about this dataset, create directory structures
if [[ $MODE != analysis* ]]
then
  # set up
  echo -e "\ncalling writeTaskIdFiles.py"
  python ${SCRIPT_DIR}/analysis/writeTaskIdFiles.py -r ${TASK_FILE_READ_DIR} -a ${ALIGN_PARDIR} -d ${DATASET_NAME} ${UFLAG}
  
  # select correct prefix name to use for bowtie index files
  if [[ $MODE = *mouse* ]]
  then
    bt_prefix="mm10"
  elif [[ $MODE = *rat* ]]
  then
    bt_prefix="rn5"
  elif [[ $MODE = *fly* ]]
  then
    bt_prefix="dm3"
  elif [[ $MODE = *pombe* ]]
  then
    bt_prefix="ASM294v2_23"
  elif [[ $MODE = *crypto* ]]
  then
    bt_prefix="cryptococcus_neoformans_grubii_h99"
  elif [[ $MODE = *cerevisiae* ]]
  then
    bt_prefix="Scer"
  elif [[ $MODE = *mikatae* ]]
  then
    bt_prefix="Smik"
  elif [[ $MODE = *bayanus* ]]
  then
    bt_prefix="Sban"
  elif [[ $MODE = *HSV* ]]
  then
    bt_prefix="KOS"
  elif [[ $MODE = *capsas* ]]
  then
    bt_prefix="capsaspora_atcc_30864_2"
  elif [[ $MODE = *rosetta* ]]
  then
    bt_prefix="salpingoeca_rosetta_1"
  else
    bt_prefix="hg19"
  fi
fi

NUM_FILES=`cat $TASK_DATA_FILE | wc -l`

# alignments
if [[ $MODE != *analysis* ]]
then
  # have to be inside the index directory so bowtie can find the indices
  if [[ $MODE = *unaligned* ]]
  then
    cd ${ALIGN_PARDIR}/${DATASET_NAME}/denovo
  else
    cd ${BT2_INDEX}
  fi
  
  for (( i=1; i<=NUM_FILES; i++ ))
  do
    READ_FILE=`awk 'FNR == '${i}' {print $1}' $TASK_DATA_FILE`
    SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
    echo "MODE is $MODE"
    if [[ $MODE = *unaligned* ]]
    then
      echo "starting junction alignment"
      ${SCRIPT_DIR}/analysis/align.sh ${READ_FILE} ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} denovo_${DATASET_NAME}_${DENOVOCIRC} denovo denovo_${DATASET_NAME}_onlycircles${DENOVOCIRC}.fa
      echo -e "junction alignment complete "`date`
      #echo -e "\nLaunched align into the background "`date`
    else
      echo "starting genome, ribosome, and transcriptome alignment"
      ${SCRIPT_DIR}/analysis/align.sh ${READ_FILE} ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${bt_prefix}_genome genome ${bt_prefix}_genome.fa
      echo "genome alignment complete "`date`
      ${SCRIPT_DIR}/analysis/align.sh ${READ_FILE} ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${bt_prefix}_ribosomal ribo ${bt_prefix}_ribosomal.fa
      echo "ribosomal alignment complete "`date`
      ${SCRIPT_DIR}/analysis/align.sh ${READ_FILE} ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${bt_prefix}_transcriptome transcriptome ${bt_prefix}_transcriptome.fa
      echo "transcriptome alignment complete "`date`
      #echo -e "\nLaunched 3 align.sh's into the background "`date`
    fi
    run_count=`ps -ealf | grep align.sh | grep ${USER} | grep -v grep | wc -l`
    while [ "$run_count" -gt 3 ]
    do
      sleep 1
      run_count=`ps -ealf | grep align.sh | grep ${USER} | grep -v grep | wc -l`
    done
  done
    
  # Wait until finished
  echo 'Done looping over bowtie, awaiting the last few align.sh to finish '`date`
  while [ "$run_count" -gt 0 ]
  do
     sleep 1
     run_count=`ps -ealf | grep align.sh | grep ${USER} | grep -v grep | wc -l`
  done
  echo 'All non-junction align.sh are done '`date`

  # now do alignments for the junction indices, calling 1 at a time but using multiple processors to speed it up
  if [[ $MODE != *unaligned* ]]
  then
    echo -e "\nstarting junction align"
    for (( i=1; i<=NUM_FILES; i++ ))
    do
      READ_FILE=`awk 'FNR == '${i}' {print $1}' $TASK_DATA_FILE`
      SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
      ${SCRIPT_DIR}/analysis/align.sh ${READ_FILE} ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${bt_prefix}_junctions_scrambled junction ${bt_prefix}_junctions_scrambled.fa
      echo "scrambled junction alignment complete "`date`
      ${SCRIPT_DIR}/analysis/align.sh ${READ_FILE} ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${bt_prefix}_junctions_reg reg ${bt_prefix}_junctions_reg.fa
      echo "regular junction alignment complete "`date`
      #echo "Launched 2 junction aligns into the background "`date`
      
      # only want to launch for 1 sample at a time and then wait for completion so we don't overwhelm the server
      run_count=`ps -ealf | grep align.sh | grep ${USER} | grep -v grep | wc -l`
      while [ "$run_count" -gt 0 ]
      do
        sleep 1
        run_count=`ps -ealf | grep align.sh | grep ${USER} | grep -v grep | wc -l`
      done
   done
    echo 'All junction align.sh are done '`date`
  fi
  # then change back so relative paths below work  
  cd ..
fi

#echo -e "\nready to preprocess"
# preprocessing
for (( i=1; i<=NUM_FILES; i++ ))
do
  SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
  ${SCRIPT_DIR}/analysis/preprocessAlignedReads.sh ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${JUNCTION_MIDPOINT} ${OVERLAP} ${JUNCTION_DIR_SUFFIX} &
  echo "Launched preprocessAlignedReads.sh into the background "`date`
  run_count=`ps -ealf | grep preprocessAlignedReads.sh | grep ${USER} | grep -v grep | wc -l`
  while [ "$run_count" -gt 3 ]
  do
    sleep 1
    run_count=`ps -ealf | grep preprocessAlignedReads.sh | grep ${USER} | grep -v grep | wc -l`
  done
done

# Wait until finished
echo 'Done looping over preprocessing, awaiting the last few preprocessAlignedReads.sh to finish '`date`
while [ "$run_count" -gt 0 ]
do
  sleep 1
  run_count=`ps -ealf | grep preprocessAlignedReads.sh | grep ${USER} | grep -v grep | wc -l`
done

# was having some I/O issues where output files weren't written even though preprocessAlignedReads.sh was complete, but did get written a few minutes later
# so just taking a little pause here to make sure the files get output before moving on
sleep 300  

# analysis
for (( i=1; i<=NUM_FILES; i++ ))
do
  SAMPLE_ID=`awk 'FNR == '${i}' {print $2}' $TASK_DATA_FILE`
  ${SCRIPT_DIR}/analysis/filterFDR.sh ${MODE} ${SAMPLE_ID} ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${READ_STYLE} ${OVERLAP} ${SCRIPT_DIR} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} &
  echo -e "\nLaunched filterFDR.sh into the background "`date`
  run_count=`ps -ealf | grep filterFDR.sh | grep ${USER} | grep -v grep | wc -l`
  while [ "$run_count" -gt 3 ]
  do
    sleep 1
    run_count=`ps -ealf | grep filterFDR.sh | grep ${USER} | grep -v grep | wc -l`
  done
done

# Wait until finished
echo 'Done looping over filterFDR, awaiting the last few filterFDR.sh to finish '`date`
while [ "$run_count" -gt 0 ]
do
  sleep 1
  run_count=`ps -ealf | grep filterFDR.sh | grep ${USER} | grep -v grep | wc -l`
done

echo -e '\n\nCompleted 1 call to findCircularRNA.sh '`date`