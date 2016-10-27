#!/bin/sh

# ./completeRun.sh /srv/gsfs0/projects/salzman/Linda/DEBUG_MISEQ/trimmed complete /srv/gsfs0/projects/salzman/Linda/alignments NUGENpipeline 10 circReads 55

## parse parameters because we need some for each of the steps. 
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

# ntrim for denovo
if [ $# -ge 8 ]
then
  NTRIM=${8}
else
  NTRIM=50
fi

# add directory for downloaded bt1 indices & junction index .fa/bt index run by first pass
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

JUNCTION_MIDPOINT=${16}
## end parse parameters

# initial alignment
echo -e "\n............\nBeginning FIRST pass of findCircularRNA.sh.\n............\n"
echo "${SCRIPT_DIR}/findCircularRNA.sh \"$@\""
${SCRIPT_DIR}/findCircularRNA.sh "$@"  

# was having some I/O issues where output files weren't written even though findCircularRNA.sh was complete, but did get written a few minutes later
# so just taking a little pause here to make sure the files get output before moving on
sleep 300

# run GLM to output reports  
if [[ $MODE != *skipGLM* ]]
then
#  echo -e "\ncalling parseForAnalysis.sh"
#  echo "${SCRIPT_DIR}/analysis/parseForAnalysis.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE}"
#  ${SCRIPT_DIR}/analysis/parseForAnalysis.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} # get info we need from sam file for linear junctions and circular junctions
  echo -e "\ncalling predictJunctions.sh"
  echo "${SCRIPT_DIR}/analysis/predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME}"
  ${SCRIPT_DIR}/analysis/predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME} ${SCRIPT_DIR}
  cd ..
fi

# generate samples alignment statistics, also generates the unaligned reads

cd ${SCRIPT_DIR}/qualityStats
echo -e "\ncalling qualityStatsAll.sh"
echo "${SCRIPT_DIR}/qualityStats/qualityStatsAll.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${NTRIM} ${JUNCTION_DIR_SUFFIX}"
${SCRIPT_DIR}/qualityStats/qualityStatsAll.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${REPORTDIR_NAME} ${NTRIM} ${JUNCTION_DIR_SUFFIX}
cd ..

 
# run denovo analysis to generate new fasta file to use for unaligned run below
if [[ $MODE != *skipDenovo* ]]
then
  
  cd ${SCRIPT_DIR}/denovo_scripts
  echo -e "\ncalling process_directory_unaligned.pl"
  echo "perl ${SCRIPT_DIR}/denovo_scripts/process_directory_unaligned.pl ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${NTRIM} ${DENOVOCIRC} ${BT1_INDEX}"
  perl ${SCRIPT_DIR}/denovo_scripts/process_directory_unaligned.pl ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${NTRIM} ${DENOVOCIRC} ${BT1_INDEX}
  cd ..
  
  # run unaligned mode
  echo -e "\n............\nBeginning SECOND pass of findCircularRNA.sh.\n............\n"
  echo "${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_unaligned ${REPORTDIR_NAME} ${NTRIM} ${BT1_INDEX} ${BT2_INDEX} ${SCRIPT_DIR} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}"
  ${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_unaligned ${REPORTDIR_NAME} ${NTRIM} ${BT1_INDEX} ${BT2_INDEX} ${SCRIPT_DIR} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}
  
fi

# run pipeline in a read2-centric manner
#cd ${SCRIPT_DIR}/analysis
#  echo "./findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_unaligned ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}"
#  ./findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_unaligned ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}
  
# run R2 analysis to update junctional counts for R2
if [[ ${MODE} != *skipR2* ]]
then
  cd ${SCRIPT_DIR}/analysis
  # set up directories and swap alignment file names to run with R2 as R1
  echo "python analysis/createSwappedDirectories.py -a ${ALIGN_PARDIR} -d ${DATASET_NAME} -m swap -v"
  python createSwappedDirectories.py -a ${ALIGN_PARDIR} -d ${DATASET_NAME} -m swap -v 
  
  UNSWAPPED_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}
  UNSWAPPED_DATASET_NAME=${DATASET_NAME}
  DATASET_NAME=${DATASET_NAME}Swapped  # from here on out, the analysis mode will operate on the swapped directory
  SWAPPED_DIR=${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME}
  
  # if we ran previously in unaligned mode, we will want to analyze that output too. But no errors thrown if you didn't run in unaligned mode previously
  DENOVO_TASK_DATA_FILE=${DATASET_NAME}_unaligned.txt  # this is the to-be-created-below file
  f=`find ${ALIGN_PARDIR}/taskIdFiles -type f -name ${UNSWAPPED_DATASET_NAME}_unaligned.txt`
  if [ ! -f "$f" ]
  then
    NUM_DENOVO_FILES=0
  else
    # make a copy of the unaligned task data file with the correct name for the swapped run
    cp ${ALIGN_PARDIR}/taskIdFiles/${UNSWAPPED_DATASET_NAME}_unaligned.txt ${ALIGN_PARDIR}/taskIdFiles/${DENOVO_TASK_DATA_FILE}  
    NUM_DENOVO_FILES=`cat ${ALIGN_PARDIR}/taskIdFiles/${DENOVO_TASK_DATA_FILE} | wc -l`
  fi
  cd ..
  echo "NUM_DENOVO_FILES is ${NUM_DENOVO_FILES}"
  
  # runs analysis mode on the now-swapped reads
  echo "${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_R2analysis ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}"
  ${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_R2analysis ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}
  
  if [ ${NUM_DENOVO_FILES} -gt 0 ]
  then
    # run unaligned mode on the swapped files. Using string "unalign" instead of "unaligned" allows for selection of correct TASK_ID_FILE but prevents attempt to re-run denovo alignment
    echo -e "\n............\nBeginning THIRD pass of findCircularRNA.sh.\n............\n"
    echo "${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_R2analysis_unalign ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}"
    ${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME} ${OVERLAP} ${MODE}_R2analysis_unalign ${REPORTDIR_NAME} ${NTRIM} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}
  fi
  
  cd ${SCRIPT_DIR}/analysis
  echo "python combineSwappedReadsNaive.py -a ${UNSWAPPED_DIR} -b ${SWAPPED_DIR} -q ${READ_STYLE}"
  python combineSwappedReadsNaive.py -a ${UNSWAPPED_DIR} -b ${SWAPPED_DIR} -q ${READ_STYLE}

  if [[ ${MODE} != *skipGLM* ]]
  then
    echo "${SCRIPT_DIR}/predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME}"
    ${SCRIPT_DIR}/predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME} ${MODE} ${REPORTDIR_NAME}
    
    echo "python combineSwappedReadsGLM.py -a ${UNSWAPPED_DIR} -b ${SWAPPED_DIR} -q ${READ_STYLE}"
    python combineSwappedReadsGLM.py -a ${UNSWAPPED_DIR} -b ${SWAPPED_DIR} -q ${READ_STYLE}
  fi

  # swap back all of the alignment files to their original directories
  echo "python createSwappedDirectories.py -a ${ALIGN_PARDIR} -d ${UNSWAPPED_DATASET_NAME} -m restore -v"
  python createSwappedDirectories.py -a ${ALIGN_PARDIR} -d ${UNSWAPPED_DATASET_NAME} -m restore -v

fi

### Will delete folowwing code once confirmation of successful run ###

#echo -e "\npython createSwappedDirectories.py -a ${ALIGN_PARDIR} -d ${DATASET_NAME} -v"
#python createSwappedDirectories.py -a ${ALIGN_PARDIR} -d ${DATASET_NAME} -v
#
#run in analysis mode on new R2 directory
#echo -e "\n............\nBeginning THIRD pass of findCircularRNA.sh.\n............\n"
#echo "${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME}Swapped ${OVERLAP} ${MODE}_analysis ${REPORTDIR_NAME} ${NTRIM} ${BT1_INDEX} ${BT2_INDEX} ${SCRIPT_DIR} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}"
#${SCRIPT_DIR}/findCircularRNA.sh ${READ_DIR} ${READ_STYLE} ${ALIGN_PARDIR} ${DATASET_NAME}Swapped ${OVERLAP} ${MODE}_analysis ${REPORTDIR_NAME} ${NTRIM} ${BT1_INDEX} ${BT2_INDEX} ${SCRIPT_DIR} ${DENOVOCIRC} ${JUNCTION_DIR_SUFFIX} ${RD1_THRESH} ${RD2_THRESH} ${JUNCTION_MIDPOINT}
#
# run GLM to output reports for read2 pass
#echo -e "\ncalling parseForAnalysis.sh"
#echo "${SCRIPT_DIR}/analysis/parseForAnalysis.sh ${ALIGN_PARDIR} ${DATASET_NAME}Swapped ${MODE}_analysis"
#${SCRIPT_DIR}/analysis/parseForAnalysis.sh ${ALIGN_PARDIR} ${DATASET_NAME}Swapped ${MODE}_analysis # get info we need from sam file for linear junctions and circular junctions
#echo -e "\ncalling predictJunctions.sh"
#echo "${SCRIPT_DIR}/analysis/predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME}Swapped ${MODE}_analysis ${REPORTDIR_NAME} ${SCRIPT_DIR}"
#${SCRIPT_DIR}/analysis/predictJunctions.sh ${ALIGN_PARDIR} ${DATASET_NAME}Swapped ${MODE}_analysis ${REPORTDIR_NAME} ${SCRIPT_DIR}
#cd ..
#
#create combined report files
#cd ${SCRIPT_DIR}/analysis
#python combineSwappedReadsGLM.py -a ${ALIGN_PARDIR}/${DATASET_NAME}/${REPORTDIR_NAME} -b ${ALIGN_PARDIR}/${DATASET_NAME}Swapped/${REPORTDIR_NAME} -q ${READ_STYLE} -v
