#!/bin/bash

# Author:  Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# -------------------------------------------------------------
# To run this script go to your project folder containing a 'raw_reads' folder with
# fastq files and execute this script via:
# bash /some/path/to/NAME_OF_THIS_SCRIPT.sh
# -------------------------------------------------------------

############################
###  RNASeq CONFIG FILE  ###
############################

# !!!! This file must be provided by the user !!!!
#CONF=provide/config/file/here
CONF=/srv/attract/Scripts/RNAseq_custompipes/config/drerio_50bp.conf #lminor_50bp.conf

# Check if proved config file exists:
if [ -f "$CONF" ]; then
	echo "$CONF exists" 2>&1 | tee -a runtime.log
	cp -f $CONF ./config #copy config file to PDIR
else
	echo "$CONF does not exists! 
Please check the provided path!
Exiting script"
	exit 0
fi
##########################


PDIR=$(pwd) # Project directory
echo "Executing: `realpath $0`" >> runtime.log
# Check if project dir contains 'raw_reads' folder
[[ ! -d raw_reads ]] && \
	printf "ERROR: \nWorking directory does not contain a 'raw_reads' directory! \nPlease check your present working dirteory (pwd) \n" && \
	sleep 15s && \
	exit 0
# Specify path to script repo
ScriptPATH="`dirname $(realpath $0)`/bin" #more advanced and more flexible
# Runtime log 
SECONDS=0
printf "`date` - START RNAseq pipeline \n" >> runtime.log



### fastp & sequence QC ### ---------------------------------------------------
printf "\n#### fastp, fastQscreen, fastQC ####" 2>&1 | tee -a runtime.log
echo "
 `date` --- START fastp and sequence QC ---
" 2>&1 | tee -a run.stout
{ time bash $ScriptPATH/fastp_seqQC_SR.sh 2>&1 | tee -a run.stout ; } 2>> runtime.log
echo "
 `date` --- END fastp and sequence QC ---
" 2>&1 | tee -a run.stout



### STAR with fastp filtered reads & gene count ### ----------------------------
printf "\n#### STAR mapping ####" 2>&1 | tee -a runtime.log
echo "
 `date` --- START STAR mapping, gene counts ---
" 2>&1 | tee -a run.stout
{ time bash $ScriptPATH/STAR_SR.sh 2>&1 | tee -a run.stout ; } 2>> runtime.log
echo "
 `date` --- END STAR mapping, gene counts ---
" 2>&1 | tee -a run.stout



### Merging GeneCount files to single CountMatrix via R --------------------------
printf "\n#### CountMatrix generation (R) ####" 2>&1 | tee -a runtime.log | tee -a run.stout
printf "\n" >> run.stout
{ time Rscript $ScriptPATH/CountMatrix_generator.R --verbose 2>&1 | tee -a run.stout ; } 2>> runtime.log



### Alignment QC (Qualimap, samtools, RSeQC) ### ------------------------------
printf "\n#### AlignmentQC - RSeQC Qualimap samtools ####" 2>&1 | tee -a runtime.log
echo "
 `date` --- START RSeQC Qualimap samtools ---
" 2>&1 | tee -a run.stout
{ time bash $ScriptPATH/alignmentQC_SR.sh 2>&1 | tee -a run.stout ; } 2>> runtime.log
echo "
 `date` --- END RSeQC Qualimap samtools ---
" 2>&1 | tee -a run.stout


### Final cleanup ### ----------------------------------------------------------
cd $PDIR/seqQC/fastp && rm -fv $(ls *.fastq.gz) && \
	cd $PDIR/STAR && rm -fv $(ls *.out.bam)
cd $PDIR

# Analysis finished
echo "
 `date` - FINISHED SUCCESFULLY :) Yay!" 2>&1 | tee -a run.stout

# Runtime log
duration=$SECONDS
echo "
`date` - END RNAseq pipeline
$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed." 2>&1 | tee -a runtime.log

# SET PERMISSION RIGHTS
chmod -R g+rw *

### END OF SCRIPT ###