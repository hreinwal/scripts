#!/bin/bash
# Author: 	Hannes Reinwald
# Contact: 	hannes.reinwald@ime.fraunhofer.de

################################
##  Seq2Fun Analysis Pipeline ##
################################

## README: ------------------------------------------------
# This script will:
# 1. Run some basic NGS quality filtering and QC (fastp --> FastQScreen, FastQC => MultiQC)
# 2. Run seq2fun on fastp filtered reads (!!! Provide Config file !!!)
# 3. Run DGEA (DESeq2) and ORA/GSEA (clusterProfiler) on the mapped result tables

# INPUT:
# - coldata.csv file (sample annotation file)
# - raw RNASeq reads in the following format:

# For SR:
# /raw_reads/
#	sampleA.fastq.gz
#	sampleB.fastq.gz

# For PE:
# /raw_reads/
#	/sampleA/
#		sampleA_R1.fastq.gz
#		sampleA_R2.fastq.gz
#	/sampleB/
#		sampleB_R1.fastq.gz
#		sampleB_R2.fastq.gz
# ---------------------------------------------------------


####################################################
## !!! PROVIDE A CONFIG FILE FOR THE ANALYSIS !!! ##
####################################################

#configFile=/some/path/to/s2f.config
configFile=$(find $(dirname $(realpath $0)) -name s2f4plants_SR_50bp.config)

####################################################


# Check if proved config file exists:
if [ -f "$configFile" ]; then
	echo "$configFile exists" 2>&1 | tee -a runtime.log
	cp -f $configFile ./run.config #copy config file to PDIR
else
	echo "$configFile does not exists! 
Please check the provided path!
Exiting script"
	sleep 5s && exit 0
fi

# Run params #
source $configFile
PDIR=$(pwd) # Project directory
ScriptPATH="`dirname $(realpath $0)`/bin" # Specify path to script repo


### Pre-RUN Check ### -----------------------------------------------------------
echo "Executing: `realpath $0`" >> runtime.log
# Check if project dir contains 'raw_reads' folder
[[ ! -d raw_reads ]] && \
	printf "ERROR: \nWorking directory does not contain a 'raw_reads' directory! \nPlease check your present working dirteory (pwd) \n" && \
	sleep 10s && exit 0

# Runtime log 
SECONDS=0
printf "`date` - START Seq2Fun Analysis pipeline \n" >> runtime.log


### fastp filtering and data QC (FastQC, FastQ-Screen) ### ---------------------
printf "\n#### fastp, fastQscreen, fastQC ####" 2>&1 | tee -a runtime.log
echo "
 `date` --- START fastp and sequence QC ---
" 2>&1 | tee -a run.stout

if [[ $readType == "SR" ]]; then
	# if SR --> run fastp for SR
	{ time bash $ScriptPATH/fastp_seqQC_SR.sh 2>&1 | tee -a run.stout ; } 2>> runtime.log
elif [[ $readType == "PE" ]]; then
	# if PE --> run fastp for PE
	{ time bash $ScriptPATH/fastp_seqQC_PE.sh 2>&1 | tee -a run.stout ; } 2>> runtime.log
else
	echo "readType parameter not defined in config file! Must be one of 'SR' or 'PE'!\n" 2>&1 | tee -a runtime.log
fi

### Run multiQC ### ------------------------------------------------------------
printf "\n#### multiQC ####" 2>&1 | tee -a runtime.log
{ time bash $ScriptPATH/multiQC.sh 2>&1 | tee -a run.stout ; } 2>> runtime.log
echo "
 `date` --- END fastp and sequence QC ---
" 2>&1 | tee -a run.stout


### Seq2Fun ### ----------------------------------------------------------------
printf "\n#### Seq2Fun ####" 2>&1 | tee -a runtime.log
{ time bash $ScriptPATH/seq2fun.sh 2>&1 | tee -a run.stout ; } 2>> runtime.log



### DGEA (DESeq2) & GSEA/ORA ### -----------------------------------------------
printf "\n#### DGEA 4 S2F output in R ####" 2>&1 | tee -a runtime.log
source activate Renv
{ time Rscript $ScriptPATH/DESeq2_pairwise.4S2F.R --verbose 2>&1 | tee -a run.stout ; } 2>> runtime.log




### CLEAN UP ### ---------------------------------------------------------------
rm -f $PDIR/run.config $PDIR/run.stout #Rmv this line to keep stout in case for the need of bug search
# rm fastp.fastq files!!!
cd $PDIR/seqQC/fastp && rm -fr $(find -name "*[.]fastq[.]gz") && cd $PDIR

## END OF SCRIPT ##