#!/bin/bash
# Author: 	Hannes Reinwald
# Contact: 	hannes.reinwald@ime.fraunhofer.de

########################
##  Seq2Fun Analysis  ##
########################

## !!! PROVIDE A CONFIG FILE FOR THE ANALYSIS !!! ##
#configFile=/some/path/to/s2f.config
#configFile=$(find $(dirname $(realpath $0)) -name s2f4plants_SR_50bp.config)
configFile=$(ls *[.]config)
# --------------------------------------------------


# Check if config file exists before running the script
[[ ! -f $configFile ]] && printf "\nERROR: \t$configFile does not exist.
A config file must be provided for the S2F run.
Please check the provided path! \n" && exit 0

# Load config file
source $configFile

### Compute suitable --minlength & --n_base_limit from readL parameter ###
#-eq # Equal
#-ne # Not equal
#-lt # Less than
#-le # Less than or equal
#-gt # Greater than
#-ge # Greater than or equal
if [[ $ntMinLen -ge $readL ]]; then (( ntMinLen=(readL-6) )); fi

# --minlength	minimum matching length of amino acid sequence in comparison with protein database
#				default value 19, for GREEDY and 13 for MEM model
# !!! This default parameter is way to strict for SE short read lengths (e.g. 50 bp)
# => Needs to be adjusted accordingly!
# e.g. 50bp reads => ~16 AA length at best!; Trimming the first 6 reads => 14 AA at best
if [[ $searchMode == "tGREEDY" ]]; then 
	N=19
else
	N=13
fi
# get floor value
(( z=($readL-6)/3-5 )) && \
## AAminLen: for short reads (e.g. 50 bp) values will be computed. Otherwise default
if [[ $z -ge $N ]]; then
	AAminLen=$N
else
	AAminLen=$z
fi

## maxNbase: Max 10% of single read
perc="10" #10%; already mulitplied by 100
#maxNbase="5" #DEFAULT
(( z=(($readL*100)*perc)/(100*100) )) && \
if [[ $z < $maxNbase ]]; then
	maxNbase=5
else
	maxNbase=$z
fi &&\
# If maxNbase > 10 set back to 10.
if [[ $maxNbase -ge 10 ]]; then maxNbase=10; fi
#echo $maxNbase



### RUN Rscript to build sampleTable file and output dir ###
printf "\n`date` - Starting Seq2Fun (S2F) analysis pipeline ...
Creating sample table for batch processing from provided raw read files ...\n\n"
# both scripts (this and the Rscript should be loacted in the same bin folder)
rPATH=$(dirname $(realpath $0))
Rscript $rPATH/sampleTable4S2F_fastp.R


### In-script settings ###
CPUs=$[$(nproc)-2]
home=$(pwd)
sampleTable=$(ls *_sampleTable.tab)
S2F=/srv/attract/git/Seq2Fun/bin/seq2fun #path to seq2fun

# Check if provided files exits:
[[ ! -f $protDbIndex ]] && printf "ERROR: \t$protDbIndex \ndoes not exist. Please check the provided path! \n" && exit 0
[[ ! -f $geneMap ]] && printf "ERROR: \t$geneMap \ndoes not exist. Please check the provided path! \n" && exit 0
[[ ! -f $sampleTable ]] && printf "ERROR: \t$sampleTable \ndoes not exist. Please check the provided path! \n" && exit 0


### run Seq2Fun ###
printf "##\n`date`\nRunning S2F with the following search settings:\n
ProtDb Index:\t(--tfmi)\t\t\t $protDbIndex
Gene Map:\t\t(--genemap)\t\t\t $geneMap
Sample Table:\t(--sampletable)\t\t $sampleTable
Search Mode:\t(--mode)\t\t\t $searchMode
Min MatchScore:\t(--minscore)\t\t $MinScore
Max Mismatches:\t(--mismatch)\t\t $MaxMismatch
MaxN Bases:\t\t(--n_base_limit)\t $maxNbase
MinNt Length:\t(--length_required)\t $ntMinLen
MinAA Length:\t(--minlength)\t\t $AAminLen
MaxAA Length:\t(--maxtranslength)\t $AAmaxLen \n##\n\n" 2>&1 | \
	tee -a $home/seq2fun/S2F_run.config


if [[ $TrimReads == "TRUE" ]]; then
## Trimming random hexamer ends ##
	if [[ $readType == "SR" ]]; then
	### For SR ### 
	$S2F --thread $CPUs --profiling --n_base_limit $maxNbase --trim_front1 6 \
	--sampletable $sampleTable --mode $searchMode \
	--tfmi $protDbIndex --genemap $geneMap --length_required $ntMinLen \
	--minscore $MinScore --minlength $AAminLen --mismatch $MaxMismatch \
	--maxtranslength $AAmaxLen -V #\
	#--allFragments #this parameter will fore s2f to use all AAseq generated for mapping which will enhance the number of mapped reads but also false-positives
	#--outputMappedCleanReads --outputReadsAnnoMap #Could be used for de-novo gene assembly e.g. Trinity... 
	else
	### For PE ###
	$S2F --thread $CPUs --profiling --n_base_limit $maxNbase --trim_front1 6 --trim_front2 6 \
	--sampletable $sampleTable --mode $searchMode \
	--tfmi $protDbIndex --genemap $geneMap --length_required $ntMinLen \
	--minscore $MinScore --minlength $AAminLen --mismatch $MaxMismatch \
	--maxtranslength $AAmaxLen --overlap_len_require $OverlLen -V #\
	#--outputMappedCleanReads --outputReadsAnnoMap
	fi
else
## No trimming of reads ##
	if [[ $readType == "SR" ]]; then
	### For SR ### 
	$S2F --thread $CPUs --profiling --n_base_limit $maxNbase \
	--sampletable $sampleTable --mode $searchMode \
	--tfmi $protDbIndex --genemap $geneMap --length_required $ntMinLen \
	--minscore $MinScore --minlength $AAminLen --mismatch $MaxMismatch \
	--maxtranslength $AAmaxLen -V #\
	#--outputMappedCleanReads --outputReadsAnnoMap 
	else
	### For PE ###
	$S2F --thread $CPUs --profiling --n_base_limit $maxNbase \
	--sampletable $sampleTable --mode $searchMode \
	--tfmi $protDbIndex --genemap $geneMap --length_required $ntMinLen \
	--minscore $MinScore --minlength $AAminLen --mismatch $MaxMismatch \
	--maxtranslength $AAmaxLen --overlap_len_require $OverlLen -V #\
	#--outputMappedCleanReads --outputReadsAnnoMap 
	fi
fi


### organize the ouput ###
mkdir $home/seq2fun/s2f_reports
mkdir $home/seq2fun/s2fID_abundance
cd $home/seq2fun && \
mv -f $(ls *_report.json *_report.html) $home/seq2fun/s2f_reports/ && \
mv -f $(ls *_s2fid_abundance.txt) $home/seq2fun/s2fID_abundance/ && \
tar -czvf s2fID_abundance.tar.gz s2fID_abundance && rm -fr s2fID_abundance && \
tar -czvf s2f_reports.tar.gz s2f_reports && rm -fr s2f_reports && \
cd $home

### END OF SCRIPT ###