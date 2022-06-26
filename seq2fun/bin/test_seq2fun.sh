#!/bin/bash
# Author: 	Hannes Reinwald
# Contact: 	hannes.reinwald@ime.fraunhofer.de

###############################
##  TEST Seq2Fun PARAMETERS  ##
###############################

## README: ----------------------------------------------------
# This script is designed to test your S2F search parameters. You can simply modify the parameters
# within this script and then check your output. For phylogenetic far distant organisms try to find
# a good compormise between accuracy and a good read mapping rate
# -------------------------------------------------------------

home=$(pwd)

### PARAMETERS ### --------------------------------------------

# Reference Database #
protDbIndex="/srv/attract/git/seq2fun_db/plants/plants_v2.0.fmi" #(--tfmi)
geneMap="/srv/attract/git/seq2fun_db/plants/plants_annotation_v2.0.txt"
# Download from: https://www.seq2fun.ca/database.xhtml
# Database version 3.0 will be soon released 

# S2F search params #
searchMode="tGREEDY" #tMEM or tGREEDY (default); Greedy mode is designed for organisms who do not have reference sequences in database.
MinScore="70" 	#--minscore; Specify minimum matching score (DEFAULT 80)
MaxMismatch="5" #--mismatch; Max Nbr of AA mismatches (DEFAULT 2)
maxNbase="5"	#--n_base_limit; max Nbr of N bases in nt seq (DEFAULT 5)
ntMinLen="40"	#--length_required; min nt seq length to be used for translation. Must be < than readL (DEFAULT 60)
AAminLen="10"	#--minlength; min length of AA seq length (DEFAULT 19 (tGreedy); 13 (tMEM)) - Must be < than $readL/3 (if front clipped: < ($readL-6)/3)
AAmaxLen="60" 	#--maxtranslength; Maximum cutoff length of translated peptides (DEFAULT 60)

#For PE only
OverlLen="22" 	#--overlap_len_require; only for paired-end, min overlap between read1 & read2

# For details on the parameters check the S2F manual under:
# https://www.seq2fun.ca/manual.xhtml 
#If your studied species has no or limited number of close-related reference species in the 
#database. You can tune parameters to obtain a better balance between sensitivity and precision.
#For example: 
#1) increasing the number of mismatches (--mismatch) from default 2 to 3 or 4; 
#2) decreasing the minimum matching length (--minlength) from default 19; 
#3) decreasing the minimum BLOSUM62 score (--minscore) from default 80;
#4) decreasing the maximum length cutoff of the translated amino acid sequences (--maxtranslength)
#for overlapped paired-end reads from default 60; will increase the mapping reads or the mapping chance for the highly divergent homologs.
# -----------------------------------------------------------------

### In-script settings ###
CPUs=$[$(nproc)-2]
S2F=/srv/attract/git/Seq2Fun/bin/seq2fun #path to seq2fun
# Check if provided files exits:
[[ ! -f $protDbIndex ]] && printf "ERROR: \t$protDbIndex \ndoes not exist. Please check the provided path! \n" && exit 0
[[ ! -f $geneMap ]] && printf "ERROR: \t$geneMap \ndoes not exist. Please check the provided path! \n" && exit 0


### run Seq2Fun ###
mkdir $home/s2f_testRun #output dir

printf "##\n`date`\nRunning S2F with the following search settings:\n
ProtDb Index:\t(--tfmi)\t\t\t $protDbIndex
Gene Map:\t\t(--genemap)\t\t\t $geneMap
Search Mode:\t(--mode)\t\t\t $searchMode
Min MatchScore:\t(--minscore)\t\t $MinScore
Max Mismatches:\t(--mismatch)\t\t $MaxMismatch
MaxN Bases:\t\t(--n_base_limit)\t $maxNbase
MinNt Length:\t(--length_required)\t $ntMinLen
MinAA Length:\t(--minlength)\t\t $AAminLen
MaxAA Length:\t(--maxtranslength)\t $AAmaxLen \n##\n\n" 2>&1 | \
	tee -a $home/s2f_testRun/S2F_testRun.config


INPUT=$(echo $(find -name *.fastq.gz) | awk '{print $1;}') # finds all fastq files in pwd and takes the first a input for test run
out=$(echo $INPUT | sed 's,^.*/,,')

$S2F -i $INPUT --prefix $home/s2f_testRun/${out/.fastq.gz/} --thread $CPUs --profiling \
	--mode $searchMode --tfmi $protDbIndex --genemap $geneMap --length_required $ntMinLen \
	--minscore $MinScore --minlength $AAminLen --mismatch $MaxMismatch \
	--maxtranslength $AAmaxLen -V #\
	#--allFragments \
	#--trim_front1 6 --trim_front2 6 \
	#--outputMappedCleanReads --outputReadsAnnoMap #Could be used for de-novo gene assembly e.g. Trinity... 

### END OF SCRIPT ###