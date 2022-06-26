#!/bin/bash

# Author:  Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# Navigate to your project folder (/srv/attract/seq_files/project_folder) then execute this script.
# This script will:
# - map raw reads via STAR and will generate a mapped reads output file in BAM format
# - generate a Count List (ReadsPerGene.out.tab) similar to htseq count with default parameters
#       col1: gene ID   col2: counts of unstranded RNAseq col3: counts for 1st read col4: counts for 2nd read

#############################################################################################

home=$(pwd)
source $home/config
source activate rnaseq

##############
##   STAR   ##
##############

## STAR parameters 
CPUs=$[$(nproc)-2] # 14
OUT=$outputFormat
TYPE=$outputType #Unsorted (Sorting requiered for Samtools)
MAXmultimap=$outFilterMultimapNmax  #Max number of multiple alignments allowed for a single read (Set to 1 for uniquely aligned reads); DEFAULT = 20
BAMsortRAM=$limitBAMsortRAM  # > 10Gb!
seedWindow=$seedPerWindowNmax
MAXmismatch=$outFilterMismatchNoverReadLmax

echo "
 `date` - Start STAR mapping with the following settings:
    - MaxMultimap  : $MAXmultimap
    - Output format: $OUT $TYPE
    - nThreads     : $CPUs
    - BAM sort RAM : $[BAMsortRAM/1000000000] GB
    - Ref. Genome  : $STARindexRefGenome
"

# Run STAR mapping with fastp filtered reads
cd $home/seqQC/fastp #INPUT fastq files
mkdir $home/STAR     #OUTPUT for STAR alignment

for SID in $(ls -d */); do
    cd $SID
    echo "
    Processing paired-end sample files in $SID"
	STAR --runThreadN $CPUs \
    --seedPerWindowNmax $seedWindow \
    --runMode alignReads \
    --genomeDir $STARindexRefGenome \
    --genomeLoad NoSharedMemory \
    --twopassMode Basic \
    --readFilesCommand gunzip -c \
    --readFilesIn $(ls *.fastq.gz) \
    --outFileNamePrefix $home/STAR/$SID \
    --outFilterMultimapNmax $MAXmultimap \
    --outSAMtype $OUT $TYPE \
    --limitBAMsortRAM $BAMsortRAM \
    --outFilterMismatchNoverReadLmax $MAXmismatch\
    --quantMode GeneCounts
    cd ..
done

echo "
 `date` - Finished STAR mapping with Max Multimap = $MAXmultimap
"
# --------------------------------------------------------

# Now organize the output in ./STAR
cd $home/STAR
# Clean up
for d in $(ls -d */); do cd $d && rm -rf SJ.out.tab _STARgenome/ _STARpass1/ && head Log.out -n 99 > ../run.info && rm Log.out && cd ..; done

# --------------------------------------------------------

### Extract gene counts ###
# The following one liner code is designed to filter STAR's --quantMode GeneCounts output (FILENAME.ReadsPerGene.out.tab).
# The filtered output txt tables can then be directly used as input for CountMatrixGenerator.R script.
# The counts from STAR coincide with those produced by htseq-count with default parameters.
# STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options:
#column 1: gene ID
#column 2: counts for unstranded RNA-seq => this is the column we are interested for single read sequencing
#column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

# CODE: -------------------------------------
# cut: filter 1st and 2nd column of input file (specified by --fields (-f) 1-2)
# tail: remove first four rows (+NUM to output starting with line NUM)
# cut -f 1-2 *out.tab | tail -n +5 > GeneCounts.txt

# Go to your project dir, creating a GeneCounts folder
mkdir $home/GeneCounts
STR=ReadsPerGene.out.tab #string to match from gene counts
for d in $(ls -d */); do cd $d && D=$(echo $d | sed 's,/,,') && cut -f 1-2 $STR | tail -n +5 > ../../GeneCounts/$D'_GeneCounts.txt' && cd ..; done

# clip file names down to sample ID
cd ../GeneCounts
for VAR in $(ls *_GeneCounts.txt); do mv -vnT $VAR ${VAR/*R/R}; done #matches everything from p to R in string and replaces by R
for VAR in $(ls *_GeneCounts.txt); do mv -vnT $VAR ${VAR/-*_GeneCounts.txt/_GeneCounts.txt}; done

cd $home
conda deactivate

#### END OF SCRIPT ####