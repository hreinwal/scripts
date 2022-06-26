#!/bin/bash

# Author:  Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# This script performs post mapping alignment QC with samtools, RSeQC and Qualimap
# First BAM files are indexed via samtools and quality checked.
# Subsequently RSeQC and Qualimap analysis is performed. (requires prior index and BAM12 file)
home=$(pwd)
CPUs=$[$(nproc)-2]
source $home/config
STR=Aligned.sortedByCoord.out.bam #String to match to find corresponding BAM files 

# Not the most elegant way but will do for now ... 
# Will create a separated tmp dir to store all BAM files in there. Easier to process that way
cd $home/STAR #BAM file directory
mkdir $home/tmp
for bam in $(ls *$STR); do 
	cp -f $bam $home/tmp/${bam/$STR/.out.bam}
done

cd $home/tmp # Directory in which BAM files for alignment QC are tmp located 
source activate rseqc # Environment to run samtools, RSeQC & Qualimap with 


###########################
### samtools & indexing ###
########################### ----------------------------------------

mkdir $home/samtools $home/samtools/flagstat $home/samtools/idxstats
STR=.out.bam
# 1st indexing
time parallel -j $CPUs --eta 'samtools index {}' ::: *$STR
# 2nd flagstat / idxstats
time parallel -j $CPUs --eta 'samtools flagstat {} >{}.flagstat | samtools idxstats {} >{}.idxstats' ::: *$STR 2>&1 | tee -a $home/samtools/samtools.report
# organize output
for f in *.flagstat; do mv -f $f $home/samtools/flagstat/${f/$STR/}; done
for f in *.idxstats; do mv -f $f $home/samtools/idxstats/${f/$STR/}; done



################
### Qualimap ###
################ ---------------------------------------------------

## PARAMETER ##
GTF=$GTF4QC

# Check if provided file exists
if [ -f "$GTF" ]; then
	echo "$GTF exists."
else
	echo "$GTF does not exists. 
	Exiting script"
	cd ..
	exit 0
fi

### qualimap ###
cd $home/tmp #input dir
mkdir $home/qualimap #output dir
time parallel -j $CPUs --eta 'qualimap rnaseq -outdir ../qualimap/{} \
	-a proportional -p non-strand-specific -bam {} -gtf '$GTF' \
	--java-mem-size=10G' ::: *$STR

# organize output in qualimap/
cd $home/qualimap && \
for d in *; do mv -v $d ${d/$STR/}; done
# remove everything not needed for multiQC report
RMV=$(find | egrep 'css|images_qualimapReport|qualimapReport.html')
rm -fr $RMV



############
##  RSeQC ##
############---------------------------------------------------------

## ṔARAMETER ##
bedPath=$BED4QC

cd $home/tmp
# with tin.py => takes tooooo long!
#time parallel -j $nCPU --eta 'geneBody_coverage.py -i {} -r '$bedPath' -o ../RSeQC/{} | tin.py -i {} -r '$bedPath' | infer_experiment.py -i {} -r '$bedPath' >{}_inferExperiment.txt | read_distribution.py -i {} -r '$bedPath' >{}_readDistribution.txt' ::: *$MATCH
# Without tin.py => GenebodyCov takes toooooo long !!!
#time parallel -j $nCPU --eta 'geneBody_coverage.py -i {} -r '$bedPath' -o {} | infer_experiment.py -i {} -r '$bedPath' >{}.infEx.txt | read_distribution.py -i {} -r '$bedPath' >{}.rDis.txt' ::: *$MATCH
# Without tin.py & geneBody_coverage.py
time parallel -j $CPUs --eta 'infer_experiment.py -i {} -r '$bedPath' >{}.infEx.txt | \
 read_distribution.py -i {} -r '$bedPath' >{}.rDis.txt' ::: *$STR

# organize output in tmp/
mkdir $home/RSeQC $home/RSeQC/inferExperiment $home/RSeQC/readDistribution
F1=$(ls *.infEx.txt) && F2=$(ls *.rDis.txt)
for f in $F1; do mv -vf $f $home/RSeQC/inferExperiment/${f/.out.bam.infEx/}; done
for f in $F2; do mv -vf $f $home/RSeQC/readDistribution/${f/.out.bam.rDis/}; done


## FINAL Cleanup ## -------------------------------------------------
cd $home && rm -fr tmp/
DIR=$(ls -d */ | egrep 'qualimap|RSeQC|samtools') #directories to move 
mkdir $home/alignQC
mv -fv $DIR $home/alignQC


###############
##  multiQC  ##
############### ------------------------------------------------------
source activate multiqc

### MultiQc ###
NAME="`basename "$(pwd)"`"
multiqc --zip-data-dir --outdir . --filename ${NAME}"_multiQC" -f .  
#the "." specificies the current wd! ~ to "./"

conda deactivate && conda deactivate # Deactivate the loaded conda env again

### END OF SCRIPT ###