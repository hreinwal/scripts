#!bin/bash
# Author:  Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# This script performs fastp dataQC and filtering. 
# Then runs fastQScreen & fastQC on fastp filtered fastq files.

# Navigate to your project folder (/srv/attract/seq_files/project_folder) then execute this script.
home=$(pwd)
CPUs=$[$(nproc)-2] # 14
source $home/config

##########################------------------------------------------------
### fastp on RAW fastq ###
##########################

## fastp PARAMETERS ##
#out=../../fastp
out=$home/fastp
Nlim=$maxNperRead  #~ 10% of total read length (e.g. for 150bp Nlim=15)
RLN=$minReadLength #minimum read length
QPh=$minQphred #minimum Quality phred score
#FQ=$(find . -name "*.fastq.gz") #files to process

# load and prepare env
source activate rnaseq
cd $home/raw_reads #input file location
mkdir $home/fastp #output file location

### fastp ###
path2adapters=$(dirname $(realpath $0))

for d in $(ls -d */); do
	cd $d && R1=$(ls *_1.fastq.gz) && R2=$(ls *_2.fastq.gz) && \
	dir=$(echo $d | sed -e 's/_lib[0-9].*$//' -e 's/NG-[0-9].*_R/R/') && \
	r1=$(echo $R1 | sed -e 's/_lib[0-9].*$/_1.fastq.gz/' -e 's/NG-[0-9].*_R/R/') && \
	r2=$(echo ${r1/_1.fastq.gz/_2.fastq.gz}) && \
	mkdir $out/$dir
	printf "\nRunning fastp for read-pair:\n$r1 \n$r2 \n... \n"
	fastp -i $R1 -I $R2 -o $out/$dir/$r1 -O $out/$dir/$r2 \
		--trim_poly_g --thread $CPUs \
		--qualified_quality_phred $QPh \
		--n_base_limit $Nlim \
		--detect_adapter_for_pe \
		--length_required $RLN \
		--adapter_fasta $path2adapters/adapters.fasta
	rm -fv fastp.html && mv -fv fastp.json $out/$dir
	cd ..
done
# to deactivate quality filter use -Q parameter!
# https://github.com/OpenGene/fastp#output-splitting 

# fix sample labeling issue for multiQC report; replace read1 string in json report
cd $home/fastp
#sed -i 's/_1.fastq.gz/.fastq.gz/g' $(find . -name fastp.json)
sed -i 's/-i NG-[0-9]*[_R]*/-i R/g' $(find . -name fastp.json)
sed -i 's/_lib[0-9].*[1-2].fastq.gz -I/.fastq.gz -I/g' $(find . -name fastp.json)


########################################----------------------------------
## FastQ-Screen post fastp for PE Seq ##
########################################

# Specify path to config file
# Config file specifies genomes to map against. Look at the file for details. 
# All bowtie2 indexed ref genomes to map against can be found under:
# /srv/attract/ref_genome/FastQ_Screen_Genomes
configFile=$FQS_configFile

# go fastp filtered reads
cd $home/fastp

# output dir
mkdir $home/fastQ_screen
touch $home/fastQ_screen/fastQscreen.log

# Check if provided file exists
if [ -f "$configFile" ]; then
	echo "$configFile exists." 2>&1 | tee -a $home/fastQ_screen/fastQscreen.log
else
	echo "$configFile does not exists. 
	Exiting script" 2>&1 | tee -a $home/fastQ_screen/fastQscreen.log
	cd $home
	exit 0
fi

### FastQ-Screen ###
source activate rnaseq
FQ=$(find . -name "*.fastq.gz")
time parallel -j 9 --eta 'fastq_screen --threads 2 --conf '$configFile' --aligner bowtie2 --subset '$sampleN' --force --outdir ../fastQ_screen ' ::: $FQ 2>&1 | tee -a ../fastQ_screen/fastQscreen.log

# Clean up
rm -fv $home/fastQ_screen/*screen.png $home/fastQ_screen/*screen.html



##################################-----------------------------------------
## FastQC post fastp for PE Seq ##
##################################

# Activate Conda Env
source activate multiqc

### FastQC ###
cd $home/fastp #input file location
mkdir $home/fastQC #output location
touch $home/fastQC/fastqc.progress
fastqc --outdir $home/fastQC/ --threads $CPUs --noextract --format fastq --kmers 7 \
	$(find . -name "*.fastq.gz") 2>&1 | tee -a $home/fastQC/fastqc.progress

# clean up
rm -fv $home/fastQC/*fastqc.html

## Organize final fastq sequence QC reports ## 
cd $home
DIR=$(ls -d */ | egrep 'fastQ_screen|fastQC|fastp') #directories to move 
mkdir $home/seqQC
mv -fv $DIR $home/seqQC

conda deactivate && conda deactivate && cd $home
##### END OF SCRIPT ###### 