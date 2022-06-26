#!/bin/bash
# Author: Hannes Reinwald
# Contact: hannes.reinwald@ime.fraunhofer.de

# README:
# Write details and general info about script here
# Desgined for gff3 file and dna.toplevel.fa
#
#
# -----------------------------------------------------

#### START OF SCRIPT ##################################

# path to dir containing ref genome fasta and gff3 file
home=$(pwd)

# activate conda env for STAR
source activate rnaseq

# Set STAR run params 
echo "
 What is your sequence read length? For 50 bp just type: 50"
read nbp
overhang=$[$nbp-1]  # read length-1 (i.e. 50 bp; overhang= 50-1 = 49) 


# PARAMETERS:
outdir="genome_index_files_${nbp}bp_GTF"
gtf=$(ls *.gtf)
fasta=$(ls *.dna.toplevel.fa)


# Compute optimal binbit size
BYT=$(wc -c $fasta | awk '{print $1}') #fasta byte size
CONT=$(grep -c "^>" $fasta) #Nbr of contigs
# log2 function
function log2 {
    local x=0
    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
        let x=$x+1
    done
    echo $x
}
z=$[$(log2 $[$BYT/$CONT])-1] #Ideal Binbit size rounded to one digit

echo "
 Give a genomeChrBinNbits size. 
 Recommend BinNbit size for your Input: " $z
read binbits 
# Ideal size for binbits = log2(genome size in bytes (unzipped!) / Number of contigs)
# i.e. genome size of 1.5 GB (1.5e9 bytes) & 993 contigs: log2(1.5e9/993) = ~ 20

# get fa file size with: wc -c *.fa | awk '{print $1}'
# get Nbr of contigs: grep -c "^>" *.fa

CPUs=$[$(nproc)-2]  # Number of available CPUs - 2

# run STAR genome index generator 
STAR --runThreadN $CPUs \
    --runMode genomeGenerate \
    --genomeDir ./$outdir \
    --genomeFastaFiles $fasta \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfile $gtf \
    --sjdbOverhang $overhang \
    --genomeChrBinNbits $binbits \
    --genomeSAindexNbases 12 #recomended by STAR for the size of the Dmagna genome 2.4
    
    # For gtf file use this instead:
    #--sjdbGTFtagExonParentTranscript transcript_id \
    #--sjdbGTFfile ./*.gtf \

### Change file permission to be accessible for everyone ### 
chmod -R g+rx $home/$outdir


echo "
	If everything worked fine, your genome index files are located in:
	$home/$outdir
	"

exit 0
### END OF SCRIPT ###