### PARAMETERS to play with ###
searchMode="tGREEDY" #tMEM or tGREEDY (default); Greedy mode is designed for organisms who do not have reference sequences in database.
MinScore="70" # minscore; Specify minimum matching score DEFAULT 80
MaxMismatch="5" # --mismatch; Max Nbr of AA mismatches (DEFAULT 2)
maxNbase="10" # --n_base_limit; max Nbr of N bases in nt seq (DEFAULT 5)
ntMinLen="60" # --length_required; min nt seq length to be used for translation. Must be < than readL (DEFAULT 60)
AAminLen="19" # --minlength; min length of AA seq length (DEFAULT 19 (tGreedy); 13 (tMEM)) - Must be < than $readL/3 (if front clipped: < ($readL-6)/3)
AAmaxLen="60" # --maxtranslength; Maximum cutoff length of translated peptides (DEFAULT 60)
#For PE only
OverlLen="22" # --overlap_len_require; only for paired-end, min overlap between read1 & read2

###############################

# Reference Database #
protDbIndex="/srv/attract/git/seq2fun_db/insects/insects_v2.0.fmi" #(--tfmi)
geneMap="/srv/attract/git/seq2fun_db/insects/insects_annotation_v2.0.txt"
S2F=/srv/attract/git/Seq2Fun/bin/seq2fun #path to seq2fun
CPUs=14 #cpus for processing


## Go to your project directory and run the following lines: 
# Picks a single paired read sample
R1=$(find -name *_1.fastq.gz | head -n 1)
R2=$(echo ${R1/_1.fastq.gz/_2.fastq.gz})

# Subsample these samples (otherwise test runs might take too long)
SR=500000 #500k subsampled reads
mkdir test_run

source activate rnaseq
for FQ in $(echo $R1 $R2); do
	fq=$(echo $FQ | sed -e 's/^.*NG-/NG-/g' -e 's/_lib[0-9].*_/_/' -e 's/^.*_R/R/')
	seqtk sample -s100 $FQ $SR > test_run/${fq/.fastq.gz/.fastq}
	gzip -vf test_run/*.fastq
done

## Depending on which read type you are working with, now run paired or single read job
R1=$(ls ./test_run/*_1.fastq.gz)
R2=$(ls ./test_run/*_2.fastq.gz)

# PE
$S2F -i $R1 -I $R2 --prefix ./test_run/testRun --trim_front1 7 --trim_front2 7 --overlap_len_require $OverlLen \
	--thread $CPUs --profiling \
	--mode $searchMode --tfmi $protDbIndex --genemap $geneMap --length_required $ntMinLen \
	--minscore $MinScore --minlength $AAminLen --mismatch $MaxMismatch \
	--maxtranslength $AAmaxLen -V #--allFragments

# SR
$S2F -i $R1 --prefix ./test_run/testRun --trim_front1 7 \
	--thread $CPUs --profiling \
	--mode $searchMode --tfmi $protDbIndex --genemap $geneMap --length_required $ntMinLen \
	--minscore $MinScore --minlength $AAminLen --mismatch $MaxMismatch \
	--maxtranslength $AAmaxLen -V #\
	#--allFragments \
	#--trim_front1 6 --trim_front2 6 \