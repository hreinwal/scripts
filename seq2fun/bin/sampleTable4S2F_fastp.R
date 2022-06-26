# This script is designed to combine coldata info with fastq file names to generate the 
# required sampleTable for downstream Seq2Fun analysis pipeline. 
# Sample table must consist of 3 columns (sample prefix name (sample01), forward reads name (sample01_R1.fq.gz), group info (control) 
# for single-reads: 3 columns (sample prefix name (sample01), forward reads name (sample01_R1.fq.gz), group info (control)
# for paired-end: 	4 columns (sample prefix name (sample01), forward reads (sample01_R1.fq.gz), reverse reads (sample01_R2.fq.gz), group info (control) for paired-end reads.
# The columns must be separated by tab!

## coldata import ## ----------------------------------------------------------------------
tmp = list.files(full.names = T, pattern = "[Cc]oldata")
if(grepl(".csv",tmp)==T){
  coldata = read.csv2(file = tmp, row.names = 1, header = T)
  # In case csv file was not differently exported run read.csv
  if (ncol(coldata) <= 1) {
    coldata = read.csv(file = tmp, row.names = 1, header = T)
  }
} else {
  # if not ending with csv import with read.delim2
  coldata = read.delim2(file = tmp, row.names = 1, header = T)
  if (ncol(coldata) <= 1) {
    # if file was differently formated use read.delim
    coldata = read.delim(file = tmp, row.names = 1, header = T)
  }
}
# Select sample group column (Condition)
Condition = grep("[Cc]ondition|[Tt]reatment",colnames(coldata))#"Condition"
stopifnot(length(Condition) == 1) # Sum check
colnames(coldata)[Condition] = "Condition"
coldata$sampleID = row.names(coldata)
coldata = droplevels(coldata[which(coldata$Condition != ""), # filtering for non empty rows
                              which(!(grepl("X.",colnames(coldata))) == T)]) # filtering cols
coldata$Condition = gsub(" ","",coldata$Condition) #removing any spaces in Condition levels
####################

fastq = list.files(path = "./seqQC/fastp", full.names = T, pattern = "[.]fastq", recursive = T)
outDir= "seq2fun" # output dir for downstream seq2fun analysis
dir.create(outDir, showWarnings = F)

ls = list()# output list to compile information in.
for(id in coldata$sampleID) {
  fa = grep(id, fastq, value = T)
  cond = coldata[id,"Condition"]
  out = paste(".",outDir,id, sep="/")
  if(length(fa) > 1){
    # paired-reads: if fa contains two file names = paired end reads --> requires different output file type
    tmp = c(out, sort(fa), cond)
    names(tmp) = c("ID","read1","read2","group")
    ls[[id]] = tmp
  } else {
    # single-read
    tmp = c(out, fa, cond)
    names(tmp) = c("ID","read1","group")
    ls[[id]] = tmp
  }
}
df = do.call(rbind.data.frame, ls)
colnames(df) = NULL

# Export sampletable.tab file
write.table(df, paste0(basename(getwd()),"_sampleTable.tab"),
            sep = "\t", quote = F, row.names = F, col.names = F)
### END OF SCRIPT ### 