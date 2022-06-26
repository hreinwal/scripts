# R_DESeq2_DEG_analysis

Repository containing R scripts for differentially expressed gene analysis with DESeq2.

Consider this analysis as a first quality assessement on your data. i.e. if differentially gene expression actually occured, correlation between biol. replicates etc.
The scripts are desgined to run the analysis on its own using a gene count matrix as input. 

**Requiered Input:**
*  **CountMatrix** (Rows = GeneIDs, Cols = SampleID; can be automatically generated from STAR's Gene count files via the *`CountMatrix_generator.R`* script. As input file types, .txt and .csv are supported.
*  **coldata** (Rows = SampleID, Cols = Metainfo: [*"Condition";"Substance";"Tank";"Conc";"SamplingDate",...*] As input file types, tab delim and .csv are supported.

**Output**
* Data Quality plots (pvalue & LFC distr, biol. sample correlation, Gene count transformation, Norm. count rlog transformation)
* PCA & Hclustering 
* MA and Vulcano Plots (based on 90% quantile LFC cut off and padj value cut off)
* Log2FC correlation analysis on DEGs (based on different filtering parameters)
* Heatmaps
* Xcel list of DESeq2 results (HighExposure vs NC, LowExposure vs NC)
* Overlapping DEGs (inkl LFC cut off) as Xcel list 
* R Session Info

**Run Script:**

Navigate to the folder containing *CountMatrix* and *coldata* with `setwd("path/to/my/files")` in R. Then execute the R script via: `source("path/to/Rscript.R")`
