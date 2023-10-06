
# ChIPseeker is an R package for annotating ChIP-seq data analysis. 
# It supports annotating ChIP peaks and provides functions to visualize ChIP peaks coverage over chromosomes and profiles of peaks binding to TSS regions. 
# Comparison of ChIP peak profiles and annotation are also supported. Moreover, it supports evaluating significant overlap among ChIP-seq datasets. 
# Currently, ChIPseeker contains 17,000 bed file information from GEO database. 
# These datasets can be downloaded and compare with user’s own data to explore significant overlap datasets for inferring co-regulation or transcription factor complex for further investigation.

# source > https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

# **Note:** Remember to create a "data/" and a "output/" folder in the current directory.

################################################################################

# Install packages.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
library(ChIPseeker)

################################################################################
# Getting the current working directory and setting paths for input and output directories.
wd <- getwd()
path_to_input= file.path(wd, "data")
path_to_output= file.path(wd, "output")
# Defining a file pattern (.txt) to search for in the input directory.
pattern=".txt" #.narrowPeak
# Loading a transcript database (TxDb.Hsapiens.UCSC.hg19.knownGene) for the human genome.
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# Listing files in the input directory that match the specified pattern.
samplefiles <- list.files(path = path_to_input, pattern = pattern, full.names = T)
samplefiles <- as.list(samplefiles)
# Extracting filenames from the full file paths and removing the file extension.
filenames <- gsub(file.path(wd, "data", "*"), "", samplefiles) 
filenames <- gsub("*.narrowPeak", "", filenames)
names(samplefiles) <- filenames
# Annotating peaks in each input file using the annotatePeak function and storing the results in a list (peakAnnoList).
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, verbose = F, tssRegion=c(-2500, 2500))
# Creating a PDF file for an annotation plot (annotations_plot.pdf) and plotting the annotations using plotAnnoBar.
pdf(file = paste0(path_to_output, "/annotations_plot.pdf"), height = 3)
plotAnnoBar(peakAnnoList)
dev.off()
# Defining a custom function SaveTable to save the annotation data as separate Excel files for each input file.
SaveTable <- function(n) {
  anno <- data.frame(peakAnnoList[[n]]@anno)
  print(str(anno))  # Print the structure of anno
  write.table(anno, file = paste0(wd, "/output/", filenames[n], "_anno.xls"), sep = "\t")
}
# Using lapply to iterate through the list of input files and save their annotation data as Excel files in the output directory.
lapply(1:length(samplefiles), SaveTable)
################################################################################
