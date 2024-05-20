# Load necessary libraries
suppressPackageStartupMessages(library('DESeq2'))
library("BiocParallel")
register(MulticoreParam(4))
library("matrixStats")

# Function to extract the first variable from the design formula
extract_first_variable <- function(input_string) {
  match <- regmatches(input_string, regexpr("(?<=~)\\w+", input_string, perl = TRUE))
  if (length(match) == 0) {
    stop("No variable found in the design formula.")
  } else {
    return(match)
  }
}

# Function to run DESeq2 analysis
run_DE <- function(indata, insample, rowm, design, outfile, reference) {
  # Read the sample information and expression data
  insample <- read.table(insample, sep='\t', header=TRUE, row.names=1)
  indata <- read.table(indata, sep='\t', header=TRUE, row.names=1)
  
  # Ensure the column names of indata match the row names of insample
  if (!all(colnames(indata) == rownames(insample))) {
    stop('indata and meta table do not have matching names')
  }
  
  # Create DESeq dataset
  dds_pre <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design=as.formula(design))
  dds_pre <- DESeq(dds_pre)
  
  # Normalize counts and filter based on row medians
  norm_counts <- counts(dds_pre, normalized=TRUE)
  indata_filter <- indata[rowMedians(norm_counts) >= 10,]
  
  # Extract the variable from the design and relevel the factor
  de_variable <- extract_first_variable(design)
  insample[[de_variable]] <- factor(insample[[de_variable]])
  insample[[de_variable]] <- relevel(insample[[de_variable]], ref=reference)
  
  # Create a new DESeq dataset with the filtered data
  dds2 <- DESeqDataSetFromMatrix(countData = indata_filter, colData = insample, design=as.formula(design))
  dds2 <- DESeq(dds2)
  
  # Get results and save to output file
  resultDESeq2 <- results(dds2)
  res <- na.omit(resultDESeq2)
  write.table(res, outfile, sep='\t', quote=FALSE)
  
  # Generate MA plot and save as PDF
  pdf('plotMA.pdf')
  plotMA(resultDESeq2)
  dev.off()
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
indata <- args[1]
insample <- args[2]
rowm <- as.numeric(args[3])
design <- args[4]
reference <- args[5]
outfile <- args[6]

# Run DESeq2 analysis
run_DE(indata, insample, rowm, design, outfile, reference)
