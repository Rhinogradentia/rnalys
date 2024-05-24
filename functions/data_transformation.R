suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('limma'))
suppressPackageStartupMessages(library('qsmooth'))

##DESEQ2 is used for variance stabalizing and rlog transformation
# Function to normalize RNA-seq data using various methods and optionally remove batch effects


normalize_remove_vsd_rlog_batchef <- function(indata, insample, remove_b, outfile, norm_type){
  
 #normalize_remove_vsd_rlog_batchef <- function(indata, insample, remove_b, norm_type) {
  # Ensure indata and insample are data frames directly
  if (is.character(indata)) {
    # Assuming 'indata' is a file path to a tab-separated file (TSV)
    indata <- read.table(indata, sep='\t', header=TRUE, row.names=1)
  
    # Check if the read data is numeric; if not, convert it
    if (!is.matrix(indata) || mode(indata) != "numeric") {
        indata <- data.matrix(indata)
        if (any(is.na(indata))) {
            stop("Some values could not be converted to numeric.")
        }
    }
} else if (is.matrix(indata)) {
    # Directly use the matrix but ensure it is numeric
    indata <- indata
    if (mode(indata) != "numeric") {
        indata <- as.numeric(indata)
    }
} else if (is.data.frame(indata)) {
    # If 'indata' is a data frame, convert to a matrix
    indata <- data.matrix(indata)
} else {
    stop("indata should be a matrix or a path to a file containing the matrix.")
}
  print(head(indata))
  
  if (is.character(insample)) {
    insample <- read.table(insample, sep='\t', header=TRUE, row.names=1)
  } else if (is.data.frame(insample)) {
    insample <- insample
  } else {
    stop("insample should be a data frame or a path to a CSV file containing the sample information.")
  }
  
  #insample <- read.table(insample, sep='\t', header=T)
  #rownames(insample) <- insample$X
  #insample$X <- NULL
  
  
  norm = 0
  #ddbatch <- insample$SeqTag
  
  dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~1)
  
  
  #print('Performing transformation..')
  if(norm_type == c('rlog')){
      df_transf <- rlog(dds)
  } else if (norm_type == c('vsd')){
      df_transf <- varianceStabilizingTransformation(dds)
  } else if (norm_type == c('normalize')){
      dds <- estimateSizeFactors(dds)
      df_transf <- counts(dds, normalized=TRUE)
      norm = 2

  } else if (norm_type == c('quantilenorm_log2')){
      library(qsmooth)
      #print(assay(dds))
      #print('qs')
      qs <- qsmooth(object = assay(dds), group_factor = insample$tissue)
      qs_data <- qsmoothData(qs) # extract smoothed quantile normalized data
      #Log2 the values
      write.table(log2(qs_data+1), paste(outfile,'_qs.tab', sep=''), sep='\t', quote = F)
      df_transf <- log2(qs_data+1)
      norm = 2
  } else{
      df_transf <- dds
  }
  #print('done')
  #print(norm_type)
  #print(remove_b)
  if (remove_b != 'nobatch'){ 
    assay(df_transf) <- limma::removeBatchEffect(assay(df_transf), df_transf[[remove_b]])
  }
  
  if (norm == 0){
    #print(head(assay(df_transf)))
    write.table(assay(df_transf),outfile, sep='\t', quote = F)
  }else if(norm == 2){
    write.table(df_transf, outfile, sep='\t', quote = F)
  }else{
    #print(head(df_transf))
    write.table(assay(df_transf), outfile, sep='\t', quote = F)
  }
}

if (!interactive()) {
    
    args <- commandArgs(trailingOnly = TRUE)
    indata = args[1]
    insample = args[2]
    remove_b = args[3]
    outfile = args[4]
    norm_type = args[5]
    #remove_b <- gsub("batch_", "", remove_b)
    normalize_remove_vsd_rlog_batchef(indata, insample, remove_b, outfile, norm_type)
}