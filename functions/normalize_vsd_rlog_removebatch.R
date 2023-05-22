suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('limma'))

args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
remove_b = args[3]
outfile = args[4]
norm_type = args[5]

remove_b <- gsub("batch_", "", remove_b)

normalize_remove_vsd_rlog_batchef <- function(indata, insample, remove_b, outfile, norm_type){
  indata = read.table(indata, sep='\t', header=T)
  rownames(indata) <- indata$X
  indata$X <- NULL
  keepRows = rowSums(indata >= 25) >= 2
  indata <- indata[keepRows,]
  
  insample <- read.table(insample, sep='\t', header=T)
  rownames(insample) <- insample$X
  insample$X <- NULL
  norm = 0
  batch <- insample$SeqTag
  
  if (remove_b == 'batch_tissue'){
    dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~SeqTag + tissue)
    remove_batch = 1
  } else if (remove_b == 'batch_notissue'){
    dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~SeqTag)
    remove_batch = 1
  } else if (remove_b == 'nobatch_tissue'){
    dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~tissue)
    remove_batch = 0
  } else if (remove_b == 'nobatch_notissue'){
    dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~1)
    remove_batch = 0
  }else{
    print('ERROR: type of performance not selected')
    remove_batch = 0
  }
  print('Performing transformation..')
  if(norm_type == c('rlog')){
      df_transf <- rlog(dds)
  } else if (norm_type == c('vsd')){
      df_transf <- varianceStabilizingTransformation(dds)
  } else if (norm_type == c('normalize')){
      dds <- estimateSizeFactors(dds)
      df_transf <- counts(dds, normalized=TRUE)
      norm = 2
      if (remove_batch == 1){
        remove_batch = 2
      }
  } else if (norm_type == c('quantilenorm_log2')){
      library(qsmooth)
      #print(assay(dds))
      print('qs')
      qs <- qsmooth(object = assay(dds), group_factor = insample$tissue)
      qs_data <- qsmoothData(qs) # extract smoothed quantile normalized data
      #Log2 the values
      write.table(log2(qs_data+1), paste(outfile,'_qs.tab', sep=''), sep='\t', quote = F)
      df_transf <- log2(qs_data+1)
      norm = 2
  } else{
      df_transf <- dds
  }
  print('done')
  print(norm_type)
  print(remove_b)
  if (remove_batch == 1){ 
    assay(df_transf) <- limma::removeBatchEffect(assay(df_transf), df_transf$SeqTag)
  }
  if (remove_batch == 2){
    print('remove_batch 2')
    df_transf <- limma::removeBatchEffect(df_transf, batch)
  }
  if (norm == 0){
    print(head(assay(df_transf)))
    write.table(assay(df_transf),outfile, sep='\t', quote = F)
  }else if(norm == 2){
    write.table(df_transf, outfile, sep='\t', quote = F)
  }else{
    #print(head(df_transf))
    write.table(assay(df_transf), outfile, sep='\t', quote = F)
  }
}

normalize_remove_vsd_rlog_batchef(indata, insample, remove_b, outfile, norm_type)
