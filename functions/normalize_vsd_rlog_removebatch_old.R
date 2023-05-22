library('DESeq2', verbose = F)
library('limma', verbose = F)

args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
remove_b = args[3]
outfile = args[4]
norm_type = args[5]



normalize_remove_vsd_rlog_batchef <- function(indata, insample, remove_b, outfile, norm_type){
  indata = read.table(indata, sep='\t', header=T)
  
  rownames(indata) <- indata$X
  indata$X <- NULL
  #print(head(indata))
  insample <- read.table(insample, sep='\t', header=T)
  #print(head(insample))
  rownames(insample) <- insample$X
  insample$X <- NULL

  batch <- insample$SeqTag
  
  if(remove_b == c('batch')){
    dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~SeqTag + tissue)
    if(norm_type == c('rlog')){
      df_transf <- rlog(dds)
    } else if (norm_type == c('vsd')){
      df_transf <- varianceStabilizingTransformation(dds)
    } else{
      df_transf <- dds
    }
    #Remove batch effect for pca
    assay(df_transf) <- limma::removeBatchEffect(assay(df_transf), df_transf$SeqTag)
  }else if (remove_b == c('batch_notissue')){
    dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~1)
    if(norm_type == c('rlog')){
      df_transf <- rlog(dds)
    } else if (norm_type == c('vsd')){
      df_transf <- varianceStabilizingTransformation(dds)
    } else{
      df_transf <- dds
    #Remove batch effect for pca
    print('remove batch effect batch_notissuse')
    print(indata)
    assay(df_transf) <- limma::removeBatchEffect(assay(df_transf), df_transf$SeqTag)
  }
  }else{
    if (remove_b == 'normalize'){
      dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~tissue)
      if(norm_type == c('rlog')){
        df_transf <- rlog(dds)
      } else{
        df_transf <- varianceStabilizingTransformation(dds)
      }
      
    } else{
      dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~1)
      if(norm_type == c('rlog')){
        df_transf <- rlog(dds)
      } else{
        df_transf <- varianceStabilizingTransformation(dds)
      }
    }
    assay(df_transf) <- limma::removeBatchEffect(assay(df_transf), df_transf$SeqTag)
  }
  write.table(assay(df_transf), outfile, sep='\t', quote = F)
}

normalize_remove_vsd_rlog_batchef(indata, insample, remove_b, outfile, norm_type)