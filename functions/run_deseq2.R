#Run deseq2
suppressPackageStartupMessages(library('DESeq2'))

args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
batch = args[3]
intype = args[4]
norm_type = args[5]
outfile = args[6]
rowm = args[7]
paired_samples = args[8]
print(indata)
print(insample)


run_DE <-  function (indata, insample, batch, intype, norm_type, outfile, rowm, paired_samples) {
  #The insample table requires the column "SeqTag" and tissue
  insample <- read.table(insample, sep='\t', header=T)
  rownames(insample) <- insample$X
  insample$X <- NULL
  
  indata <- read.table(indata, sep='\t', header=T)
  rownames(indata) <- indata$X
  indata$X <- NULL


  keepRows = rowSums(indata >= rowm) >= 2
  #indata <- indata[rowSums(indata >= 25) >= 0.9*dim(indata)[2],]
  
    indata <- indata[keepRows,]
  print('indata after fitlering')
  print(head(indata))
  
  if (all(colnames(indata) == rownames(insample))){
    print("Both the column names matches in position and names")
  } else{
    stop('indata and meta table does not have matching names')
  }
  
  if(batch == c('batch')){
    print('batch')
    if(intype == c('tissue')){
      dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~SeqTag + tissue)
    } else {
      dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~SeqTag + type)
    }
  } else{
    print('no batch')
    if(intype == c('tissue')){
      if(paired_samples == 'TRUE'){
        insample$sll_id <- as.factor(insample$sll_id)
        dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~ sll_id + tissue)
      } else {
      dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~tissue)
      }
      
    } else {
      dds <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design= ~type)
    }
  }
  #if(norm_type == c('rlog')){
  #  dds <- rlog(dds)
  #} else if (norm_type == c('vsd')){
  #  dds <- varianceStabilizingTransformation(dds)
  #} else if (norm_type == c('normalize')){
  #  dds <- estimateSizeFactors(dds)
  #  dds <- counts(dds, normalized=TRUE)
  #  norm = 1
  #} else if (norm_type == c('quantilenorm_log2')){
  #  print('Dont run quantile norm with DESeq')
  #  stop('break')
  #  }
  
  print('DESeq')
  #print(dds)
  dds <- DESeq(dds)
  print('done')
  resultDESeq2 = results(dds)
  #alpha <- 0.05 # Threshold on the p-value
  res <- na.omit(resultDESeq2)
  #res <- res[res$padj <= alpha,]
  write.table(res, outfile, sep='\t', quote = F)
  
}


run_DE(indata, insample, batch, intype, norm_type, outfile, rowm, paired_samples)
