library('edgeR')
library('DESeq2')

args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
remove_b = args[3]
outfile = args[4]



normalize_remove_batchef <- function(indata, insample, remove_b, outfile){
  indata = read.table(indata, sep='\t', header=T)
  
  rownames(indata) <- indata$X
  indata$X <- NULL
  #print(head(indata))
  insample <- read.table(insample, sep='\t', header=T)
  #print(head(insample))
  rownames(insample) <- insample$X
  insample$X <- NULL
  
  y <- DGEList(counts = indata,genes = rownames(indata))
  #y$samples
  #print(y)
  # Filtering lowly expressed genes
  keep <- rowSums(cpm(y)>2) >= 2 # A gene should be expressed in both samples from same patient
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  batch <- insample$SeqTag
  
  if(remove_b == c('batch')){
    logCPM <- cpm(y,log = TRUE, prior.count = 2)
    log2CPM_norm <- removeBatchEffect(logCPM, batch = batch)
  } else{
    #Just normalize
    log2CPM_norm <- cpm(y,log = TRUE, prior.count = 2)
  }
  write.table(log2CPM_norm, outfile, sep='\t', quote = F)
}

normalize_remove_batchef(indata, insample, remove_b, outfile)