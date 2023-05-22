#Remove batcheffect

library('edgeR')

args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
outfile = args[3]

remove_batchef <- function(indata, insample, remove_b, outfile){
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
  print(y)
  # Filtering lowly expressed genes
  keep <- rowSums(cpm(y)>2) >= 2 # A gene should be expressed in both samples from same patient
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  logCPM <- cpm(y,log = TRUE,prior.count = 2)
  batch <- insample$SeqTag
  #print(batch)
  print(dim(logCPM))
  print(dim(insample))
  logCPMAdj <- removeBatchEffect(logCPM,batch = batch)
  #print(head(logCPMAdj))
  write.table(logCPMAdj, outfile, sep='\t', quote = F)
}

remove_batchef(indata, insample, outfile)
