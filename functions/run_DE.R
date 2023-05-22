#DE analysis
source('./functions/DE_analysis_nodesign.R')
library('edgeR')
args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
outfile = args[3]


run_DE <-  function (indata, insample, outfile) {
  insample <- read.table(insample, sep='\t', header=T)
  indata <- read.table(indata, sep='\t', header=T)
  rownames(indata) <- indata$X
  indata$X <- NULL
  rownames(insample) <- insample$id_tissue
  insample$X <- NULL
  #print(insample)
  #print(head(indata))
  #print(dim(insample))
  #print(dim(indata))
  y1= DGEList(counts=indata, genes=rownames(indata), group = factor(insample$tissue))
  keep2 <- rowSums(cpm(y1)>2) >= 2
  y1 <- y1[keep2, , keep.lib.sizes=FALSE]
  y1 <- calcNormFactors(y1)
  y1 = estimateCommonDisp(y1, verbose=TRUE)
  group <- factor(insample$tissue)
  y1 <- estimateDisp(y1)

  fit <- glmFit(y1)
  DEgenes_LV_RV<- DE_analysis_nod(fit,y = y1,prefix = "--",foldchange = 0.5,output = outfile) 

  write.table(DEgenes_LV_RV, outfile, sep='\t', quote = F)
  print(DEgenes_LV_RV)
}
run_DE(indata, insample, outfile)