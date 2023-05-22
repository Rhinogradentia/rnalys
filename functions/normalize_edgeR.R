
library('edgeR')
#library('dplyr')
args <- commandArgs(trailingOnly = TRUE)
#indata countdata
#Insample containing the table of patients
indata = args[1]
insample = args[2]
outfile = args[3]

#indata = "../H_2019/merged_counts_sll.tab" 
#indata = 'LV_RV_raw_count.tab'
#insample = "LV_RV_sample.tab"               
#outfile = "test.tab"               



normalize_deseq <- function (indata, insample, outfile) {
  insample <- read.table(insample, sep='\t', header=T)
  #print(head(insample))
  print(insample$id_tissue)
  indata <- read.table(indata, sep='\t', header=T)
  rownames(indata) <- indata$X
  indata$X <- NULL
  indata <- indata[,insample$id_tissue]

  #print(head(indata2))
  y= DGEList(counts=indata, genes=rownames(indata), group = factor(insample$id_tissue))
  print(colnames(indata))
  #filter lowly expressed
  keep2 <- rowSums(cpm(y)>2) >= 2
  y <- y[keep2, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y = estimateCommonDisp(y, verbose=TRUE)
  
  logCPM2_2 <- cpm(y,log = TRUE,prior.count = 2)
  write.table(logCPM2_2, outfile, sep='\t', quote = F)
}


normalize_deseq(indata, insample, outfile)


