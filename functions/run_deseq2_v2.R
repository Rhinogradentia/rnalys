#Run deseq2
suppressPackageStartupMessages(library('DESeq2'))
library("BiocParallel")
register(MulticoreParam(4))

args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
rowm = args[3]
design = args[4]
outfile = args[5]
reference = args[6]

#print(indata)
#print(insample)
print(paste('design:', design, sep=' '))
print(design)

extract_variable <- function(formula_str) {
  pattern <- ".*\\+\\s*([^\\s]+)\\s*$|^\\s*~\\s*([^\\s]+)\\s*$"
  matches <- regmatches(formula_str, regexpr(pattern, formula_str, perl = TRUE))
  
  if (length(matches) > 0) {
    matched_vars <- matches[[1]][-1]
    variable_str <- matched_vars[!is.na(matched_vars)][1]
  } else {
    variable_str <- ""
  }
  
  return(variable_str)
}

run_DE <-  function (indata, insample, rowm, design, outfile, reference) {
  #The insample table requires the column "SeqTag" and tissue
  insample <- read.table(insample, sep='\t', header=T)
  rownames(insample) <- insample$X
  insample$X <- NULL
 

  if (grepl('bin4', design, fixed = TRUE) == TRUE){
    insample$age_bin4 <- cut(insample$age, b=4)
  }
  
  if (grepl('bin5', design, fixed = TRUE) == TRUE){
    insample$age_bin5 <- cut(insample$age, b=5)
  }
  
  
  indata <- read.table(indata, sep='\t', header=T)
  rownames(indata) <- indata$X
  indata$X <- NULL
  
  
  #keepRows = rowSums(indata >= rowm) >= 3
  #keepRows <- rowSums(indata) >= rowm
  #indata <- indata[rowSums(indata >= 25) >= 0.9*dim(indata)[2],]
  #indata <- indata[keepRows,]
  print(1) 
  dds_pre <- DESeqDataSetFromMatrix(countData = indata, colData = insample, design=as.formula(design))
  dds_pre <- DESeq(dds_pre)
  indata_filter = indata[rowMedians(counts(dds_pre, normalize=T)) >= 10,]
  print(2)
  #indata = indata[apply(indata,1,median)>=20,]
  
  print('indata after fitlering')
  print('before filter:')
  print(dim(indata))
  print('after filter:')
  print(dim(indata_filter))
  #print(head(indata_filter))
  
  
  if (all(colnames(indata) == rownames(insample))){
    print("Both the column names matches in position and names")
  } else{
    stop('indata and meta table does not have matching names')
  }
  
  #de_variable <- gsub("^(~\\s*)(\\w+)(.*)$", "\\2", design)
  
  #de_variable <- gsub(".*\\+(\\w+).*", "\\1", design)
  #de_variable <- gsub(".*\\+\\s*([^\\s]+)\\s*", "\\1", design)
  
  if (!grepl("\\+", design)) {
      de_variable <- strsplit(design, "~")[[1]][2]
  } else {
      de_varriable <- tail(strsplit(design, "\\+")[[1]], 1)
  } 
  
  print('de_variable')
  print(de_variable)
  insample[[de_variable]] <- as.factor(insample[[de_variable]])
  insample[[de_variable]] <- relevel(factor(insample[[de_variable]]), ref=reference)

  #if (grepl('tissue', design, fixed = TRUE) == TRUE){
  # insample$tissue = relevel(factor(insample$tissue), ref='SF')
  #}
  
  
  dds2 <- DESeqDataSetFromMatrix(countData = indata_filter, colData = insample, design= as.formula(design))
  
  print('DESeq')

  dds2 <- DESeq(dds2)
  
  print('done')
  #resultDESeq2 = results(dds2, contrast=c('tissue', 'LV', 'RV'))
  resultDESeq2 = results(dds2)
  pdf('plotMA.pdf')
  plotMA(resultDESeq2)
  dev.off()

  #dds2_vsd = vst(dds2)
  #pdf('plotPCA.pdf')
  #print(insample$type_HF)
  #plotPCA(dds2_vsd, intgroup='type_HF')
  #dev.off()
  
  res <- na.omit(resultDESeq2)

  
  write.table(res, outfile, sep='\t', quote = F)
  
}


run_DE(indata, insample, rowm, design, outfile, reference)
