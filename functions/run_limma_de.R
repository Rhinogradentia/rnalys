#Run limma voom on smooth quantile normalized data
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(qsmooth))
args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
batch = args[3]
intype = args[4]
outfile = args[5]
print(indata)
print(insample)

run_DE <-  function (indata, insample, batch, intype, outfile) {
  #Indata should be data generated from qsmooth 
  
  insample <- read.table(insample, sep='\t', header=T)
  rownames(insample) <- insample$X
  insample$X <- NULL
  
  indata <- read.table(indata, sep='\t', header=T)
  rownames(indata) <- indata$X
  indata$X <- NULL
  #qs_data <- indata
  #prefix for plots
  plotfolder = 'plots/'
  prefix = strsplit(outfile, '.tab')
  #print(head(indata))
  
  if (intype=='tissue'){
    qs <- qsmooth(object = indata, group_factor = insample$tissue)
  }else{
    qs <- qsmooth(object = indata, group_factor = insample$type)
  }
  print('Performing qsmooth')
  qs_data <- qsmoothData(qs) # extract smoothed quantile normalized data
  print('done')
  #Log2 the values
  #qs_data <- log2(qs_data+1)
  d0 <- DGEList(qs_data)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  dim(d) # number of genes left
  
  #d0 <- calcNormFactors(d0)
  #SeqTag <- factor(insample$SeqTag)
  #type <- factor(insample$intype)
  if(batch == c('batch')){
    print('batch')
    if(intype == c('tissue')){
      mm <- model.matrix(~SeqTag + tissue, insample)
    } else {
      mm <- model.matrix(~SeqTag + type, insample)
    }
  } else{
    print('no batch')
    if(intype == c('tissue')){
      mm <- model.matrix(~0 + tissue, insample)
    } else {
      mm <- model.matrix(~0 + type, nsample)
    }
  }
  #mm <- model.matrix(~0 + insample$intype)
  
  y <- voom(d, mm, plot = T)
  fit <- lmFit(y, mm)
  contr <- makeContrasts(contrasts = paste(colnames(fit)[1],'-',colnames(fit)[2], sep=''), levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  print('top')
  #print(head(top.table, 20))
  print(dim(top.table[which(top.table$adj.P.Val <0.05),]))

  #Rename to make comapatable with deseq and edgeR output
  names(top.table)[5] <- 'padj'
  names(top.table)[4] <- 'pvalue'
  names(top.table)[1] <- 'log2FoldChange'
  names(top.table)[2] <- 'baseMean'
  write.table(top.table, outfile, sep='\t', quote = F)
  print('de done')
}

run_DE(indata, insample, batch, intype, outfile)
