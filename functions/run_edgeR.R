#DE analysis edgeR
suppressPackageStartupMessages(library('edgeR'))

args <- commandArgs(trailingOnly = TRUE)
indata = args[1]
insample = args[2]
design = args[3]
reference = args[4]
outfile = args[5]

options(error = function() traceback(3))
extract_variable <- function(input_string) {
  match <- regmatches(input_string, regexpr("(?<=\\+|~)\\w+$", input_string, perl = TRUE))
  
  if (match == "") {
    return(NULL)
  } else {
    return(match)
  }
}


extract_first_variable <- function(input_string) {
  match <- regmatches(input_string, regexpr("(?<=~)\\w+(?=\\+)", input_string, perl = TRUE))
  
  if (match == "") {
    return(NULL)
  } else {
    return(match)
  }
}

run_DE <-  function (indata, insample, design, reference, outfile) {
  insample <- read.table(insample, sep='\t', header=T)
  rownames(insample) <- insample$id_tissue
  insample$X <- NULL
  
  indata <- read.table(indata, sep='\t', header=T)
  rownames(indata) <- indata$X
  indata$X <- NULL

  #print(dim(indata))
  #print( dim(insample))
  #print(outfile)
  plotfolder = 'data/generated/plots/'
  prefix = strsplit(outfile,'/')[[1]][4]
  prefix = strsplit(prefix, '.tab')
   
  #if (intype == 'tissue'){
  #  y1= DGEList(counts=indata, genes=rownames(indata), group = factor(insample$tissue))
  #}else{
  #  y1= DGEList(counts=indata, genes=rownames(indata), group = factor(insample$type))
  #}

  #de_variable <- extract_variable(design)
  #print(design) 
  var = extract_variable(design)

  group <- factor(insample[[var]]) 

  #print(var)
  #print(insample[[var]])
  
  y1= DGEList(counts=indata, genes=rownames(indata), group = group)
  
  print(dim(indata)) 
  #Filter out genenes with low expression
  keepRows <- rowSums(cpm(y1)>2) >= 2
  #keepRows = rowMeans(indata) > 10
  indata <- indata[keepRows,]
  y1 <- y1[keepRows, , keep.lib.sizes=FALSE]
  y1 <- calcNormFactors(y1)
  #y1 = estimateCommonDisp(y1, verbose=TRUE)
  
  if(!dir.exists(plotfolder)) {
    dir.create(plotfolder)
  }

  pdf(paste(plotfolder, prefix,'_MDS.pdf', sep=''))
  plotMDS(y1)
  dev.off()
  
  #print(y1)

  #print(y1$samples)
  #group <- factor(insample$tissue)
  
  #design <- model.matrix(~group)

  #print(
  design <- model.matrix(as.formula(design), data=insample)
  
  
  #print(design)
  y1 <- estimateGLMCommonDisp(y1, design=design)
  y1 <- estimateGLMTrendedDisp(y1, design=design)
  y1 <- estimateGLMTagwiseDisp(y1, design=design)
  y1 <- estimateDisp(y1, design=design, robust=TRUE)
  y1 <- estimateGLMCommonDisp(y1, design=design) 
  pdf(paste(plotfolder, prefix, '_BCV.pdf', sep=''))
  plotBCV(y1)
  dev.off()
  
  fit <- glmFit(y1, design)
  lrt <- glmLRT(fit)
  #qlf <- glmQLFTest(fit)
  toptags <- topTags(lrt, n=length(rownames(indata)))
  degenes <- toptags$table
  #Quasi
  #degenes <- qlf$table
  pdf(paste(plotfolder,prefix,'_DEs.pdf', sep=''))
  deGenes <- decideTestsDGE(lrt, p=0.05)
  deGenes <- rownames(lrt)[as.logical(deGenes)]
  plotSmear(lrt, de.tags=deGenes)
  abline(h=c(-1, 1), col=2)
  
  dev.off()
    
  write.table(degenes, outfile, sep='\t', quote = F)
}


#indata = './data/generated/SLL62_LV_SLL64_RV_SLL53_RV_LV_RV_@P11701_@HWI-D0048_Normal_HFPEF_HFrEF_edgeR_counts.tab'
##insample = './data/generated/SLL62_LV_SLL64_RV_SLL53_RV_LV_RV_@P11701_@HWI-D0048_Normal_HFPEF_HFrEF_edgeR_meta.tab'
#batch = 'batch'
#outfile = 'test'
run_DE(indata, insample, design, reference, outfile)
