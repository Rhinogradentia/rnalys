#DE analysis edgeR
suppressPackageStartupMessages(library('edgeR'))
suppressPackageStartupMessages(library('BiocParallel'))

# Adjust BiocParallel based on the OS
if (.Platform$OS.type == "windows") {
  register(SnowParam(4))
} else {
  register(MulticoreParam(4))
}

extract_last_variable <- function(input_string) {
  # This pattern matches variable names that consist of alphanumeric characters and underscores
  matches <- strsplit(input_string, "[^[:alnum:]_]+")[[1]]
  
  # Filter out empty strings and the tilde
  valid_matches <- matches[nchar(matches) > 0 & matches != "~"]
  
  if (length(valid_matches) == 0) {
    stop("No variable found in the design formula.")
  } else {
    # Return the last variable component from the formula
    return(tail(valid_matches, n = 1))
  }
}



run_DE <-  function (indata, insample, design, reference, outfile) {
  #insample <- read.table(insample, sep='\t', header=T)
  #rownames(insample) <- insample$id_tissue
  #insample$X <- NULL
  
  #indata <- read.table(indata, sep='\t', header=T)
  #rownames(indata) <- indata$X
  #indata$X <- NULL

  if (is.character(indata)) {
    countData <- indata <- read.table(indata, sep='\t', header=TRUE, row.names=1)
  } else if (is.matrix(indata)) {
    countData <- indata
  } else {
    stop("indata should be a matrix or a path to a CSV file containing the matrix.")
  }
  
  if (is.character(insample)) {
    insample <- read.table(insample, sep='\t', header=TRUE, row.names=1)
  } else if (is.data.frame(insample)) {
    insample <- insample
  } else {
    stop("insample should be a data frame or a path to a CSV file containing the sample information.")
  }


  plotfolder = 'data/generated/plots/'
  prefix = strsplit(outfile,'/')[[1]][4]
  prefix = strsplit(prefix, '.tab')
  
  var = extract_last_variable(design)
  group <- factor(insample[[var]]) 
  
  y1= DGEList(counts=indata, genes=rownames(indata), group = group)
  #print(y1)
  #Filter out genenes with low expression
  #keepRows <- rowSums(cpm(y1)>2) >= 2
  #indata <- indata[keepRows,]
  #print(paste("Total genes before filtering:", nrow(countData)))
  keepRows <- rowSums(cpm(y1) > 1) >= 2
  #print(paste("Genes after CPM > 1 filtering:", sum(keepRows)))
  
  # Adjust filtering here if needed based on the results of print statements
  if (sum(keepRows) == 0) {
    warning("No genes pass the filtering criteria. Adjust the filter or check data quality.")
    return(NULL)
  }
  
  y1 <- y1[keepRows, , keep.lib.sizes=FALSE]
  y1 <- calcNormFactors(y1)
  #y1 = estimateCommonDisp(y1, verbose=TRUE)
  
  if(!dir.exists(plotfolder)) {
    dir.create(plotfolder, recursive = TRUE)
  }


  pdf(paste0(plotfolder, prefix, '_MDS.pdf'))
  plotMDS(y1)
  dev.off()
  
  design <- model.matrix(as.formula(design), data=insample)
  
  #print(design)
  y1 <- estimateGLMCommonDisp(y1, design=design)
  print(paste("Genes after common dispersion:", nrow(y1$counts)))
  y1 <- estimateGLMTrendedDisp(y1, design=design)
  y1 <- estimateGLMTagwiseDisp(y1, design=design)
  y1 <- estimateDisp(y1, design=design, robust=TRUE)
  #y1 <- estimateGLMCommonDisp(y1, design=design) 


  pdf(paste0(plotfolder, prefix, '_BCV.pdf'))
  plotBCV(y1)
  dev.off()
  
  fit <- glmFit(y1, design)
  lrt <- glmLRT(fit)

  print(lrt$table)
  # Check if there are genes left to output
  if (nrow(lrt$table) > 0) {
    toptags <- topTags(lrt, n = Inf)  # Retrieve all rows
    degenes <- toptags$table
    
    pdf(paste0(plotfolder, prefix, '_DEs.pdf'))
    deGenes <- decideTests(lrt, p=0.05)
    deGenes <- rownames(lrt)[as.logical(deGenes)]
    
    write.table(degenes, outfile, sep='\t', quote = FALSE)
  } else {
    warning("No genes available after filtering to perform DE analysis. Check filtering criteria or data variability.")
  }

  #qlf <- glmQLFTest(fit)
  #print(lrt)
}

if (!interactive()) {

  args <- commandArgs(trailingOnly = TRUE)
  indata = args[1]
  insample = args[2]
  design = args[3]
  reference = args[4]
  outfile = args[5]

  run_DE(indata, insample, design, reference, outfile)
  #options(error = function() traceback(3))
}

