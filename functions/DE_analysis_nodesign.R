DE_analysis_nod<- function (fit, y, prefix, foldchange=1, output,pval=0.05) {
  ## fit is the glmfit
  ## contrast is contrast for sample combination
  ## y DEoject object
  ## perfix is the info paste in the output file
  ## output is the project folder of the output
  library(ggplot2)
  library(grid)
  library(edgeR)
  library(colorRamps)
  library(biomaRt)
  library(org.Hs.eg.db)
  
  # Function for coverting the gene symbol
  convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
    stopifnot( inherits( db, "AnnotationDb" ) )
    ifMultiple <- match.arg( ifMultiple )
    suppressWarnings( selRes <- AnnotationDbi::select( 
      db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
    if( ifMultiple == "putNA" ) {
      duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
      selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
    return( selRes[ match( ids, selRes[,1] ), 2 ] )
  }
  lrt <- glmLRT(fit)
 
  tt<-topTags(lrt)
  allGeneList <- topTags(lrt,n = nrow(lrt$table))
  print('TopTags:')
  #print(tt)
  print(summary(de<-decideTestsDGE(lrt,p.value = pval,lfc = foldchange)))
  
  detags<-rownames(y)[as.logical(de)]
  # highlight some of the DE genes
  gt<-lrt$table
  gt$threshold<-(rep('Nosig',times = nrow(gt)))
  gt[abs(gt$logFC) >foldchange ,5]<-'Nosig-highFC'
  gt[abs(gt$logFC) >foldchange & rownames(gt) %in% detags ,5]<-'Sig'
  gt$threshold<-factor(gt$threshold,levels=rev(unique(gt$threshold)))

  gt$GeneSymbol <- NA
  fnal<-cpm(y)[rownames(cpm(y)) %in% detags,]
  
  # Extracting DE genes
  DEgenes<-lrt$table[rownames(lrt$table) %in% detags,]
  write.table(DEgenes, outfile, sep='\t', quote = F)

  #}
}



DE_analysis_lncRNA<- function (fit, coef, contrast,y, prefix, foldchange=1, output,pval=0.05) {
  ## fit is the glmfit
  ## contrast is contrast for sample combination
  ## y DEoject object
  ## perfix is the info paste in the output file
  ## output is the project folder of the output
  library(ggplot2)
  library(grid)
  library(edgeR)
  library(colorRamps)
  library(biomaRt)
  library(org.Hs.eg.db)
  
  
  convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
    stopifnot( inherits( db, "AnnotationDb" ) )
    ifMultiple <- match.arg( ifMultiple )
    suppressWarnings( selRes <- AnnotationDbi::select( 
      db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
    if( ifMultiple == "putNA" ) {
      duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
      selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
    return( selRes[ match( ids, selRes[,1] ), 2 ] )
  }
  
  
  if(missing(contrast)){
    lrt<-glmLRT(glmfit = fit,coef = coef)
  }
  else{
    lrt<-glmLRT(glmfit = fit,contrast = contrast)  
  }
  tt<-topTags(lrt)
  print('TopTags:')
  #print(tt)
  
  print(summary(de<-decideTestsDGE(lrt,p.value = pval)))
  
  detags<-rownames(y)[as.logical(de)]

  
  
  
  # Extracting DE genes
  DEgenes<-lrt$table[rownames(lrt$table) %in% detags,]
  cat("DEgenes:", nrow(DEgenes))
  cat("DETag:", length(detags))
  
  
  #if(nrow(DEgenes) > 1){
  #  DEgenes$GeneSymbol <- NA
  #  #DEgenes$GeneSymbol <- convertIDs( rownames(DEgenes), "ENSEMBL", "SYMBOL", org.Hs.eg.db,ifMultiple = "putNA")
  #  
  #  write.table(DEgenes,file = file.path(output,paste("DEgenes",paste(prefix,'txt',sep = '.'),sep = '_')),quote = F,sep = '\t',row.names = T)
  #  #head(DEgenes)
  #  DEgenesFiltered<-DEgenes[which(abs(DEgenes$logFC) > foldchange) ,]
  #  cat('Number of Total DE Genes:', nrow(DEgenesFiltered), '\n')
  ##  
  #  upRegGenes<-DEgenesFiltered[which(DEgenesFiltered$logFC > foldchange ),]
  #  dnRegGenes<-DEgenesFiltered[which(DEgenesFiltered$logFC < foldchange ),]
  #  cat('Number of Upregulated Genes:', nrow(upRegGenes), '\n')
  #  cat('Number of Downregulated Genes:', nrow(dnRegGenes), '\n')
  #  
  #  # Fetching gene description
  #  # Getting gene info for all the genomes
  #  #ensembl_hs<-useMart('ensembl',dataset = 'hsapiens_gene_ensembl',host = "www.ensembl.org")
  #  #if(nrow(upRegGenes) > 1){
  #  #  #upRegGenes_info<-getBM(attributes = c('ensembl_gene_id','wikigene_description','hgnc_symbol'), filters = 'ensembl_gene_id',values = rownames(upRegGenes), mart = ensembl_hs)
  #  #  #upRegGenes_final<-merge(upRegGenes_info,upRegGenes,by.x = 'ensembl_gene_id',by.y = "row.names",all.y = T)
  #  #  
  #  #  upGeneFile<-paste('Upregulated_Genes',prefix,sep = '_')
  #  #  write.table(upRegGenes,file = file.path(output,paste(upGeneFile,'txt',sep = '.')),quote = F,sep = '\t',row.names = F)
  #  #}
    
  #  #if(nrow(dnRegGenes) > 1){
  #  #  #dnRegGenes_info<-getBM(attributes = c('ensembl_gene_id','wikigene_description','hgnc_symbol'), filters = 'ensembl_gene_id',values = rownames(dnRegGenes), mart = ensembl_hs)
  #  #  
  #  #  #dnRegGenes_final<-merge(dnRegGenes_info,dnRegGenes,by.x = 'ensembl_gene_id',by.y = "row.names",all.y = T)
  #  #  dnGeneFile<-paste('Downregulated_Genes',prefix,sep = '_')
  #  #  #cat(output)
  #  #  write.table(dnRegGenes,file = file.path(output,paste(dnGeneFile,'txt',sep = '.')),quote = F,sep = '\t',row.names = F)
  #    
  #  #}
  #  #upRgenesAllSamples<-fnal[which(rownames(fnal) %in% rownames(upRegGenes)),]
  #  #dnRgenesAllSamples<-fnal[which(rownames(fnal) %in% rownames(dnRegGenes)),]
    return(DEgenesFiltered)
  #}
}