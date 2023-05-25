install.packages("ggplot2")
install.packages("dplyr")
#install.packages("DESeq2")
#install.packages("edgeR")
#install.packages("limma")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")


BiocManager::install(c("DESeq2", "edgeR", 'limma'))
