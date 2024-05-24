
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.19")


BiocManager::install(c("DESeq2", "edgeR", 'limma', 'statmod', 'ggplot2', 'dplyr', 'qsmooth'))
