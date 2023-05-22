library('EnhancedVolcano')
library(gridExtra)
library(grid)

infile = "infile_holder"

df <- read.table(infile, sep='\t', header=T)

rownames(df) <- df$hgnc

#devtools::install_github('kevinblighe/EnhancedVolcano')
df_sig <- df[which(df$padj<=0.05), ]

df_select <- df_sig[order(df_sig$log2FoldChange), ]
selected <- head(df_select, 8)

df_select <- df_sig[order(df_sig$log2FoldChange, decreasing = T), ]
selected2 <- head(df_select, 8)

selected3 <- rbind(selected, selected2)
selected4 <- df_select <- df_sig[order(df_sig$pvalue, decreasing = F), ]
selected4 <- head(selected4, 5)
selected5 <- rbind(selected4, selected3)


outfile = '.pdf'
pdf(file=outfile, width = 7, height = 8)
EnhancedVolcano(df,
                lab = rownames(df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title='',
                subtitle='',
                selectLab = selected5$hgnc,
                xlim = c(-10.5,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-2,
                FCcutoff = 1.5,
                pointSize = 1.0,
                labSize = 4.0,
                labCol = '#F8766D',
                #labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'topright',
                #legendLabSize = 12,
                #legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                col=c('#00CED1', '#00CED1', '#00CED1', '#F8766D'),
                colConnectors = '#F8766D')


dev.off()



