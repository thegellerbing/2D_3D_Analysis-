library(pacman)
p_load(dplyr,tidyverse, DESeq2, data.table,survival, survminer, readxl, magrittr, biomaRt, writexl, apeglm, ggplot2, vsn, EnhancedVolcano, ggbeeswarm, RColorBrewer)

#Read Count Data
cd = read_excel("rawcount.xlsx") %>%
  as.data.frame()
cd[is.na(cd)] = 0
cd %<>% mutate_if(is.numeric, round, 0)
pcg = read_excel('protein_coding.xlsx') %>% as.data.frame()
cd %<>% filter(cd$gene_id %in% pcg$id)

#Read Metadata
md = read_excel('metadata.xlsx') %>%
  as.data.frame()
md[is.na(md)] = 0
rm(pcg)

#DDS 
dds = DESeqDataSetFromMatrix(countData = cd,
                             colData = md,
                             design = ~condition, tidy = TRUE)
as.data.frame(colData(dds))
dds$condition <- relevel(dds$condition, "Normal Brain Sample") 

#Running DESEQ ####
rundds = DESeq(dds)
res = results(rundds) 
head(res, tidy = TRUE)
summary(res)
res %>% as.data.frame() -> results

#Adding Gene Names ####
results$ensembl <- sapply(strsplit(rownames(results),split="\\+"), "[",1)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = c("ensembl_gene_id"),
                 values = results$ensembl,
                 mart = ensembl)
results %<>% filter(results$ensembl %in% genemap$ensembl_gene_id)
genemap %<>% filter(genemap$ensembl_gene_id %in% results$ensembl)
colnames(genemap)[1] <- 'ensembl'
results <- merge(results, genemap, by = 'ensembl')

#Filtering Results###
results = results[which(results$padj < 0.05 ), ] 
resUp = results[results$log2FoldChange > 1, ] 
resDown = results[results$log2FoldChange < -1, ] 
sigresults = rbind(resUp, resDown)
rm(resUp, resDown)
list = c('ITGA3','PROM1','VIM','CD44','CDH1','EPHA3','ABCB1','MMP2','ITGA6','HIF1A','MMP9','IKBKB','PLAT','RELB','ABCA2','DKK1','CCND1',
         'NANOG','VEGFA','FN1','CDC20','EPCAM','TWIST1','SNAI1','SLC17A3','GFAP','NES','NOTCH2','MSI1','EDNRB','MYC','YAP1','CYP1A1','FZD7','MAML1','ABCA1','FOS','CDH2','SOX2','MMP1')
sigresults %>% filter(sigresults$hgnc_symbol %in% list) %>% as.data.frame() -> commongenes
write_xlsx(sigresults, './deseq_sig.xlsx')
write_xlsx(commongenes, './commongenes_deseq.xlsx')
rm(list)
save.image("./DDS.RData")

#MA Plot ####
resLFC = lfcShrink(rundds, coef = "idh_Mutant_vs_WT", type = "apeglm")
resLFC$ensembl <- sapply(strsplit(rownames(resLFC),split="\\+"), "[",1)
idx2 <- match(resLFC$ensembl, genemap$ensembl_gene_id)
resLFC$hgnc_symbol <- genemap$hgnc_symbol[idx2]
plotMA(res,ylim = c(-2,2))
plotMA(resLFC, ylim = c(-15,15),
       xlab = 'Mean of Normalised Counts',
       ylab = 'Log2 Fold Change')

#Volcano Plot ####
EnhancedVolcano(results, 
                lab = results$hgnc,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'GBM vs Normal Brain Samples',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.0,
                subtitle = NULL,
                xlim = c(10,-10),
                legendLabels = c('Insignificant', 'LFC > 1', 'p < 0.05', 'p < 0.05 & LFC > 1')
                ) 

#Plotting Gene Counts ####
par(mfrow=c(2,3)) # of plots per pic
plotCounts(rundds,gene="ENSG00000026508",intgroup="condition")
gplot <- plotCounts(rundds, gene="ENSG00000005884", intgroup = "condition", returnData = TRUE)
ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
  geom_boxplot(width = 0.7) +
  geom_beeswarm(cex = 0.7, size = 1.2) +
  scale_y_log10(labels = scales::comma) +
  xlab("Tissue Type")+
  ylab("Normalised Counts")+
  scale_color_brewer(palette = "Set1") +
  ggtitle("ITGA3") +
  theme_pubclean() + 
  theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = 'none',
        axis.line= element_line(colour='black', size = 0.8, linetype = 'solid'),
        axis.title = element_text(size = 16), axis.text = element_text(size = 15)) 

####Important: There are two lists that you need to first key in before running the loop (genelist & list). genelist has to be the Ensembl ID so that the loop knows what gene to pull for the plot, and the list has to be the name of the gene so that the loop knows what name to plot on the graph.

for (i in unique(1:length(genelist)))
{
  gp <- genelist[i]
  gplot <- plotCounts(rundds, gene= gp , intgroup = "condition", returnData = TRUE)
  ggraph = 
    ggplot(gplot, aes(x = condition, y = count, colour = condition)) +
    geom_boxplot(width = 0.5) +
    geom_beeswarm(cex = 0.9, size = 1.3) +
    scale_y_log10(labels = scales::comma) +
    xlab("Tissue Type")+
    ylab("Normalised Counts")+
    scale_color_brewer(palette = "Set1") +
    ggtitle(list[i]) +
    theme_pubclean() + 
    theme(plot.title = element_text(hjust = 0.5, size = 18), legend.position = 'none',
          axis.line= element_line(colour='black', size = 0.8, linetype = 'solid'),
          axis.title = element_text(size = 16), axis.text = element_text(size = 15))
  fn = list[i]
  ggsave(ggraph, filename = paste(fn, ".png", sep = ""))
}

#Data Transformations & Visualisation ####
vsd = vst(rundds, blind=FALSE) 
head(assay(vsd), 3)
ntd = normTransform(rundds, f = log2, pc = 1)
meanSdPlot(assay(ntd))
rld = rlog(rundds, blind = FALSE)

#Heatmap of Count Matrix ####
library(pheatmap)
select <- order(rowMeans(counts(rundds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(rundds)[,"condition"])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#Heatmap of Sample-to-Sample Distances ####
sampleDists = dist(t(assay(vsd)))
samplematrix = as.matrix(sampleDists)
rownames(samplematrix) <- paste(vsd$condition, sep="-")
colnames(samplematrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9,"Blues"))) (255)
pheatmap(samplematrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#Expression Heatmap of Specific Genes ####
sample = assay(vsd)[genelist,]
sample <- (sample - rowMeans(sample))/rowSds(sample)
collabel <- as.data.frame(colData(vsd)[c('condition')])
colnames(collabel)[1] <- 'Condition'
rowlabel <- commongenes %>% dplyr::select(ensembl, hgnc_symbol)
rowlabel <- rowlabel %>% remove_rownames() %>% column_to_rownames(var = 'ensembl')
colour <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
paletteLength <- 50
myBreaks <- c(seq(min(sample), 0, length.out=ceiling(paletteLength/2) + 1), 
             seq(max(sample)/paletteLength, max(sample), length.out=floor(paletteLength/2)))
pheatmap(sample, annotation_col = collabel, show_colnames = FALSE, cluster_cols = FALSE, labels_row = rowlabel$hgnc_symbol, color = colour, breaks = myBreaks)

#PCA Plot
plotPCA(vsd, intgroup = "tissue_type")

#Filtering Data Set ####
keep = rowSums(counts(rundds)) > 1 
rundds <- rundds[keep,] 

#Exporting Data Frame ####
write.csv(res,"C:/Users/brand/Documents/R/TCGA\\results.csv",row.names=TRUE)

#Checking for NA ####
any(is.na(cd))

#Converting DDS Results into Dataframe ####
dfres = as.data.frame(res)