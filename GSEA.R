library(pacman)
p_load(clusterProfiler, enrichplot, ggplot2, org.Hs.eg.db, ggnewscale, ggupset)
organism = 'org.Hs.eg.db'
library(organism, character.only = TRUE)

#Read Data from DESEq2
genelist <- results$log2FoldChange
names(genelist) <- results$hgnc
genelist <- na.omit(genelist)
genelist = sort(genelist, decreasing = TRUE)

#Gene Set Enrichment Analysis
gse <- gseGO(geneList = genelist,
             ont = "ALL",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             OrgDb = organism, 
             verbose = TRUE,
             pAdjustMethod = 'BH',
             eps = 1e-10
             )

gsekegg <- gseKEGG(geneList = kegggene,
                   keyType = "kegg",
                   exponent = 1,
                   minGSSize = 10,
                   maxGSSize = 500,
                   eps = 1e-10,
                   pvalueCutoff = 0.1,
                   pAdjustMethod = 'BH',
                   use_internal_data = FALSE
                   )

dotplot(gse, showCategory = 10, split = '.sign') + facet_grid(.~.sign)
emapplot(pairwise_termsim(gse), color = "p.adjust", min_edge = 0.2, cex_label_category = 0.8, cex_category = 0.8, cex_line = 0.8, shadowtext = TRUE,
         label_style = "shadowtext", with_edge = TRUE, nWords = 4, label_format = 30, clusterFunction = stats::kmeans)
gseaplot(gse, by = 'all', title = gse$Description[74], geneSetID = 23)

#Over Representation Analysis
ego <- enrichGO(gene = genes,
                universe = names(genelist),
                OrgDb = organism,
                keyType = "SYMBOL",
                ont = "BP",
                pvalueCutoff = 0.01,
                pAdjustMethod = 'BH',
                qvalueCutoff = 0.05,
                minGSSize = 10,
                maxGSSize = 500,
                readable = TRUE)

upsetplot(ego)
barplot(ego, drop = TRUE, showCategory = 30, font.size = 8)
emapplot(pairwise_termsim(ego), color = "p.adjust", min_edge = 0.2, cex_label_category = 0.8, cex_category = 0.8, cex_line = 0.8, shadowtext = TRUE,
         label_style = "shadowtext", with_edge = TRUE, nWords = 4, label_format = 30, clusterFunction = stats::kmeans)
goplot(ego, showCategory = 20)

##KEGG Over Represenation
ids <- bitr(names(genelist), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = organism)
deids = ids[!duplicated(ids[c('SYMBOL')]),] #Run it if there are any duplicates
results %>% filter(results$hgnc %in% deids$SYMBOL) %>% as.data.frame() -> df
df$entrez =deids$ENTREZID
kegguniverse = as.character(results$log2FoldChange)
names(kegguniverse) = df$entrez
kegguniverse = na.omit(kegguniverse)
kegguniverse = sort(kegguniverse, decreasing = TRUE)
deids %>% filter(deids$SYMBOL %in% commongenes$gene) %>% as.data.frame() -> commongenes
kegggene = as.character(commongenes$ENTREZID)

kgo <- enrichKEGG(gene = kegggene, 
                  universe = names(kegguniverse),
                  organism = 'hsa',
                  pvalueCutoff = 0.05,
                  keyType = 'kegg')
barplot(kgo, showCategory = 20, font.size = 8)
cnetplot(kgo, categorySize = 'pvalue', foldChange = NULL)
