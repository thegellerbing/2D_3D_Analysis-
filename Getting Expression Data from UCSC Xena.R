library(pacman)
p_load(UCSCXenaTools, data.table, R.utils, dplyr, magrittr)

###Generate an object tracking all datasets from Xena Data Hubs ###
data(XenaData)
write.csv(XenaData, "00_tblXenaHubInfo.csv")

### Target = TCGA Clinical Data
cohort = "TCGA Glioblastoma"
dataset = "TCGA.GBM.sampleMap/GBM_clinicalMatrix"
clin_TCGA = XenaGenerate(subset = XenaHostNames == 'tcgaHub') %>% 
  XenaFilter(filterCohorts = cohort) %>% 
  XenaFilter(filterDatasets = dataset)
XenaQuery(clin_TCGA) %>% 
  XenaDownload(destdir = "./")

### Target = TCGA Survival Data
surv_TCGA = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TCGA_survival_data")
XenaQuery(surv_TCGA) %>%
  XenaDownload(destdir = "./")

### Target = GTEx Phenotype Data
pheno_GTEX = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype")
XenaQuery(pheno_GTEX) %>%
  XenaDownload(destdir = "./")

###Retrieve IDs for GTEx normal samples of desired tissue type
filterGTEX = fread("TcgaTargetGTEX_phenotype.txt.gz")
names(filterGTEX) = gsub("\\_", "", names(filterGTEX))
GTEXstudy = "GTEX"
GTEXsite = "Brain"
GTEXtissue = "Brain - Cortex"
GTEXtissue2 = "Brain - Frontal Cortex \\(Ba9\\)"
filterGTEX2 = subset(filterGTEX, study == GTEXstudy & primarysite == GTEXsite & grepl(GTEXtissue, filterGTEX$`primary disease or tissue`)) 
filterGTEX3 = subset(filterGTEX, study == GTEXstudy & primarysite == GTEXsite & grepl(GTEXtissue2, filterGTEX$`primary disease or tissue`)) 

### Retrieve IDs for TCGA primary tumour samples of desired histological type
filterTCGA = fread(dataset)
names(filterTCGA) = gsub("\\_", "", names(filterTCGA))

####Subset expression data to include only desired samples
filterexpr = c(GTEX_metadata$sample, filterTCGA$sampleID)
filterensembl = 'sample'
exprsubsetbysamp = fread("TcgaTargetGtex_gene_expected_count.csv", select = c(filterensembl,filterexpr))

###Subset expression data to include only protein coding genes
probemap = fread('probeMap_gencode.v23.annotation.gene.probemap')
expression = merge(probemap, exprsubsetbysamp, by.x = 'id', by.y = 'sample')
pcgenes = read_excel('protein_coding_genes.xlsx') %>% as.data.table()
finalexpression = subset(expression, gene%in% pcgenes$Gene_Symbol)

###Subset Desired Clinical Data
clinvar = c('sampleID', 'CDEsurvivaltime', 'sampletype', 'CDEvitalstatus', 'PATIENT', 'GeneExpSubtype', 'CDEtherapy', 'GeneExpSubtype')
clindf = as.data.frame(do.call(cbind, filterTCGA))
md = clindf[clinvar]
