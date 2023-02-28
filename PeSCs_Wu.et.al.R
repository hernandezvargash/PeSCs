
# main --------------------------------------------------------------------

# November 2022

# Analyses for Cell Report manuscript (Wu et al.):

# Identification of a CD106+ pericyte stem cell leading to Ly6G+ cell accumulation
# responsible for resistance to immunotherapy in pancreatic cancer



# reviewer's comments -----------------------------------------------------

# specific reviewer's comments:

# 7) Figure 2 was a key figure to show the PeSC in single-cell RNA-seq data. 
# However, this result was only superficially presented without in-depth analyses or adequate data presentation to support the main conclusion.
# As mentioned by Reviewer 2 Q1 in Fig. 2A, the surface markers used for PanIN-enriched population (containing CD106+ PeSCs) sorting are 
# DAPI-CD45-CD31-Lectin PNA-EpCAM+. 
# It is understandable that after EpCAM sorting the total cell mixture can still have EpCAM-low clusters such as the Cluster7. 
# However, it is still confusing whether the cells in Figure 2B refer to only the DAPI-CD45-CD31-Lectin PNA-EpCAM+ 
# (the lower panel of Figure 2A, PanIN cell enriched fraction), 
# or from a combination of both fibroblast-enriched fraction (the upper panel of Figure 2A) and PanIN cell enriched fraction. 
# Figure 2C heatmap is just the top genes from the cell clusters automatically defined by algorithm, without showing the key genes 
# such as epithelial cell genes (EpCAM, Krt8, Krt18, Krt19), fibroblast genes (Pdgfra, Col1a1, Dcn, Pdpn), 
# pericyte genes (Cspg4, Rgs5, Pdgfrb), or stem cell genes (Cd24a, Cd44).
# Only a few of these genes were shown in supplementary figure 4 as UMAP format but not shown in heatmap or violin plot format.

# 8) Note also that the cell clustering of Figure 2 was merely defined by the algorithm, 
# without biological relevance (defining which cell subcluster is what cell subtype). 
# Figure 2(B-C) may include an error with cell cluster definition as the algorithm defines 10 total cell clusters (from cluster 0 to cluster 9). 
# But the legend defined the 10 colored dots into 9 clusters (1-9).

# 9) The author mentioned that the single-cell RNA-seq was done on 5 KC mice with a total cell number of 2000. 
# The total cell number is too low to support any conclusion on the identification of Cluster 7 PeSC (which has perhaps only 100 cells or so). 
# Such low number of cells for Cluster 7 (as well as all other cell types combined) will cause the issue of inaccuracy in cell clustering. 
# Based on in supplementary figure 4, the potential PeSC cluster (Cluster 7) has very low expression of Cspg4 and Rgs5, 
# with certain expression of Pdgfrb, and with robust expression of Col1a1/2, Col3a1, Cd34, Cxcl12. 
# It seems that the Cluster 7 is more fibroblast-like, than pericyte-like.

# 10) The analysis on human PDAC single-cell RNA-seq data does not appear to support the conclusion regarding the presence of PeSC cluster. 
# The usage of PeSC score has no specific relevance to the so-called PeSC, 
# but rather reflects a combination of mesenchymal genes associated with all fibroblasts and pericytes combined. 
# This is likely why Figure 3H right panel showed the high PeSC score in the entire fibroblast cluster. 
# The analysis of this human PDAC single-cell RNA-seq data did not show supporting data for the cell clustering or identification 
# (only superficially shown in supplementary figure 5). 
# Together with the point raised above, the identification of the PeSC in PDAC is not robustly validated, 
# and appears disconnected from the remaining functional studies using isolated cells, the identiy of which remain obscure.



# Strategy ----------------------------------------------------------------


# add a published list of Perycyte genes. Two candidates: GSE17014 and https://pubmed.ncbi.nlm.nih.gov/34018679/

# add a new larger human dataset: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9283889/


# Specifically

# point #7 & #8:

## better description of the analysis and cluster/gene selection
## new heatmap with the relevant markers that she suggests will help in the interpretation
## a violin and a dotplot with the same genes so we can choose the best representation
## check and correct the cluster labels.
## add new PeSC scores

# point #9:
## check the exact number of cells, but acknowledge a low number of cells cannot be conclusive and agree with the reviewer
## use an extra list of pericyte-specific genes to rule out that these are fibroblasts (new PeSC scores)

# point #10:
## reanalyze PDAC with new PeSC scores
## add new Atlas dataset: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9283889/
## subcluster fibroblasts and re-cluster to more precisely select the putative PeSCs.
## same thing with stellate cells.
## compare tumor vs normal.



# libraries ---------------------------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(Nebulosa)
  library(xlsx)
  library(viridis)
  library(oligo)
  library(stringr)
  library(biomaRt)
  library(limma)
  
})

set.seed(444)



# GSE17014 ----------------------------------------------------------------

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17014

# human dermal pericyte microarray expression originally published here: https://www.jci.org/articles/view/38535#sd

## Table S1 -----

# a shorter (146 genes) genelist was described here: https://pubmed.ncbi.nlm.nih.gov/26589433/
# after crossing the original data with ovarian CAFs (Table S1 below)

T.S1 <- read.csv("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/GSE17014_Table.S1.csv", row.names = 1)
head(T.S1)

table(T.S1$fc > 0) # 635
table(T.S1$p.value < 0.05) # 196

T.S1 <- T.S1[T.S1$p.value < 0.05, ]
T.S1 <- T.S1[T.S1$fc > 0, ]

length(unique(T.S1$NAME)) # 104
table(!is.na(T.S1$NAME)) # 149

# the list doesn't add to 146 after trying different filters
# not clear where the 146 figure is coming from

## Raw data -----

cel <- read.celfiles(list.files("GSE17014_RAW", full.names = T))
cel
sampleNames(cel) <- str_sub(sampleNames(cel), 1, 9)

eset <- rma(cel)
data <- exprs(eset)
dim(data) # 54675 genes and 8 samples
par(mar=c(10,5,5,5), mfrow = c(1,1))
boxplot(log10(data), las=2, main = "GSE17014")

pdata <- pData(eset)
pdata$group <- c(rep("bri", 4), rep("dim", 4)) # bri = dermal pericytes, dim = remaining dermal cells
pData(eset) <- pdata

#featureData(eset) <- getNetAffx(eset,"transcript")

ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl) 
keytypes(ensembl)
anno <- getBM(attributes = c('affy_hg_u133_plus_2','hgnc_symbol','entrezgene_id'),
              filters = 'affy_hg_u133_plus_2',
              values = rownames(eset),
              mart=ensembl)
head(anno)

anno.filt <- anno[!anno$hgnc_symbol == "", ]
anno.filt <- anno.filt[!duplicated(anno.filt$affy_hg_u133_plus_2), ] # this keeps one value for each duplicated probe ID
rownames(anno.filt) <- anno.filt$affy_hg_u133_plus_2
head(anno.filt)

common.probes <- intersect(rownames(eset), rownames(anno.filt))
anno.filt <- anno.filt[common.probes, ]

eset <- eset[common.probes, ]

fdata <- featureData(eset)
head(fdata@data)
fdata@data$EntrezID <- anno.filt$entrezgene_id
fdata@data$Symbol <- anno.filt$hgnc_symbol
featureData(eset) <- fdata


save(eset, file = "eset_GSE17014.RData")


## DEGs -----

load("eset_GSE17014.RData")

pData <- pData(eset)
exprs<-exprs(eset)

mod<-model.matrix(~0+as.factor(group),data=pData)
fit = lmFit(eset,mod)
colnames(fit)
contrast.matrix<-c(1,-1)
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
top<-topTable(eb, adjust="fdr", p=0.05, number=10000, lfc=1)
head(top) # 2197

hist(top$P.Value)
table(top$logFC > 0) # 888 

write.csv(top, file="DEGs_bri.vs.dim_GSE17014.csv")


intersect(top$Symbol, T.S1$NAME)
setdiff(top$Symbol, T.S1$NAME)
setdiff(T.S1$NAME, top$Symbol)



# other gene lists --------------------------------------------------------


# mouse brain specific pericyte expressed genes were downloaded from:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763650/ (Table S1)


pegs <- read.csv("mouse.brain.PEGs.csv", row.names = 1)
head(pegs)

intersect(top$Symbol, toupper(pegs$Gene))


# the full dataset was downloaded from: http://betsholtzlab.org/VascularSingleCells/database.html
# and the same group also published single cell data on lung and brain vascular cells
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6103262/
# check als: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7414220/ from the same group
# and its linked data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150294


# this one is quite relevant as for pancreatic pericytes:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7664344/
# they express ACE2 and are therefore sensitive to SARS-Cov2 infection
# https://pubmed.ncbi.nlm.nih.gov/33281748/
# also https://pubmed.ncbi.nlm.nih.gov/35452595/
# refer to this PANC-DB database: https://hpap.pmacs.upenn.edu/

# Pericytes in the CNS, heart, and pancreas express ACE2 strongly, 
# as do perineurial and adrenal fibroblasts, 
# whereas endothelial cells do not at any location analyzed.

# https://www.sciencedirect.com/science/article/pii/S0012160618307942?via%3Dihub
# Nkx3.2 TF activity as an indicator of pericytes?
# while Nkx3.2 is expressed in early embryonic pancreatic tissues, its expression declines at later embryonic ages and is absent from the adult tissue
# We focused on nine genes previously shown to expressed by pericytes (Armulik et al., 2011):
# Abcc9, Anpep, Cd248, Cspg4, Dlk1, Eng, Kcnj8, Pdgfrb, and Rgs5
# this dataset should be useful for pericyte transcriptome using the later time points:
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7448
# To determine changes in pancreatic mesenchymal cells during embryonic development, 
# we analyzed the transcriptome of these cells at three embryonic days: 13.5 ('L1'), 15.5 ('L2'), and 17.5 ('L3'). 
# Pancreatic mesenchymal cells were isolated by FACS from Nkx3.2-Cre;LSL-YFP pancreatic tissues based on these cells fluorescent labeling. 
# Two groups of mice were analyzed for each embryonic day ('P1', 'P2').



# Ductal dataset ----------------------------------------------------------

rm(list=ls())

load(file="~/Dropbox/BioInfo/Colabs/PeSCs/D.RData")

DimPlot(D) # 886 cells
DimPlot(D, pt.size = NULL, label = T, label.box = T, label.size = 7)

table(Idents(D))
# 0   1   2   3   4   5   6   7   8   9 
# 161 149 137 106  83  77  62  48  39  24 


head(D[])

f1 <- FeaturePlot(D, features = "PeSC.score1",cols = viridis(10), order = T) + ggtitle("PeSC Score") + ggtitle("PeSC Score from curated genelist")

sc.all <- D

#DefaultAssay(sc.all) <- "SCT"

sel.genes <- c("PDGFRB","ACE2","ABCC9","ANPEP","CD248",
               "CSPG4","DLK1","ENG","KCNJ8","RGS5",
               "PDGFRA","FN1","COL1A1","DCN","PDPN",
               "CD34", "CD24A","CD44",
               "DES","NES","ACTA2","NKX3-2", 
               "TGFBI","VEGFA","NGF","ACVR1B",
               "CALD1","MCAM","ANGPT1","ANGPT2",
               "EPCAM","KRT8","KRT18","KRT19"
               )
sel.genes <- str_to_title(sel.genes)
setdiff(sel.genes, rownames(sc.all))
int.1 <- intersect(sel.genes, rownames(sc.all))
jpeg("pericyte.markers_sc.all_feature.plot.jpg", width = 1500, height = 1500, quality = 100)
FeaturePlot(sc.all, sel.genes, order = T,  cols = viridis(10))
dev.off()
jpeg("pericyte.markers_sc.all_density.plot.jpg", width = 1500, height = 1000, quality = 100)
plot_density(sc.all, int.1, joint = F, pal = "inferno")
dev.off()

for(i in 1:length(sel.genes)){
  
  jpeg(paste0("cluster7_feature.plot_", sel.genes[i], ".jpg"), 
       width = 600, height = 500, quality = 100)
  p1 <- FeaturePlot(sc.all, sel.genes[i], cols = viridis(10), order = T, pt.size = 1)
  plot(p1)
  dev.off()
  
}

# plot density changes the UMAP distances ?

for(i in 1:length(int.1)){
  
  jpeg(paste0("cluster7_density.plot_", int.1[i], ".jpg"), 
       width = 600, height = 500, quality = 100)
  p1 <- plot_density(sc.all, int.1[i], joint = F, pal = "inferno")
  plot(p1)
  dev.off()
  
}


DotPlot(D, features = int.1,
        assay = "SCT",
        cols = c("blue", "red"),) + RotatedAxis()



## PeSC score Anca-curated  ---------


# list of PeSC markers curated by Anca using literature and enriched with Loupe analysis
genes <- read.xlsx("~/Dropbox/BioInfo/Colabs/PeSCs/Upregulated genes.xlsx", sheetIndex = 1)
head(genes)
genes <- unique(str_to_title(genes$gene_name))

#DefaultAssay(sc.all) <- "SCT"

sc.all <- AddModuleScore(sc.all, features = list(genes), name = "PeSC.score", ctrl = 100)
# The following features are not present in the object: ACVR1B, NANOG, NA, SAA3, NKX3-2, ABCC9, CSPG4, DLK1, NGF, FOXD1, CD274, COL4A6, COL6A4, TGFB1, PDGFB, CTGF
FeaturePlot(sc.all, features = "PeSC.score1",cols = viridis(10), order = T) + ggtitle("PeSC Score") + ggtitle("PeSC Score from curated genelist")

## PeSC score Loupe  ---------

# 2nd list (based on cluster 6 from loupe -> Anca)

genes.2 <- read.csv("~/Dropbox/BioInfo/Colabs/PeSCs/Cluster 6.csv", row.names = 1)
head(genes.2)
genes.2 <- unique(str_to_title(genes.2$FeatureName))

sc.all <- AddModuleScore(sc.all, features = list(genes.2), name = "PeSC.score2.", ctrl = 100)
FeaturePlot(sc.all, features = "PeSC.score2.1",cols = viridis(10), order = T) +
  ggtitle("PeSC Score from Cluster.6 scRNAseq")



## PeSC score 146 genelist  ---------

# a shorter (146 genes) genelist was described here: https://pubmed.ncbi.nlm.nih.gov/26589433/
# after crossing the original data with ovarian CAFs (Table S1 below)

T.S1 <- read.csv("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/GSE17014_Table.S1.csv", row.names = 1)
head(T.S1)

table(T.S1$fc > 0) # 635
table(T.S1$p.value < 0.05) # 196

T.S1 <- T.S1[T.S1$p.value < 0.05, ]
T.S1 <- T.S1[T.S1$fc > 0, ]

T.S1 <- unique(str_to_title(T.S1$NAME)) # 104

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.score.TS1.", ctrl = 100)
FeaturePlot(sc.all, features = "PeSC.score.TS1.1", cols = viridis(10), order = T) + 
  ggtitle("PeSC Score from Table S1 ovarian cancer paper")



## PeSC score dermal pericytes  ---------

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17014

T.S1 <- read.csv("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/DEGs_bri.vs.dim_GSE17014.csv", row.names = 1)
head(T.S1)

table(T.S1$logFC > 0) # 888
table(T.S1$adj.P.Val < 0.05) # 2197
T.S1 <- T.S1[T.S1$logFC > 0, ]
T.S1 <- unique(str_to_title(T.S1$Symbol)) # 606

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.dermal.", ctrl = 100)
FeaturePlot(sc.all, features = "PeSC.dermal.1",cols = viridis(10), order = T) +
  ggtitle("Dermal Pericyte Score (GSE17014)")



## Brain Pericytes score  ---------

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763650/

T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/mouse.brain.PEGs.csv", sep = ",")
head(T.S1) # 200

table(T.S1$Highest.expressionb)
T.S1 <- T.S1[T.S1$Highest.expressionb == "Pericytes", ]
T.S1 <- unique(str_to_title(T.S1$Gene)) # 62  or 200

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.brain.", ctrl = 100)
FeaturePlot(sc.all, features = "PeSC.brain.1",cols = viridis(10), order = T) +
  ggtitle("Brain Pericyte Score")



## Pericyte score Nkx3-2 paper  ---------

# https://www.sciencedirect.com/science/article/pii/S0012160618307942?via%3Dihub

# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7448
T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L2vsL3.genes.edgeR-tag.txt")
rownames(T.S1) <- make.names(T.S1$id, unique = T)
T.S1 <- na.omit(T.S1)
head(T.S1) # 29864 -> 15060

T.S1 <- T.S1[T.S1$EdgeR.tag.logFC > 0, ] # 7930
T.S1 <- T.S1[T.S1$edgeR.tag.adjusted.p.value < 0.05, ] # 1815
T.S1 <- unique(str_to_title(T.S1$id)) # 2197

genes.present <- intersect(T.S1, rownames(sc.all)) # 330

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.L2L3.Up", ctrl = 100)
FeaturePlot(sc.all, features = "PeSC.Nkx.L2L3.Up1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L2.vs.L3.Up")

T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L1vsL3.genes.edgeR-tag.txt")
rownames(T.S1) <- make.names(T.S1$id, unique = T)
T.S1 <- na.omit(T.S1)
head(T.S1) # 29864 -> 15060

T.S1 <- T.S1[T.S1$EdgeR.tag.logFC > 0, ] # 7930
#T.S1 <- T.S1[T.S1$EdgeR.tag.logFC < 0, ] # 7930
T.S1 <- T.S1[T.S1$edgeR.tag.adjusted.p.value < 0.05, ] # 1815
T.S1 <- unique(str_to_title(T.S1$id)) # 2197

genes.present <- intersect(T.S1, rownames(sc.all)) # 330

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.", ctrl = 100) # L1 vs L3 Up
FeaturePlot(sc.all, features = "PeSC.Nkx.1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L3.Up")

T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L1vsL2.genes.edgeR-tag.txt")
rownames(T.S1) <- make.names(T.S1$id, unique = T)
T.S1 <- na.omit(T.S1)
head(T.S1) # 29864 -> 15060

T.S1 <- T.S1[T.S1$EdgeR.tag.logFC > 0, ] # 7930
#T.S1 <- T.S1[T.S1$EdgeR.tag.logFC < 0, ] # 7930
T.S1 <- T.S1[T.S1$edgeR.tag.adjusted.p.value < 0.05, ] # 1815
T.S1 <- unique(str_to_title(T.S1$id)) # 2197

genes.present <- intersect(T.S1, rownames(sc.all)) # 330

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.L1L2.Up", ctrl = 100)
FeaturePlot(sc.all, features = "PeSC.Nkx.L1L2.Up1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L2.Up")


head(sc.all[])
saveRDS(sc.all, file = "D.rds")




# GSA-PDAC ----------------------------------------------------------------

# this dataset is included in the PDAC scRNASeq Atlas (see next sections)

# PRJCA001063 / CRA001160
# https://ngdc.cncb.ac.cn/gsa/browse/CRA001160

# single cell PDAC data corresponding to manuscripts:
# https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00108-0/fulltext
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6796938/




# PDAC single cell Atlas --------------------------------------------------

rm(list=ls())

load(file="~/Dropbox/BioInfo/freelance/StromaCare/scPADC/sc.all.RData")

head(sc.all[]) # 38919 features across 136163 samples
length(unique(sc.all$Patient)) # 72
table(sc.all$Patient)
table(sc.all$Type)
# Normal  Tumor 
# 22926 113237 
table(sc.all$Project)
# CA001063  GSE111672  GSE154778  GSE155698 GSM4293555     OUGS 
#   57423       3625       8000      50943       4874     11298 
table(Idents(sc.all))
#Fibroblast cell  Stellate cell  Macrophage cell  Endothelial cell  T cell  B cell  Ductal cell type 2  Endocrine cell  Ductal cell type 1  Acinar cell 
#     11003           7624            23703            10059        27569     5838         33607              1298            11605             3857 

Idents(sc.all)
DimPlot(sc.all, label = T, label.size = 6, repel = T, label.color = "darkblue") + NoLegend()

DimPlot(sc.all, group.by = "Project")
DimPlot(sc.all, group.by = "Patient", label = T, repel = T)+ NoLegend()
DimPlot(sc.all, group.by = "Type")

DefaultAssay(sc.all) <- "SCT"

sel.genes <- c("PDGFRB","ACE2","ABCC9","ANPEP","CD248",
               "CSPG4","DLK1","ENG","KCNJ8","RGS5",
               "DES","NES","ACTA2","NKX3-2", 
               "TGFBI","VEGFA","NGF","ACVR1B",
               "PDGFRA","FN1","COL1A1","DCN","PDPN",
               "CD34", "CD24","CD44",
               "EPCAM","KRT8","KRT18","KRT19",
               "CALD1","MCAM","ANGPT1","ANGPT2")
setdiff(sel.genes, rownames(sc.all))
int.1 <- intersect(sel.genes, rownames(sc.all))
jpeg("pericyte.markers_sc.all_feature.plot.jpg", width = 1500, height = 1500, quality = 100)
FeaturePlot(sc.all, sel.genes, order = T,  cols = viridis(10))
dev.off()
jpeg("pericyte.markers_sc.all_density.plot.jpg", width = 1500, height = 1000, quality = 100)
plot_density(sc.all, int.1, joint = F, pal = "inferno")
dev.off()

for(i in 1:length(sel.genes)){
  
  jpeg(paste0("sc.all_feature.plot_", sel.genes[i], ".jpg"), 
       width = 600, height = 500, quality = 100)
  p1 <- FeaturePlot(sc.all, sel.genes[i], cols = viridis(10), order = T, pt.size = 0.5)
  plot(p1)
  dev.off()
  
}

for(i in 1:length(int.1)){
  
  jpeg(paste0("sc.all_density.plot_", int.1[i], ".jpg"), 
       width = 600, height = 500, quality = 100)
  p1 <- plot_density(sc.all, int.1[i], joint = F, pal = "inferno")
  plot(p1)
  dev.off()
  
}



## PeSC score Anca-curated  ---------


# list of PeSC markers curated by Anca using literature and enriched with Loupe analysis
genes <- read.xlsx("~/Dropbox/BioInfo/Colabs/PeSCs/Upregulated genes.xlsx", sheetIndex = 1)
head(genes)
genes <- unique(toupper(genes$gene_name))

DefaultAssay(sc.all) <- "integrated"
#DefaultAssay(sc.all) <- "SCT" # there is no SCT assay

sc.all <- AddModuleScore(sc.all, features = list(genes), name = "PeSC.score", ctrl = 100)
# The following features are not present in the object: ACVR1B, NANOG, NA, SAA3, NKX3-2, ABCC9, CSPG4, DLK1, NGF, FOXD1, CD274, COL4A6, COL6A4, TGFB1, PDGFB, CTGF
f1 <- FeaturePlot(sc.all, features = "PeSC.score1",cols = viridis(10), order = T) + ggtitle("PeSC Score") + ggtitle("PeSC Score from curated genelist")
f2 <- plot_density(sc.all, features = "PeSC.score1", joint = F, pal = "inferno") + ggtitle("PeSC Score from curated genelist")
f1+f2

## PeSC score Loupe  ---------

# 2nd list (based on cluster 6 from loupe -> Anca)

genes.2 <- read.csv("~/Dropbox/BioInfo/Colabs/PeSCs/Cluster 6.csv", row.names = 1)
head(genes.2)
genes.2 <- unique(toupper(genes.2$FeatureName))

sc.all <- AddModuleScore(sc.all, features = list(genes.2), name = "PeSC.score2.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.score2.1",cols = viridis(10), order = T) +
  ggtitle("PeSC Score from Cluster.6 scRNAseq")
f2 <- plot_density(sc.all, features = "PeSC.score2.1", joint = F, pal = "inferno") + 
  ggtitle("PeSC Score from Cluster.6 scRNAseq")
f1+f2



## PeSC score 146 genelist  ---------

# a shorter (146 genes) genelist was described here: https://pubmed.ncbi.nlm.nih.gov/26589433/
# after crossing the original data with ovarian CAFs (Table S1 below)

T.S1 <- read.csv("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/GSE17014_Table.S1.csv", row.names = 1)
head(T.S1)

table(T.S1$fc > 0) # 635
table(T.S1$p.value < 0.05) # 196

T.S1 <- T.S1[T.S1$p.value < 0.05, ]
T.S1 <- T.S1[T.S1$fc > 0, ]

T.S1 <- unique(toupper(T.S1$NAME)) # 104

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.score.TS1.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.score.TS1.1", cols = viridis(10), order = T) + 
  ggtitle("PeSC Score from Table S1 ovarian cancer paper")
f2 <- plot_density(sc.all, features = "PeSC.score.TS1.1", joint = F, pal = "inferno") + 
  ggtitle("PeSC Score from Table S1 ovarian cancer paper")
f1+f2


## PeSC score dermal pericytes  ---------

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17014

T.S1 <- read.csv("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/DEGs_bri.vs.dim_GSE17014.csv", row.names = 1)
head(T.S1)

table(T.S1$logFC > 0) # 888
table(T.S1$adj.P.Val < 0.05) # 2197
T.S1 <- T.S1[T.S1$logFC > 0, ]
T.S1 <- unique(toupper(T.S1$Symbol)) # 606

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.dermal.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.dermal.1",cols = viridis(10), order = T) +
  ggtitle("Dermal Pericyte Score (GSE17014)")
f2 <- plot_density(sc.all, features = "PeSC.dermal.1", joint = F, pal = "inferno") + 
  ggtitle("Dermal Pericyte Score (GSE17014)")
f1+f2



## Brain Pericytes score  ---------

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763650/

T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/mouse.brain.PEGs.csv", sep = ",")
head(T.S1) # 200

table(T.S1$Highest.expressionb)
T.S1 <- T.S1[T.S1$Highest.expressionb == "Pericytes", ]
T.S1 <- unique(toupper(T.S1$Gene)) # 62  or 200

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.brain.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.brain.1",cols = viridis(10), order = T) +
  ggtitle("Brain Pericyte Score")
f2 <- plot_density(sc.all, features = "PeSC.brain.1", joint = F, pal = "inferno") + 
  ggtitle("Brain Pericyte Score")
f1+f2



## Pericyte score Nkx3-2 paper  ---------

# https://www.sciencedirect.com/science/article/pii/S0012160618307942?via%3Dihub

# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7448
T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L2vsL3.genes.edgeR-tag.txt")
rownames(T.S1) <- make.names(T.S1$id, unique = T)
T.S1 <- na.omit(T.S1)
head(T.S1) # 29864 -> 15060

T.S1 <- T.S1[T.S1$EdgeR.tag.logFC > 0, ] # 7930
#T.S1 <- T.S1[T.S1$EdgeR.tag.logFC < 0, ] # 7930
T.S1 <- T.S1[T.S1$edgeR.tag.adjusted.p.value < 0.05, ] # 1815
T.S1 <- unique(toupper(T.S1$id)) # 2197

genes.present <- intersect(T.S1, rownames(sc.all)) # 330

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.L2L3.Up", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.Nkx.L2L3.Up1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L2.vs.L3.Up")
f2 <- plot_density(sc.all, features = "PeSC.Nkx.1", joint = F, pal = "inferno") + 
  ggtitle("Pericyte Score Nkx3-2 paper L2.vs.L3.Up")
f1+f2

T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L1vsL3.genes.edgeR-tag.txt")
rownames(T.S1) <- make.names(T.S1$id, unique = T)
T.S1 <- na.omit(T.S1)
head(T.S1) # 29864 -> 15060

T.S1 <- T.S1[T.S1$EdgeR.tag.logFC > 0, ] # 7930
#T.S1 <- T.S1[T.S1$EdgeR.tag.logFC < 0, ] # 7930
T.S1 <- T.S1[T.S1$edgeR.tag.adjusted.p.value < 0.05, ] # 1815
T.S1 <- unique(toupper(T.S1$id)) # 2197

genes.present <- intersect(T.S1, rownames(sc.all)) # 330

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.", ctrl = 100) # L1 vs L3 Up
f1 <- FeaturePlot(sc.all, features = "PeSC.Nkx.1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L3.Up")
f2 <- plot_density(sc.all, features = "PeSC.Nkx.1", joint = F, pal = "inferno") + 
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L3.Up")
f1+f2

T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L1vsL2.genes.edgeR-tag.txt")
rownames(T.S1) <- make.names(T.S1$id, unique = T)
T.S1 <- na.omit(T.S1)
head(T.S1) # 29864 -> 15060

T.S1 <- T.S1[T.S1$EdgeR.tag.logFC > 0, ] # 7930
#T.S1 <- T.S1[T.S1$EdgeR.tag.logFC < 0, ] # 7930
T.S1 <- T.S1[T.S1$edgeR.tag.adjusted.p.value < 0.05, ] # 1815
T.S1 <- unique(toupper(T.S1$id)) # 2197

genes.present <- intersect(T.S1, rownames(sc.all)) # 330

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.L1L2.Up", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.Nkx.L1L2.Up1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L2.Up")
f2 <- plot_density(sc.all, features = "PeSC.Nkx.1", joint = F, pal = "inferno") + 
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L2.Up")
f1+f2


head(sc.all[])


head(sc.all[])
save(sc.all, file="sc.all.RData")





# PDAC single cell Atlas: Stellate Cells --------------------------------------------------

rm(list=ls())

sc.all <- readRDS("~/Dropbox/BioInfo/freelance/StromaCare/scPADC/stella.rds")

head(sc.all[]) # 38919 features across 136163 samples
length(unique(sc.all$Patient)) # 72
table(sc.all$Patient)
table(sc.all$Type)
table(sc.all$Project)
table(Idents(sc.all))

d1 <- DimPlot(sc.all, label = T, label.size = 6, repel = T, label.color = "darkblue") + NoLegend()

DimPlot(sc.all, group.by = "Project")
DimPlot(sc.all, group.by = "Patient", label = T, repel = T)+ NoLegend()
d2 <- DimPlot(sc.all, group.by = "Type")

d1+d2

DefaultAssay(sc.all) <- "SCT"

sel.genes <- c("PDGFRB","ACE2","ABCC9","ANPEP","CD248",
               "CSPG4","DLK1","ENG","KCNJ8","RGS5",
               "DES","NES","ACTA2","NKX3-2", 
               "TGFBI","VEGFA","NGF","ACVR1B",
               "PDGFRA","FN1","COL1A1","DCN","PDPN",
               "CD34", "CD24","CD44",
               "EPCAM","KRT8","KRT18","KRT19",
               "CALD1","MCAM","ANGPT1","ANGPT2")
setdiff(sel.genes, rownames(sc.all))
int.1 <- intersect(sel.genes, rownames(sc.all))
jpeg("pericyte.markers_stella_feature.plot.jpg", width = 1500, height = 2000, quality = 100)
FeaturePlot(sc.all, sel.genes, order = T,  cols = viridis(10), ncol = 5)
dev.off()
jpeg("pericyte.markers_stella_density.plot.jpg", width = 1500, height = 1500, quality = 100)
plot_density(sc.all, int.1, joint = F, pal = "inferno")
dev.off()

for(i in 1:length(sel.genes)){
  
  jpeg(paste0("stella_feature.plot_", sel.genes[i], ".jpg"), 
       width = 600, height = 500, quality = 100)
  p1 <- FeaturePlot(sc.all, sel.genes[i], cols = viridis(10), order = T, pt.size = 0.5)
  plot(p1)
  dev.off()
  
}

for(i in 1:length(int.1)){
  
  jpeg(paste0("stella_density.plot_", int.1[i], ".jpg"), 
       width = 600, height = 500, quality = 100)
  p1 <- plot_density(sc.all, int.1[i], joint = F, pal = "inferno")
  plot(p1)
  dev.off()
  
}



## PeSC score Anca-curated  ---------


# list of PeSC markers curated by Anca using literature and enriched with Loupe analysis
genes <- read.xlsx("~/Dropbox/BioInfo/Colabs/PeSCs/Upregulated genes.xlsx", sheetIndex = 1)
head(genes)
genes <- unique(toupper(genes$gene_name))

#DefaultAssay(sc.all) <- "integrated"
DefaultAssay(sc.all) <- "SCT"

sc.all <- AddModuleScore(sc.all, features = list(genes), name = "PeSC.score", ctrl = 100)
# The following features are not present in the object: ACVR1B, NANOG, NA, SAA3, NKX3-2, ABCC9, CSPG4, DLK1, NGF, FOXD1, CD274, COL4A6, COL6A4, TGFB1, PDGFB, CTGF
f1 <- FeaturePlot(sc.all, features = "PeSC.score1",cols = viridis(10), order = T) + ggtitle("PeSC Score") + ggtitle("PeSC Score from curated genelist")
f2 <- plot_density(sc.all, features = "PeSC.score1", joint = F, pal = "inferno") + ggtitle("PeSC Score from curated genelist")
f1+f2

## PeSC score Loupe  ---------

# 2nd list (based on cluster 6 from loupe -> Anca)

genes.2 <- read.csv("~/Dropbox/BioInfo/Colabs/PeSCs/Cluster 6.csv", row.names = 1)
head(genes.2)
genes.2 <- unique(toupper(genes.2$FeatureName))

sc.all <- AddModuleScore(sc.all, features = list(genes.2), name = "PeSC.score2.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.score2.1",cols = viridis(10), order = T) +
  ggtitle("PeSC Score from Cluster.6 scRNAseq")
f2 <- plot_density(sc.all, features = "PeSC.score2.1", joint = F, pal = "inferno") + 
  ggtitle("PeSC Score from Cluster.6 scRNAseq")
f1+f2



## PeSC score 146 genelist  ---------

# a shorter (146 genes) genelist was described here: https://pubmed.ncbi.nlm.nih.gov/26589433/
# after crossing the original data with ovarian CAFs (Table S1 below)

T.S1 <- read.csv("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/GSE17014_Table.S1.csv", row.names = 1)
head(T.S1)

table(T.S1$fc > 0) # 635
table(T.S1$p.value < 0.05) # 196

T.S1 <- T.S1[T.S1$p.value < 0.05, ]
T.S1 <- T.S1[T.S1$fc > 0, ]

T.S1 <- unique(toupper(T.S1$NAME)) # 104

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.score.TS1.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.score.TS1.1", cols = viridis(10), order = T) + 
  ggtitle("PeSC Score from Table S1 ovarian cancer paper")
f2 <- plot_density(sc.all, features = "PeSC.score.TS1.1", joint = F, pal = "inferno") + 
  ggtitle("PeSC Score from Table S1 ovarian cancer paper")
f1+f2


## PeSC score dermal pericytes  ---------

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17014

T.S1 <- read.csv("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/DEGs_bri.vs.dim_GSE17014.csv", row.names = 1)
head(T.S1)

table(T.S1$logFC > 0) # 888
table(T.S1$adj.P.Val < 0.05) # 2197
T.S1 <- T.S1[T.S1$logFC > 0, ]
T.S1 <- unique(toupper(T.S1$Symbol)) # 606

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.dermal.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.dermal.1",cols = viridis(10), order = T) +
  ggtitle("Dermal Pericyte Score (GSE17014)")
f2 <- plot_density(sc.all, features = "PeSC.dermal.1", joint = F, pal = "inferno") + 
  ggtitle("Dermal Pericyte Score (GSE17014)")
f1+f2



## Brain Pericytes score  ---------

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763650/

T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/mouse.brain.PEGs.csv", sep = ",")
head(T.S1) # 200

table(T.S1$Highest.expressionb)
T.S1 <- T.S1[T.S1$Highest.expressionb == "Pericytes", ]
T.S1 <- unique(toupper(T.S1$Gene)) # 62  or 200

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.brain.", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.brain.1",cols = viridis(10), order = T) +
  ggtitle("Brain Pericyte Score")
f2 <- plot_density(sc.all, features = "PeSC.brain.1", joint = F, pal = "inferno") + 
  ggtitle("Brain Pericyte Score")
f1+f2



## Pericyte score Nkx3-2 paper  ---------

# https://www.sciencedirect.com/science/article/pii/S0012160618307942?via%3Dihub

# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7448
#T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L2vsL3.genes.edgeR-tag.txt")
#T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L1vsL3.genes.edgeR-tag.txt")
T.S1 <- read.delim("~/Dropbox/BioInfo/freelance/StromaCare/PeSCs/E-MTAB-7448/L1vsL2.genes.edgeR-tag.txt")
rownames(T.S1) <- make.names(T.S1$id, unique = T)
T.S1 <- na.omit(T.S1)
head(T.S1) # 29864 -> 15060

T.S1 <- T.S1[T.S1$EdgeR.tag.logFC > 0, ] # 7930
#T.S1 <- T.S1[T.S1$EdgeR.tag.logFC < 0, ] # 7930
T.S1 <- T.S1[T.S1$edgeR.tag.adjusted.p.value < 0.05, ] # 1815
T.S1 <- unique(toupper(T.S1$id)) # 2197

genes.present <- intersect(T.S1, rownames(sc.all)) # 330

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.L2L3.Up", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.Nkx.L2L3.Up1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L2.vs.L3.Up")
f2 <- plot_density(sc.all, features = "PeSC.Nkx.1", joint = F, pal = "inferno") + 
  ggtitle("Pericyte Score Nkx3-2 paper L2.vs.L3.Up")
f1+f2

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.", ctrl = 100) # L1 vs L3 Up
f1 <- FeaturePlot(sc.all, features = "PeSC.Nkx.1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L3.Up")
f2 <- plot_density(sc.all, features = "PeSC.Nkx.1", joint = F, pal = "inferno") + 
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L3.Up")
f1+f2

sc.all <- AddModuleScore(sc.all, features = list(T.S1), name = "PeSC.Nkx.L1L2.Up", ctrl = 100)
f1 <- FeaturePlot(sc.all, features = "PeSC.Nkx.L1L2.Up1",cols = viridis(10), order = T) +
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L2.Up")
f2 <- plot_density(sc.all, features = "PeSC.Nkx.1", joint = F, pal = "inferno") + 
  ggtitle("Pericyte Score Nkx3-2 paper L1.vs.L2.Up")
f1+f2


head(sc.all[])

saveRDS(sc.all, file = "stella.rds")




# GEO counts --------------------------------------------------------------

rm(list=ls())

load(file="~/Dropbox/BioInfo/Colabs/PeSCs/D.RData")

DimPlot(D)
head(D[[]])

#counts.normal <- GetAssayData(sc4, assay = "SCT", slot = "counts")
#counts.normal <- GetAssayData(sc4, assay = "RNA", slot = "counts")
counts.D <- GetAssayData(D, slot = "data")
head(counts.D)
counts.D[1:5,1:5]
counts.D <- as.data.frame(counts.D)
test <- rbind(as.character(D$seurat_clusters), counts.D)
row.names(test)[1] <- "Identity"
test[1:5,1:10]
write.csv(test, "counts_ductal.csv")



# end ---------------------------------------------------------------------
sessionInfo()




