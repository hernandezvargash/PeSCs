
# Description -------------------------------------------------------------

# scRNAseq on mouse CAFs and Ducts
# two progressing mouse models of PDAC (KC and 4KC) 
# preliminary analyses did not show high PeSC score in KC CAFs
# so, only using KC Duct samples to identify a putative pericyte stem cell (PeSC) subpopulation

# associated manuscript: https://pubmed.ncbi.nlm.nih.gov/36802267/
# associated dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220687


# Libraries ---------------------------------------------------------------

rm(list=ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(sctransform)
  library(SingleR)
  library(xlsx)
  library(viridis)
  library(gridExtra)
  library(RColorBrewer)
  library(enrichR)
  library(pathfindR)
  library(hrbrthemes)
  library(reshape)
  library(stringr)
  library(celldex)
  library(scRNAseq)
  library(scuttle)
  library(Nebulosa)
  library(oligo)
  library(biomaRt)
  library(limma)
  
})

setwd()
list.files()

set.seed(1234)



# Pre-processing (Fig2B) ----------------------------------------------------------

raw <- Read10X(data.dir = "filtered_feature_bc_matrix/") # from GSE220687
rawD <- CreateSeuratObject(raw) # 1169 cells
rawD$condition <- "Ducts"

D <- PercentageFeatureSet(rawD, pattern = "^mt-", col.name = "percent.mt") %>%
  PercentageFeatureSet(pattern = "^Rp[sl][[:digit:]]", col.name = "percent.ribo") %>%
  subset(subset = percent.mt < 10) %>% # PeSC cluster seems to have a lower RNA content, so no lower filter is applied at this point
#  subset(subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10) %>%
  SCTransform(vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:20) %>% RunTSNE()

# 886 cells left after a percent.mt filter of < 10%
VlnPlot(D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 2)

plot1 <- FeatureScatter(D, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(D, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(D, feature1 = "percent.mt", feature2 = "percent.ribo", pt.size = 1.5)
plot4 <- FeatureScatter(D, feature1 = "percent.ribo", feature2 = "nFeature_RNA", pt.size = 1.5)
plot_grid(plot1, plot2, plot3, plot4)

jpeg("Fig2B.jpeg", height = 960, width = 960, quality = 100)
DimPlot(D, pt.size = NULL, label = T, label.box = T)
dev.off()

table(D$seurat_clusters)
#  0   1   2   3   4   5   6   7   8   9 
# 161 149 137 106  83  77  62  48  39  24 



# Markers (Fig2C) -----------------------------------------------------------------

all.markers <- FindAllMarkers(D, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
all.markers <- all.markers[all.markers$p_val_adj<0.05, ]
head(all.markers)
tail(all.markers)

write.xlsx(all.markers, file="all.markers.ducts.xlsx")
write.csv(all.markers, file="all.markers.ducts.csv")

top.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
jpeg("FeaturePlot.top.markers.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = top.markers$gene)
dev.off()
jpeg("VlnPlot.top.markers.jpeg", height = 960, width = 960, quality = 100)
VlnPlot(D, features = top.markers$gene)#, cols = show_col(hue_pal()(16))) # DiscretePalette(12, "glasbey"))
dev.off()
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
jpeg("Fig2C.jpeg", height = 960, width = 960, quality = 100)
DoHeatmap(D, features = top10$gene)#, group.colors =  rainbow(11)) #+ NoLegend()
dev.off()




# Pathways (Fig2D) ----------------------------------------------------------------

rm(list=ls())

markers <- read.csv("all.markers.ducts.csv",row.names = 1)
head(markers)

# pathways cluster 7 (putative PeSCs)
genes <- markers[markers$cluster==7, ]
genes <- genes[genes$avg_log2FC > 0, ]
head(genes)
genes <- unique(genes$gene)


# enrichr

listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "BioCarta_2016")

enriched <- enrichr(genes, dbs)

head(enriched[["KEGG_2019_Mouse"]])
head(enriched[["WikiPathways_2019_Mouse"]])
head(enriched[["BioCarta_2016"]])

K <- enriched[["KEGG_2019_Mouse"]]
write.xlsx(K, "KEGG_2019_Mouse_cluster.7.xlsx")
W <- enriched[["WikiPathways_2019_Mouse"]] # Fig 2D
write.xlsx(W, "WikiPathways_2019_Mouse_cluster.7.xlsx")
B <- enriched[["BioCarta_2016"]]
write.xlsx(B, "BioCarta_2016_cluster.7.xlsx")


# save all cluster 7 markers
head(markers)
markers.clus <- markers[markers$cluster==7, ]
write.csv(markers.clus, "markers.cluster.7.csv")




# Scores ----

## New PeSC Score (Fig3H) --------------------------------------------------------------

# list of curated PeSC markers from the literature and bulk RNAseq analysis
genes <- read.xlsx("Upregulated genes.xlsx", sheetIndex = 1)
head(genes)
genes <- unique(genes$gene_name)

jpeg("FeaturePlots.1.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = genes[1:12])
dev.off()
jpeg("FeaturePlots.2.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = genes[13:24])
dev.off()
jpeg("FeaturePlots.3.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = genes[25:36])
dev.off()
jpeg("FeaturePlots.4.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = genes[37:48])
dev.off()
jpeg("FeaturePlots.5.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = genes[49:60])
dev.off()

D <- AddModuleScore(D, features = list(genes), name = "PeSC.score", ctrl = 100)
names(x = D[[]])
jpeg("Fig3H.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = "PeSC.score1",cols = viridis(4), order = T, pt.size = 1) + ggtitle("PeSC Score")
dev.off()


# PeSC markers (cluster 7)

cluster.7 <- all.markers[all.markers$cluster == "7", "gene"]
intersect(genes, cluster.7)
setdiff(genes, cluster.7)
intersect(genes.2$FeatureName, cluster.7)
setdiff(genes.2$FeatureName, cluster.7)


FeaturePlot(D, feature = cluster.7[1:9])



## Ovarian PeSC score (FigEV2 A)  ---------

# a shorter (146 genes) genelist was described here: https://pubmed.ncbi.nlm.nih.gov/26589433/
# after crossing the original data with ovarian CAFs (Table S1 below)

T.S1 <- read.csv("GSE17014_Table.S1.csv", row.names = 1)
head(T.S1)

table(T.S1$fc > 0) # 635
table(T.S1$p.value < 0.05) # 196

T.S1 <- T.S1[T.S1$p.value < 0.05, ]
T.S1 <- T.S1[T.S1$fc > 0, ]

T.S1 <- unique(toupper(T.S1$NAME)) # 104

D <- AddModuleScore(D, features = list(T.S1), name = "PeSC.ovarian.", ctrl = 100)
f1 <- FeaturePlot(D, features = "PeSC.ovarian.1", cols = viridis(10), order = T) + 
  ggtitle("ovarian cancer pericyte score")



## Brain Pericytes score (FigEV2 B)  ---------

load(file="D.RData")

# mouse brain specific pericyte expressed genes were downloaded from:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8763650/ (Table S1)

T.S1 <- read.delim("mouse.brain.PEGs.csv", sep = ",")
head(T.S1) # 200

table(T.S1$Highest.expressionb)
T.S1 <- T.S1[T.S1$Highest.expressionb == "Pericytes", ]
T.S1 <- unique(toupper(T.S1$Gene)) # 62  or 200

D <- AddModuleScore(D, features = list(T.S1), name = "PeSC.brain.", ctrl = 100)
f2 <- FeaturePlot(D, features = "PeSC.brain.1",cols = viridis(10), order = T) +
  ggtitle("Brain cancer pericyte score")

grid.arrange(f1, f2, ncol=2)



save(D, file="D.RData")




# GSA-PDAC ----------------------------------------------------------------

# PRJCA001063 / CRA001160
# https://ngdc.cncb.ac.cn/gsa/browse/CRA001160

# single cell PDAC data corresponding to manuscripts:
# https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00108-0/fulltext
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6796938/

## preprocessing ----

samples <- list.files(path = "PRJCA001063/PDAC")

list.samples <- list()
for(i in 1:length(samples)){
  raw1 <- Read10X(data.dir = paste0("PRJCA001063/PDAC/", samples[i]))
  raw2 <- CreateSeuratObject(raw1)
  raw2$sample <- samples[i]
  list.samples[[i]] <- raw2
}
names(list.samples) <- samples
lapply(list.samples, dim)


check.one.sample <- PercentageFeatureSet(list.samples[[30]], pattern = "^MT", col.name = "percent.mt")
check.one.sample <- PercentageFeatureSet(check.one.sample, pattern = "^RP[SL]", col.name = "percent.ribo")
VlnPlot(check.one.sample, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), ncol = 4)


# Integrate
# following: https://satijalab.org/seurat/archive/v3.0/integration.html

bm40k.list <- lapply(X = list.samples, FUN = function(x) {
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
  x <- PercentageFeatureSet(x, pattern = "^RP[SL]", col.name = "percent.ribo")
  x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25) # increase percent.mt to 25%, see Neou et al Cancer Cell 2020
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = bm40k.list)
bm40k.list <- lapply(X = bm40k.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = bm40k.list, reference = c(1), reduction = "rpca", # selected sample 1 as reference for no particular reason
                                  dims = 1:50) # "increased the dimensionality to 50 to reflect the increased cell number and diversity."
sc.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
# 34738 features across 84433 samples

sc.integrated <- ScaleData(sc.integrated, verbose = FALSE)
sc.integrated <- RunPCA(sc.integrated, verbose = FALSE)
sc.integrated <- RunUMAP(sc.integrated, dims = 1:50)
sc.integrated <- FindNeighbors(sc.integrated, reduction = "pca", dims = 1:30)
sc.integrated <- FindClusters(sc.integrated, resolution = 0.5)



## PeSC score (Fig3I) ----

p1 <- DimPlot(sc.integrated, reduction = "umap", group.by = "sample")
p2 <- DimPlot(sc.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
p1 + p2

table(sc.integrated$sample)
sc.integrated$condition <- str_sub(sc.integrated$sample, 1, 1)
table(sc.integrated$condition)
#   N     T 
# 32002 52431 
head(sc.integrated[[]])

DimPlot(sc.integrated, reduction = "umap", split.by = "condition")

pesc.genes <- read.csv("markers.cluster.7.csv")
pesc.genes <- toupper(pesc.genes$gene)

common.genes <- intersect(pesc.genes, rownames(sc.integrated))
# write.xlsx(common.genes, file = "PeSC.genes.in.PDAC.xlsx")

DefaultAssay(sc.integrated) <- "RNA"
sc.integrated <- AddModuleScore(sc.integrated, features = list(common.genes), name = "PeSC_Score")
head(sc.integrated[[]])

DefaultAssay(sc.integrated) <- "integrated"

FeaturePlot(sc.integrated, "PeSC_Score1", order = T)
FeaturePlot(sc.integrated, "PeSC_Score1", split.by = "condition", order = T) # Fig 3I



## cell type identification (FigEV2 C-F) ----

ref <- MuraroPancreasData() # human pancreas dataset: https://pubmed.ncbi.nlm.nih.gov/27693023/
ref <- ref[,!is.na(ref$label)]
ref <- logNormCounts(ref)
rownames(ref) <- str_split_fixed(rownames(ref), "__", 2)[, 1]
table(ref$label)

pred.R <- SingleR(test = as.SingleCellExperiment(sc.integrated), ref = ref, labels = ref$label, de.method="wilcox")
pred.R
table(pred.R$labels)

par(mfrow = c(1,1), mar=c(10,4,4,4))
barplot(sort(table(pred.R$labels), decreasing = T), col = viridis(10), las = 2, main = "SingleR - MuraroPancreasData", border = 0)

sc.integrated <- AddMetaData(sc.integrated, pred.R$labels, "MuraroPancreasData")
head(sc.integrated[[]])

DimPlot(sc.integrated, group.by = "MuraroPancreasData", label = T, repel = T)# + NoLegend()
DimPlot(sc.integrated, group.by = "MuraroPancreasData", label = T, repel = T, split.by = "condition") + NoLegend()


save(sc.integrated, file = "sc.integrated_PRJCA001063.RData")
#load("sc.integrated_PRJCA001063.RData")



# GEO counts --------------------------------------------------------------

rm(list=ls())

load(file="D.RData")

DimPlot(D)
head(D[[]])

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

