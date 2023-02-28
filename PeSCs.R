
# Description -------------------------------------------------------------

# 121120 / 160322
# Mouse CAFs and Ducts
# two progressing mouse models of PDAC (KC and 4KC) 
# only using KC Duct samples to identify a putative pericyte stem cell (PeSC) subpopulation
# preliminary analyses did not show high PeSC score in KC CAFs
# run on @nanotower

# associated manuscript:
# CD106+ pericyte stem cells lead to Ly6G+ suppressor cell accumulation and sustained resistance to PD-1 immunotherapy in pancreatic cancer

# update for paper revision November 2022


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
  
})

setwd("~/Dropbox/BioInfo/Colabs/PeSCs")
list.files()

set.seed(1234)


# Pre-processing ----------------------------------------------------------

raw <- Read10X(data.dir = "raw/D/filtered_feature_bc_matrix/")
rawD <- CreateSeuratObject(raw) # 1169 cells
rawD$condition <- "Ducts"

#load("D.RData")

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

jpeg("DimPlot.after.filtering.jpeg", height = 960, width = 960, quality = 100)
DimPlot(D, pt.size = NULL, label = T, label.box = T)
dev.off()

table(D$seurat_clusters)
#  0   1   2   3   4   5   6   7   8   9 
# 161 149 137 106  83  77  62  48  39  24 



# Markers -----------------------------------------------------------------

all.markers <- FindAllMarkers(D, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
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
jpeg("Heatmap.top.markers.jpeg", height = 960, width = 960, quality = 100)
DoHeatmap(D, features = top10$gene)#, group.colors =  rainbow(11)) #+ NoLegend()
dev.off()



# PeSC Score --------------------------------------------------------------

# list of PeSC markers curated by Anca using literature and enriched with Loupe analysis
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
jpeg("PeSC.score.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = "PeSC.score1",cols = viridis(4), order = T, pt.size = 1) + ggtitle("PeSC Score")
dev.off()
#FeaturePlot(sc.merged, features = "PeSC.score1",cols = rev(brewer.pal(n = 11, name = "Spectral")), sort.cell = T, pt.size = 1)

save(D, file="D.RData")


# 2nd list (based on cluster 6 from loupe -> Anca)

genes.2 <- read.csv("Cluster 6.csv", row.names = 1)
head(genes.2)
tail(genes.2)

D <- AddModuleScore(D, features = list(genes.2$FeatureName), name = "PeSC.score.2", ctrl = 100)
head(x = D[[]])
jpeg("PeSC.score.2.jpeg", height = 960, width = 960, quality = 100)
FeaturePlot(D, features = "PeSC.score.21",cols = viridis(4), order = T, pt.size = 1) + ggtitle("PeSC Score 2")
dev.off()

save(D, file="D.RData")



cluster.7 <- all.markers[all.markers$cluster == "7", "gene"]
intersect(genes, cluster.7)
setdiff(genes, cluster.7)
intersect(genes.2$FeatureName, cluster.7)
setdiff(genes.2$FeatureName, cluster.7)


FeaturePlot(D, feature = cluster.7[1:9])



# Pathways ----------------------------------------------------------------

rm(list=ls())

listEnrichrDbs()
dbs <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "BioCarta_2016")

markers <- read.csv("all.markers.ducts.csv",row.names = 1)
head(markers)

# pathways cluster 7 (putative PeSCs)

genes <- markers[markers$cluster==7, ]
genes <- genes[genes$avg_log2FC > 0, ]
head(genes)
genes <- unique(genes$gene)


# enrichr

enriched <- enrichr(genes, dbs)
head(enriched[["KEGG_2019_Mouse"]])
head(enriched[["WikiPathways_2019_Mouse"]])
head(enriched[["BioCarta_2016"]])

K <- enriched[["KEGG_2019_Mouse"]]
write.xlsx(K, "KEGG_2019_Mouse_cluster.7.xlsx")
W <- enriched[["WikiPathways_2019_Mouse"]]
write.xlsx(W, "WikiPathways_2019_Mouse_cluster.7.xlsx")
B <- enriched[["BioCarta_2016"]]
write.xlsx(B, "BioCarta_2016_cluster.7.xlsx")


# pathfindR

pathways <- c("mmu_KEGG","Reactome","BioCarta")

head(markers)
markers.clus <- markers[markers$cluster==7, ]
input_df <- markers.clus[,c(2,5)]
#input_df$Gene.symbol <- toupper(rownames(input_df))
input_df$Gene.symbol <- rownames(input_df)
input_df <- input_df[,c(3,1,2)]
colnames(input_df) <- c("Gene.symbol","logFC","adj.P.Val")
head(input_df)

for(p in pathways){
  output_df <- run_pathfindR(input_df, output_dir = paste0(p, "_ducts.markers_cluster_7"), gene_sets = p)
  if(dim(output_df)[1] > 0){
    jpeg(filename = paste0(p, "_ducts.markers_cluster_7", ".jpeg"), width = 600, height = 900, quality = 100)
    g <- enrichment_chart(result_df = output_df, top_terms = 10)
    plot(g)
    dev.off()  
  }
  rm(output_df)
}



write.csv(markers.clus, "markers.cluster.7.csv")




# GSE165399 --------------------------------------------

rm(list=ls())

setwd("~/Dropbox/BioInfo/Colabs/scCAFs")
list.files()

#Sample 1_IPMN_ZX-2
#tissue type: Intraductal papillary mucinous neoplasm
#age: 74 years
#gender: Male
#cell number: 6276

#Sample 2_PASC_ZX-3
#tissue type: Pancreatic adenosquamous carcinoma
#age: 59 years
#gender: Male
#cell number: 2215

#Sample 3_normal pancreas_ZX-4
#tissue type: normal pancreas sample
#age: 50 years
#gender: Male
#cell number: 1855

raw1 <- read.delim("GSE165399_RAW/GSM5032771_ZX-2.expression_matrix.txt.gz")
raw1[1:5,1:10]
dim(raw1)
raw2 <- read.delim("GSE165399_RAW/GSM5032772_ZX-3.expression_matrix.txt.gz")
raw2[1:5,1:10]
dim(raw2)
raw3 <- read.delim("GSE165399_RAW/GSM5032773_ZX-4.expression_matrix.txt.gz")
raw3[1:5,1:10]
dim(raw3)

setwd("~/Dropbox/BioInfo/Colabs/PeSCs")

ipmn <- CreateSeuratObject(raw1)
ipmn$sample <- "IPMN"
pasc <- CreateSeuratObject(raw2)
pasc$sample <- "PASC"
norm <- CreateSeuratObject(raw3)
norm$sample <- "Norm"

sc.merged <- merge(ipmn, y = list(pasc,norm), add.cell.ids = c("IPMN", "PASC", "Normal"), project = "GSE165399")
# 38034 features across 10346 samples

sc.merged <- PercentageFeatureSet(sc.merged, pattern = "^MT-", col.name = "percent.mt") %>%
  PercentageFeatureSet(pattern = "^RP[SL]", col.name = "percent.ribo")

VlnPlot(sc.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

sc.merged <-  subset(sc.merged, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15) %>%
  SCTransform(vars.to.regress = "percent.mt") %>% RunPCA() %>% FindNeighbors(dims = 1:5, k.param = 30) %>% 
  FindClusters(resolution = 0.5) %>% RunUMAP(dims = 1:5) %>% RunTSNE()
# 61352 features across 8186 samples within 2 assays

head(sc.merged[[]])

table(sc.merged$sample)
# IPMN Norm PASC 
# 5552  636 1998

p1 <- DimPlot(sc.merged)
p2 <- DimPlot(sc.merged, reduction = "umap", group.by = "sample")#, pt.size = 1)
p1+p2
DimPlot(sc.merged, split.by = "sample")

pesc.genes <- read.csv("~/Dropbox/BioInfo/Colabs/PeSCs/markers.cluster.7.csv")

sc.merged <- AddModuleScore(sc.merged, features = list(pesc.genes$gene), name = "PeSC_Score")
head(sc.merged[[]])

FeaturePlot(sc.merged, "PeSC_Score1", order = T)
FeaturePlot(sc.merged, "PeSC_Score1", split.by = "sample", order = T)


save(sc.merged, file = "sc.merged_GSE165399.RData")




# GSA-PDAC ----------------------------------------------------------------

# PRJCA001063 / CRA001160
# https://ngdc.cncb.ac.cn/gsa/browse/CRA001160

# single cell PDAC data corresponding to manuscripts:
# https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(21)00108-0/fulltext
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6796938/

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




# not done this time but consider for new analyses
# slower (but more accurate because of SCT normalization?) integration method:
# https://satijalab.org/seurat/archive/v3.0/integration.html

for(p in 1:length(list.samples)){
  list.samples[[p]] <- PercentageFeatureSet(list.samples[[p]], pattern = "^MT-", col.name = "percent.mt") %>%
    PercentageFeatureSet(pattern = "^RP[SL]", col.name = "percent.ribo") %>%
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 25) %>% # increase percent.mt to 25%, see Neou et al Cancer Cell 2020
    SCTransform(variable.features.n = 3000) %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(dims = 1:20, reduction = "pca", k.param = 30) %>%
    FindClusters(resolution = 0.5) %>% 
    RunUMAP(dims = 1:30)
}

DimPlot(list.samples[[10]])

features <- SelectIntegrationFeatures(object.list = list.samples, nfeatures = 3000)
list.samples <- PrepSCTIntegration(object.list = list.samples, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list.samples, normalization.method = "SCT", anchor.features = features,
                                  reference = 1) # added a reference as suggested here: https://github.com/satijalab/seurat/issues/1029
sc.integrated <- IntegrateData(anchorset = anchors, dims = 1:30, normalization.method = "SCT")
#DefaultAssay(sc.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering

sc.integrated <- RunPCA(sc.integrated, npcs = 30, verbose = FALSE)
sc.integrated <- RunUMAP(sc.integrated, reduction = "pca", dims = 1:30)
sc.integrated <- FindNeighbors(sc.integrated, reduction = "pca", dims = 1:30)
sc.integrated <- FindClusters(sc.integrated, resolution = 0.5)



# Inspect

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


pesc.genes <- read.csv("~/Dropbox/BioInfo/Colabs/PeSCs/markers.cluster.7.csv")
pesc.genes <- toupper(pesc.genes$gene)

common.genes <- intersect(pesc.genes, rownames(sc.integrated))
# write.xlsx(common.genes, file = "PeSC.genes.in.PDAC.xlsx")

DefaultAssay(sc.integrated) <- "RNA"
sc.integrated <- AddModuleScore(sc.integrated, features = list(common.genes), name = "PeSC_Score")
head(sc.integrated[[]])

DefaultAssay(sc.integrated) <- "integrated"

FeaturePlot(sc.integrated, "PeSC_Score1", order = T)
FeaturePlot(sc.integrated, "PeSC_Score1", split.by = "condition", order = T)


# SingleR

#ref <- HumanPrimaryCellAtlasData()

ref <- MuraroPancreasData() # human pancreas dataset: https://pubmed.ncbi.nlm.nih.gov/27693023/
# One should normally do cell-based quality control at this point, but for
# brevity's sake, we will just remove the unlabelled libraries here.
ref <- ref[,!is.na(ref$label)]
# SingleR() expects reference datasets to be normalized and log-transformed.
ref <- logNormCounts(ref)
rownames(ref) <- str_split_fixed(rownames(ref), "__", 2)[, 1]

table(ref$label)

# We then run SingleR() as described previously but with a marker detection mode that considers the variance of expression across cells. 
# Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. 
# This is slower but more appropriate for single-cell data compared to the default marker detection algorithm 
# (which may fail for low-coverage data where the median is frequently zero).
pred.R <- SingleR(test = as.SingleCellExperiment(sc.integrated), ref = ref, labels = ref$label, de.method="wilcox")
pred.R
table(pred.R$labels)

# pred.R <- SingleR(test = as.SingleCellExperiment(sc.integrated), ref = ref, labels = ref$label.main)
#sc.integrated <- AddMetaData(sc.integrated, pred.R$labels, "HumanPrimaryCellAtlasData")

par(mfrow = c(1,1), mar=c(10,4,4,4))
#barplot(sort(table(pred.R$labels), decreasing = T), col = viridis(10), las = 2, main = "SingleR - HumanPrimaryCellAtlasData", border = 0)
barplot(sort(table(pred.R$labels), decreasing = T), col = viridis(10), las = 2, main = "SingleR - MuraroPancreasData", border = 0)

sc.integrated <- AddMetaData(sc.integrated, pred.R$labels, "MuraroPancreasData")
head(sc.integrated[[]])

DimPlot(sc.integrated, group.by = "MuraroPancreasData", label = T, repel = T)# + NoLegend()
DimPlot(sc.integrated, group.by = "MuraroPancreasData", label = T, repel = T, split.by = "condition") + NoLegend()



save(sc.integrated, file = "sc.integrated_PRJCA001063.RData")



# downloaded cell type annotations for the same dataset

cell.types <- read.delim("PRJCA001063/all_celltype.txt")
head(cell.types)
tail(cell.types)
table(cell.types$cluster)

int1 <- intersect(colnames(sc.integrated), paste0(str_sub(cell.types$cell.name, 4, -1), "-1_1"))
diff1 <- setdiff(colnames(sc.integrated), paste0(str_sub(cell.types$cell.name, 4, -1), "-1_1"))





# end ---------------------------------------------------------------------
sessionInfo()

