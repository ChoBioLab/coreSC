library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)

#load datasets
S1.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S1raw_feature_bc_matrix/")
S2.data <- Read10X(data.dir = "/Users/shikhanayar/Dropbox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/Matrices/S2raw_feature_bc_matrix/")

#initialize each Seurat object
S1 <- CreateSeuratObject(counts = S1.data, project = "WTctrl1X", min.cells = 3, min.features = 200)
S2 <- CreateSeuratObject(counts = S2.data, project = "nod2ctrl1X", min.cells = 3, min.features = 200)

#mitochondrial cutoffs
S1[["percent.mt"]] <- PercentageFeatureSet(S1, pattern = "^mt-")
S2[["percent.mt"]] <- PercentageFeatureSet(S2, pattern = "^mt-")

#Visualize QC metrics as violin plot
VlnPlot(S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Correlation between feature, count, mito
plot1 <- FeatureScatter(S1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(S1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#set up each object after choosing QC parameters
S1$group <- "WT"
S1 <- subset(S1, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S1 <- NormalizeData(S1, normalization.method = "LogNormalize")
S1 <- FindVariableFeatures(S1, selection.method = "vst", nfeatures = 2000)

S2$group <- "nod2"
S2 <- subset(S2, subset = nFeature_RNA > 150 & nFeature_RNA < 4000 & percent.mt < 60)
S2 <- NormalizeData(S2, normalization.method = "LogNormalize")
S2 <- FindVariableFeatures(S2, selection.method = "vst", nfeatures = 2000)

#perform integration of all samples
immune.anchors <- FindIntegrationAnchors(object.list = list(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10), dims = 1:20)
new.fish <- IntegrateData(anchorset = immune.anchors)

#perform integrated analysis
DefaultAssay(new.fish) <- "integrated"
new.fish <- ScaleData(new.fish, verbose = FALSE)
new.fish <- RunPCA(new.fish, verbose = FALSE)
ElbowPlot(new.fish)
JackStrawPlot(new.fish, dims = 1:15)

#tSNE and clustering
new.fish <- RunUMAP(new.fish, reduction = "pca", dims = 1:15)
new.fish <- FindNeighbors(new.fish, reduction = "pca", dims = 1:15)
new.fish <- FindClusters(new.fish, resolution = 0.8)

#visualization
p1 <- DimPlot(new.fish, reduction = "umap", group.by = "group")
p2 <- DimPlot(new.fish, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(new.fish, reduction = "umap", split.by = "group")

#savefile
saveRDS(new.fish, file = "/Users/shikhanayar/Ddrobox/Mamta_Shikha/Zebrafish_scRNAseq_summer2020/new.fish.rds")


