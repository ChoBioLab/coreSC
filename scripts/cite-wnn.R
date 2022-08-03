# !/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(future) # parallelization

plan(multicore) # parallelization
options(future.globals.maxSize = 2000 * 1024^2)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")

load(paste0(out_path, "tmp/preamble_image.RData"))

for (i in 1:nrow(samples)) {
  raw <- Read10X( # pulling data with no filter
    data.dir = samples$dir[i]
  )
  str_section_head("Raw Object") # logging

  rna <- raw$`Gene Expression`
  adt <- raw$`Antibody Capture`

  x <- CreateSeuratObject(
    counts = rna,
    project = samples$project[i]
  )

  x[["ADT"]] <- CreateAssayObject(counts = adt)

  x[["percent.mt"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^MT-"
  )

  # raw qc plots
  p1 <- VlnPlot(
    x,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3
  )

  save_figure(
    p1,
    paste0(samples$name[i], "_rna_unfilt_vln"),
    width = 12,
    height = 6
  )

  p1 <- VlnPlot(
    x,
    features = c(
      "nFeature_ADT",
      "nCount_ADT",
    ),
    ncol = 2
  )

  save_figure(
    p1,
    paste0(samples$name[i], "_adt_unfilt_vln"),
    width = 12,
    height = 6
  )

  p1 <- FeatureScatter(
    x,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )

  p2 <- FeatureScatter(
    x,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

  p3 <- FeatureScatter(
    x,
    feature1 = "nCount_ADT",
    feature2 = "nFeature_ADT"
  )

  save_figure(
    (p1 + p2),
    paste0(samples$name[i], "_unfilt_scatter"),
    width = 12,
    height = 6
  )

  str_section_head("Base Seurat Object") # logging

  x <- subset(
    x,
    nFeature_RNA > params["min.features", ] &
      nFeature_RNA < params["max.features", ] &
      nCount_RNA > params["min.count", ] &
      nCount_RNA < params["max.count", ] &
      percent.mt < params["max.percent.mt", ] &
      percent.mt > params["min.percent.mt", ]
  )

  # norm, dimred, and clustering
  # clustering is performed on individual samples for RNA QC
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(
    x,
    vst.flavor = "v2",
    verbose = FALSE
  ) %>%
    RunPCA(
      npcs = d,
      verbose = FALSE
    )

  DefaultAssay(x) <- "ADT"
  VariableFeatures(x) <- rownames(x[["ADT"]])
  x <- NormalizeData(
    x,
    normalization.method = "CLR",
    margin = 2
  ) %>%
    ScaleData() %>%
    RunPCA(reduction.name = "apca")

  # neighborhood formation and clustering
  x <- FindMultiModalNeighbors(
    x,
    reduction.list = list("pca", "apca"),
    dims.list = list(1:d, 1:18),
    modality.weight.name = "RNA.weight"
  )

  x <- RunUMAP(
    x,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_"
  )

  x <- FindClusters(
    x,
    graph.name = "wsnn",
    algorithm = 3,
    resolution = params["res", ],
    verbose = FALSE
  )

  str_section_head("Subset, SCT Normd, Red, Clust") # logging

  p1 <- DimPlot(
    x,
    reduction = "wnn.umap",
    label = TRUE,
    repel = TRUE,
    label.size = 2.5
  ) + NoLegend()

  # clustering by modality
  x <- RunUMAP(
    x,
    reduction = "pca",
    dims = 1:d,
    assay = "RNA",
    reduction.name = "rna.umap",
    reduction.key = "rnaUMAP_"
  )

  x <- RunUMAP(
    x,
    reduction = "apca",
    dims = 1:18,
    assay = "ADT",
    reduction.name = "adt.umap",
    reduction.key = "adtUMAP_"
  )

  p3 <- DimPlot(
    x,
    reduction = "rna.umap",
    label = TRUE,
    repel = TRUE,
    label.size = 2.5
  ) + NoLegend()

  p4 <- DimPlot(
    x,
    reduction = "adt.umap",
    label = TRUE,
    repel = TRUE,
    label.size = 2.5
  ) + NoLegend()

  top10 <- head(VariableFeatures(x), 10)
  p1 <- VariableFeaturePlot(x)
  p2 <- LabelPoints(
    plot = p1,
    points = top10,
    repel = TRUE
  )

  save_figure(
    (p1 + p2),
    paste0(samples$name[i], "_var_features"),
    width = 12,
    height = 6
  )

  x@meta.data$object <- samples$name[i]
  x@meta.data$group <- samples$group[i]

  # plotting clusters by specified group
  p1 <- DimPlot(
    x,
    reduction = "wnn.umap",
    group.by = "group",
    label = T
  )

  p2 <- DimPlot(
    x,
    reduction = "wnn.umap",
    label = T
  )

  save_figure(
    (p1 + p2),
    paste0(samples$name[i], "individual_dimplot"),
    width = 12,
    height = 6
  )

  assign( # giving names to objects
    samples$name[i],
    x
  )
}

if (length(samples$name) == 1) {
  save_object(x, "individual_clustered")
} else { # integrate
  # create and save list of seurat objects
  objects <- list()
  for (i in samples$name) {
    objects <- c(
      objects,
      get(i) # need get() to call object instead of string
    )
  }

  save_object(objects, "individual_clustered")
}

sessionInfo()

print("End of cite-wnn.R")