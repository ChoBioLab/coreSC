# !/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(future) # parallelization

args <- commandArgs(trailingOnly = TRUE)
out_path <- paste0(args[1], "/")
objects <- list()
load(paste0(out_path, "tmp/preamble_image.RData"))

plan(
  multicore,
  workers = params["future.workers", ]
) # parallelization
options(future.globals.maxSize = params["future.mem", ] * 1024^2)

for (i in seq_len(nrow(samples))) {
  if (dir.exists(samples$dir[i])) {
    name <- samples$name[i]
    x <- Read10X( # pulling data with no filter
      data.dir = samples$dir[i]
    )
    str_section_head("Raw Object") # logging

    x <- CreateSeuratObject( # certain data will gen null matrix sans filters
      counts = x,
      project = samples$project[i],
      min.cells = params["min.cells", ],
      min.features = params["min.features", ]
    )

    x@meta.data$object <- samples$name[i]
    x@meta.data$group <- samples$group[i]
  } else {
    # TODO add compatibility with hdf5 format

    # it's assumed preformed objs have group and name metadata vars
    x <- readRDS(samples$dir[i])
    x
  }

  objects <- c(objects, x)
}

for (i in seq_len(length(objects))) {
  x <- objects[[i]]
  name <- unique(x$object)

  x[["percent.mt"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^MT-"
  )

  x[["percent.rb"]] <- PercentageFeatureSet(
    x,
    pattern = "RPS|RPL"
  )

  x[["percent.hb"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^HB[^(P)]"
  )

  p1 <- VlnPlot(
    x,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent.rb",
      "percent.hb"
    ),
    ncol = 5
  )

  save_figure(
    p1,
    paste0(name, "_unfilt_vln"),
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

  save_figure(
    (p1 + p2),
    paste0(name, "_unfilt_scatter"),
    width = 12,
    height = 6
  )

  str_section_head("Base Seurat Object") # logging

  x <- subset(
    x,
    nFeature_RNA > params["min.features", ] &
      nFeature_RNA < params["max.features", ] &
      nCount_RNA > params["min.count.rna", ] &
      nCount_RNA < params["max.count.rna", ] &
      percent.mt < params["max.percent.mt", ] &
      percent.mt > params["min.percent.mt", ]
  )

  # norm, dimred, and clustering

  # clustering is performed on individual samples for QC
  x <- NormalizeData(object = x) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(
      npcs = d,
      verbose = FALSE
    ) %>%
    RunUMAP(
      reduction = "pca",
      dims = 1:d,
      verbose = FALSE
    ) %>%
    FindNeighbors(
      reduction = "pca",
      dims = 1:d,
      verbose = FALSE
    ) %>%
    FindClusters(
      resolution = params["res", ],
      verbose = FALSE
    )

  str_section_head("Subset, SCT Normd, Red, Clust") # logging

  top10 <- head(VariableFeatures(x), 10)
  p1 <- VariableFeaturePlot(x)
  p2 <- LabelPoints(
    plot = p1,
    points = top10,
    repel = TRUE
  )

  save_figure(
    (p1 + p2),
    paste0(name, "_var_features"),
    width = 12,
    height = 6
  )

  p1 <- DimPlot(
    x,
    reduction = "umap",
    group.by = "group"
  )

  p2 <- DimPlot(
    x,
    reduction = "umap",
    label = TRUE
  )

  save_figure(
    (p1 + p2),
    paste0(name, "_individual_dimplot"),
    width = 12,
    height = 6
  )

  objects[[i]] <- x
}

save_object(objects, "individual_clustered")

print("End of create-object-norm.R")
