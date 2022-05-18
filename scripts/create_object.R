# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(future) # parallelization

plan(multicore) # parallelization
options(future.globals.maxSize = 2000 * 1024^2)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")

load(paste0(out_path, "tmp/preamble_image.RData"))

for (i in 1:nrow(samples)) {
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

  x[["percent.mt"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^MT-"
  )

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
    paste0(samples$name[i], "_unfilt_vln"),
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

  x <- NormalizeData(x)
  str_section_head("Subset, Normalized") # logging

  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = params["max.features", ]
  )

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

  genes <- rownames(x)
  x <- ScaleData(
    x,
    verbose = FALSE,
    features = genes
  )

  x@meta.data$object <- samples$name[i]
  x@meta.data$group <- samples$group[i]

  assign( # giving names to objects
    samples$name[i],
    x
  )
  str_section_head("Scaled") # logging
}

if (length(samples$name) == 1) {
  save_object(x, "individual")
} else { # integrate
  # create and save list of seurat objects
  objects <- list()
  for (i in samples$name) {
    objects <- c(
      objects,
      get(i) # need get() to call object instead of string
    )
  }

  save_object(objects, "individual")
}

print("End of create_object.R")
