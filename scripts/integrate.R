# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)

load("./tmp/base_image.RData")
objects <- readRDS("./tmp/objects.RDS")

for (i in 1:nrow(samples)) {
  d <- params["dims", ]
  x <- FindIntegrationAnchors(
    object.list = objects,
    dims = 1:d
  )
  x <- IntegrateData(
    anchorset = x,
    dims = 1:d
  )
  DefaultAssay(x) <- "integrated"
  message(strrep("=", 80))
  message(paste(samples$name[i], "Integrated"))
  message(strrep("=", 80))
  str(x)
  all.genes <- rownames(x)
  x <- ScaleData(
    x,
    verbose = FALSE,
    features = all.genes
  )
  x <- RunPCA(
    x,
    npcs = d,
    verbose = FALSE
  )
  x <- RunUMAP(
    x,
    reduction = "pca",
    dims = 1:d
  )
  message(strrep("=", 80))
  message(paste(samples$name[i], "Reduced"))
  message(strrep("=", 80))
  str(x)
  x <- FindNeighbors(
    x,
    reduction = "pca",
    dims = 1:d
  )
  x <- FindClusters(
    x,
    resolution = 0.6
  )
  message(strrep("=", 80))
  message(paste(samples$name[i], "Clustered"))
  message(strrep("=", 80))
  str(x)
}

save_object(x, "combined_integrated")
