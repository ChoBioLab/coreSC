# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)

load("./tmp/base_image.RData")
objects <- readRDS("./tmp/individual_objects")

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
str_section_noloop("Integrated") # logging
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
str_section_noloop("Reduced") # logging
x <- FindNeighbors(
  x,
  reduction = "pca",
  dims = 1:d
)
x <- FindClusters(
  x,
  resolution = params["res", ]
)
str_section_noloop("Clustered") # logging

save_object(x, "combined_integrated")

print("End of integrate.R")
