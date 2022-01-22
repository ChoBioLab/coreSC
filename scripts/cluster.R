# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(future) # parallelization

plan(multicore) # parallelization
options(future.globals.maxSize = 2000 * 1024^2)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")

load(paste0(out_path, "tmp/preamble_image.RData"))
x <- read_object("individual")

message("object check")
str(x)
# run check for single sample
if (length(samples$name) == 1) {
  message("Single sample detected - skipping integration")
} else { # integrate
  features <- SelectIntegrationFeatures(
    object.list = objects
  )

  x <- FindIntegrationAnchors(
    object.list = objects,
    dims = 1:d,
    anchor.features = features
  )

  x <- IntegrateData(
    anchorset = x,
    dims = 1:d
  )

  DefaultAssay(x) <- "integrated"
  str_section_noloop("Integrated") # logging

  genes <- rownames(x)
  x <- ScaleData(
    x,
    verbose = FALSE,
    features = genes
  )
}

x <- RunPCA(
  x,
  npcs = d,
  verbose = FALSE
)

x <- JackStraw(
  x,
  num.replicate = 100
)

x <- ScoreJackStraw(
  x,
  dims = 1:d
)

p1 <- JackStrawPlot(
  x,
  dims = 1:d
)

p2 <- ElbowPlot(x)

save_figure(
  (p1 + p2),
  "dimensionality",
  width = 12,
  height = 6
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

p1 <- DimPlot(
  x,
  reduction = "umap",
  group.by = "group"
)

p2 <- DimPlot(
  x,
  reduction = "umap",
  label = T
)

save_figure(
  (p1 + p2),
  "combined_dimplot_red",
  width = 12,
  height = 6
)

save_object(x, "clustered")
print("End of cluster.R")
