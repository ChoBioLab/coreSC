# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")

load(paste0(out_path, "tmp/preamble_image.RData"))
objects <- read_object("individual")
xPCA <- readRDS(paste0(out_path, "tmp/xPCA.RDS"))
clustered <- read_object("clustered")

message("object check 2")
str(objects)
# individual sample qc
for (i in 1:nrow(samples)) {
  plot <- VlnPlot(
    objects[i],
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3
  )

  plot1 <- FeatureScatter(
    objects[i],
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )

  plot2 <- FeatureScatter(
    objects[i],
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

  save_figure(
    plot,
    paste0(samples$name[i], "_unfilt_vln"),
    width = 12,
    height = 6
  )

  save_figure(
    (plot1 + plot2),
    paste0(samples$name[i], "_unfilt_scatter"),
    width = 12,
    height = 6
  )
}

xPCA <- JackStraw(xPCA)
xPCA <- ScoreJackStraw(
  xPCA,
  dims = 1:d
)

jack <- JackStrawPlot(
  xPCA,
  dims = 1:d
)

elbow <- ElbowPlot(xPCA)

save_figure(
  (jack + elbow),
  "dimensionality",
  width = 12,
  height = 6
)

dim1 <- DimPlot(
  clustered,
  reduction = "umap",
  group.by = "group"
)

dim2 <- DimPlot(
  clustered,
  reduction = "umap",
  label = T
)

save_figure(
  (dim1 + dim2),
  "combined_dimplot_red",
  width = 12,
  height = 6
)

print("End of visualize.R")
