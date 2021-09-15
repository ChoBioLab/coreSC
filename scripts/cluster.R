# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)

load("./tmp/base_image.RData")
x <- readRDS(paste0(out_path, "individual_objects.rds"))
d <- params["dims", ]

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

x <- JackStraw(x)

x <- ScoreJackStraw(
  x,
  dims = 1:d
)

plot1 <- JackStrawPlot(
  x,
  dims = 1:d
)

plot2 <- ElbowPlot(x)

save_figure(
  (plot1 + plot2),
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

save_object(x, "clustered")
print("End of cluster.R")
