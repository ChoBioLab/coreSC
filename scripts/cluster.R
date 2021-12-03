# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(future) # parallelization

plan(multicore) # parallelization
options(future.globals.maxSize = 2000 * 1024^2)

load("./tmp/preamble_image.RData")
x <- read_object("individual")
d <- params["dims", ]

# run check for single sample
if (length(samples$name) == 1) {
  message("Single sample detected - skipping integration")
} else { # integrate
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
}

all.genes <- rownames(x)

x <- ScaleData(
  x,
  verbose = FALSE,
  features = all.genes
)

xPCA <- RunPCA(
  x,
  npcs = d,
  verbose = FALSE
)

saveRDS(xPCA, file = "./tmp/xPCA.RDS") # saving for dimensionality

x <- RunUMAP(
  xPCA,
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
