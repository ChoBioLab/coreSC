# !/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(future) # parallelization
library(harmony)
library(limma)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")
load(paste0(out_path, "tmp/preamble_image.RData"))
objects <- read_object("individual_clustered")

plan(
  multicore,
  workers = params["future.workers", ]
) # parallelization
options(future.globals.maxSize = params["future.mem", ] * 1024^2)

# ref: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#perform-integration-using-pearson-residuals-1
# run check for single sample
if (length(samples$name) == 1) {
  message("Single sample detected - skipping integration")
  x <- objects[[1]]
} else { # integrate
  x <- RunHarmony(
    object = objects,
    group.by.vars = "object",
    reduction = "pca",
    assay.use = "SCT",
    reduction.save = "harmony"
  )

  DefaultAssay(x) <- "integrated"
  str_section_noloop("Integrated") # logging
}

# clustering
RunUMAP(
  x,
  reduction = "pca",
  dims = 1:d
) %>%
  FindNeighbors(
    reduction = "pca",
    dims = 1:d
  ) %>%
  FindClusters(
    resolution = params["res", ]
  )

str_section_noloop("Reduced & Clustered") # logging

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

x <- PrepSCTFindMarkers(x)

y <- FindAllMarkers(
  x,
  assay = "SCT",
  verbose = FALSE
)

save_h5(x, "integrated")
write.csv(y, paste0(out_path, "all_markers.csv"))

print("End of integrated.R")

sessionInfo()