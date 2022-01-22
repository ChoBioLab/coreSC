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
  str_section_head("Base Seurat Object") # logging

  x <- subset(
    x,
    nCount_RNA > params["min.count", ] &
      nCount_RNA < params["max.count", ] &
      percent.mt < params["max.percent.mt", ] &
      percent.mt > params["min.percent.mt", ]
  )

  x <- NormalizeData(x)
  str_section_head("Subset, Normalized") # logging

  genes <- rownames(x)

  x <- ScaleData(
    x,
    verbose = FALSE,
    features = genes
  )

  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures = params["max.features", ]
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
