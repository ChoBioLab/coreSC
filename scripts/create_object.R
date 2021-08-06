# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)

load("./tmp/base_image.RData")

for (i in 1:nrow(samples)) {
  x <- Read10X(
    data.dir = samples$dir[i]
  )
  message(strrep("=", 80))
  message(paste(samples$name[i], "Raw Object"))
  message(strrep("=", 80))
  str(x)
  x <- CreateSeuratObject(
    counts = x,
    project = samples$project[i],
    min.cells = params["min.cells", ],
    min.features = params["min.features", ]
  )
  x[["percent.mt"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^MT-"
  )
  message(strrep("=", 80))
  message(paste(samples$name[i], "Base Seurat Object"))
  message(strrep("=", 80))
  str(x)
  x <- subset(
    x,
    nCount_RNA > params["min.count", ] &
      nCount_RNA < params["max.count", ] &
      percent.mt < params["percent.mt", ]
  )
  x <- NormalizeData(x)
  message(strrep("=", 80))
  message(paste(samples$name[i], "Subset, Normalized"))
  message(strrep("=", 80))
  str(x)
  genes <- rownames(x)
  x <- ScaleData(
    x,
    verbose = F,
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
  message(strrep("=", 80))
  message(paste(samples$name[i], "Scaled"))
  message(strrep("=", 80))
  str(x)
}

# create and save list of seurat objects
objects <- list()
for (i in samples$name) {
  objects <- c(
    objects,
    get(i)
  )
  save_object(
    get(i),
    paste0(i, "_normed_scaled")
  )
}

saveRDS(objects, "./tmp/objects.RDS")
