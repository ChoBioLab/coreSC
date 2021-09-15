# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)

load("./tmp/base_image.RData")

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
      percent.mt < params["percent.mt", ]
  )
  x <- NormalizeData(x)
  str_section_head("Subset, Normalized") # logging
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
  str_section_head("Scaled") # logging
}

# run check for single sample
if (length(samples$name == 1)) {
  message("Single sample detected - skipping integration")
  saveRDS(x, file = paste0(out_path, "individual_objects.rds"))
} else {

  # create and save list of seurat objects
  objects <- list()
  for (i in samples$name) {
    objects <- c(
      objects,
      get(i) # need get() to call object instead of string
    )
  }

  x <- FindIntegrationAnchors(
    object.list = objects,
    dims = 1:d
  )

  x <- IntegrateData(
    anchorset = x,
    dims = 1:d
  )

  DefaultAssay(x) <- "integrated"
  saveRDS(x, file = paste0(out_path, "individual_objects.rds"))

  str_section_noloop("Integrated") # logging
}

print("End of create_object.R")
