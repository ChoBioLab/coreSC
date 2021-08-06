# !/usr/bin/env Rscript

################################################################################
# PREAMBLE
################################################################################

# capture bash vars
args <- commandArgs(trailingOnly = T)

out_path <- paste0("./output/", args[1], "/")
params <- read.csv("./config/params.csv", row.names = 1)
samples <- read.csv("./config/samples.csv")

# package install check and load
packages <- c(
  "Seurat",
  "tidyverse",
  "tidyseurat",
  "cowplot",
  "patchwork"
)

package_check <- lapply(
  packages,
  function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# convenience functions
save_figure <- function(plots, name, type = "png", width, height, res) {
  if (type == "png") {
    png(paste0(out_path, name, ".", type),
      width = width, height = height, units = "in", res = 200
    )
  } else {
    pdf(paste0(out_path, name, ".", type),
      width = width, height = height
    )
  }
  print(plots)
  dev.off()
}

save_object <- function(object, name) {
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

################################################################################
# OBJECT CREATION
################################################################################

for (i in 1:nrow(samples)) {
  x <- Read10X(
    data.dir = samples$dir[i]
  ) %>%
    CreateSeuratObject(
      project = samples$project[i],
      min.cells = params["min.cells", ],
      min.features = params["min_nFeatures", ]
    )
  x[["percent.mt"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^MT-"
  )
  plot <- VlnPlot(
    x,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3
  )
  plot1 <- FeatureScatter(
    x,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )
  plot2 <- FeatureScatter(
    x,
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
  x <- subset(
    x,
    nCount_RNA > params["min_nCount", ] &
      nCount_RNA < params["max_nCount", ] &
      percent.mt < params["percent.mt", ]
  ) %>%
    NormalizeData()
  all.set <- rownames(x)
  x <- ScaleData(
    x,
    verbose = F,
    features = all.set
  ) %>%
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = params["max_nFeatures", ]
    )
  x@meta.data$object <- samples$name[i]
  x@meta.data$group <- samples$group[i]
  assign(
    samples$name[i],
    x
  )
}

# create list of seurat objects
objects <- list()
for (i in samples$name) {
  objects <- c(objects, get(i))
}

################################################################################
# INTEGRATION AND CLUSTERING
################################################################################

# sample integration
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
  all.genes <- rownames(x)
  integrated <- ScaleData(
    x,
    verbose = FALSE,
    features = all.genes
  ) %>%
    RunPCA(
      npcs = d,
      verbose = FALSE
    ) %>%
    RunUMAP(
      reduction = "pca",
      dims = 1:d
    ) %>%
    FindNeighbors(
      reduction = "pca",
      dims = 1:d
    ) %>%
    FindClusters(
      resolution = 0.6
    )
}

################################################################################
# VISUALIZATION
################################################################################

## visualzing linear dim reduction
# for (i in objects) {
# VizDimLoadings(i, dims = 1:2, reduction = "pca")
# DimPlot(i, reduction = "pca")
# DimHeatmap(i, dims = 1:15, cells = 500, balanced = TRUE)
# JackStrawPlot(i, dims = 1:15)
# ElbowPlot(i)
# }
#
#
## visualzing non-linear dim reduction
# dimheat <- DimHeatmap(
#  integrated,
#  dims = 1:params["dims", ],
#  cells = 500,
#  balanced = TRUE
# )
# JackStrawPlot(integrated, dims = 1:15)
# ElbowPlot(integrated)
#
# save_figure(
#            dimheat,
#            "integrated_dimheat_red",
#            width = 12,
#            height = 6
#            )

plot1 <- DimPlot(
  integrated,
  reduction = "umap",
  group.by = "group"
)

plot2 <- DimPlot(
  integrated,
  reduction = "umap",
  label = T
)

save_figure(
  (plot1 + plot2),
  "integrated_dimplot_red",
  width = 12,
  height = 6
)

saveRDS(integrated, file = paste0(out_path, "integrated_output.rds"))
# save.image(file = paste0("./output/", file_dir, "workspace.Rdata"))
