# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)

load("./tmp/base_image.RData")
objects <- read_object("individual_objects")
combined <- read_object("combined_integrated")

# individual sample qc
for (i in 1:nrow(samples)) {
  plot <- VlnPlot(
    objects[[i]],
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3
  )
  plot1 <- FeatureScatter(
    objects[[i]],
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )
  plot2 <- FeatureScatter(
    objects[[i]],
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
  combined,
  reduction = "umap",
  group.by = "group"
)

plot2 <- DimPlot(
  combined,
  reduction = "umap",
  label = T
)

save_figure(
  (plot1 + plot2),
  "combined_dimplot_red",
  width = 12,
  height = 6
)

print("End of visualize.R")
