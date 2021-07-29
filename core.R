# !/usr/bin/env Rscript

################################################################################
# PREAMBLE
################################################################################

# package install check and load
packages <- c(
  "Seurat",
  "tidyverse",
  "tidyseurat",
  "cowplot",
  "patchwork"
)

# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
package_check <- lapply(
  packages,
  function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

params <- read.csv("./config/params.csv", row.names = 1)
samples <- read.csv("./config/samples.csv")

# convenience functions
# https://support.parsebiosciences.com/hc/en-us/articles/360053078092-Seurat-Tutorial-65k-PBMCs
save_figure <- function(plots, name, type = "png", width, height, res) {
  if (type == "png") {
    png(paste0(fig_path, name, ".", type),
      width = width, height = height, units = "in", res = 200
    )
  } else {
    pdf(paste0(fig_path, name, ".", type),
      width = width, height = height
    )
  }
  print(plots)
  dev.off()
}

save_object <- function(object, name) {
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

read_object <- function(name) {
  readRDS(paste0(data_path, name, ".RDS"))
}

# create list of seurat objects
objects <- list()
for (i in samples$name) {
  objects <- c(objects, get(i))
}

# object creation
for (i in 1:nrow(samples)) {
  x <- Read10X(
    data.dir = samples$dir[i]
  ) %>%
    CreateSeuratObject(
      project = samples$project[i],
      min.cells = params["min.cells", ]
    ) %>%
    subset(
      nCount_RNA > params["min_nCount", ] &
        nCount_RNA < params["max_nCount", ]
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
  assign(samples$name[i], x)
  testlist <- list(testlist, samples$name[i])
}

# sample integration
  for (i in 1:nrow(samples)) {
  d <- params["dims",]
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
  x <- ScaleData(
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

# my_Dimplot <- function(data) {
#  p1 <- DimPlot(data, reduction = "umap", group.by = "group")
#  p2 <- DimPlot(data, reduction = "umap", label = TRUE)
#  p123 <- plot_grid(p1, p2)
#  png(filename = "Dim_PLot.png")
#  return(plot(p123))
#  dev.off()
# }
#
# my_dimheat <- function(data, markers) {
#  a2 <- DoHeatmap(object = data, features = markers)
#  png(filename = "heatmap.png")
#  return(plot(a2))
#  dev.off()
# }
#
# my_vlnplot <- function(data, markers) {
#  a3 <- VlnPlot(object = data, features = markers)
#  png(filename = "VlnplotCombinedall.png")
#  return(plot(a3))
#  dev.off()
# }
#
# my_tsne <- function(data) {
#  p1 <- TSNEPlot(data, pt.size = 0.5, group.by = "group")
#  p2 <- TSNEPlot(data, pt.size = 0.8)
#  p3 <- plot_grid(p1, p2)
#
#  png(filename = "TSNE_Combinedall.png")
#  return(plot(p3))
#  dev.off()
# }

##############################
# Script, change the PATH
##############################

# files <- list.files("/Volumes/Back_up2_MG/Naiyun/UC/Count")
# p.dir <- "/Volumes/Back_up2_MG/Naiyun/UC/Count"
#
## you can create group label as per the variable names below, this list is manual as of now ###
# new.group <- c("inflamed", "uninflamed", "inflamed", "uninflamed", "inflamed", "uninflamed", "inflamed", "uninflamed", "uninflamed")
# names <- files


## Calling create seurat object function on each file ####
j <- 0
for (i in names) {{ j <- j + 1
  assign(i, seurat_object(i, p.dir, paste0("S", j), paste0("mygene.", j))) }}


### assigning the group names created above, comment out if not needed ####
gp <- 0
for (i in names) {
  gp <- gp + 1
  Object <- get(paste0(i))
  Object@meta.data$group <- new.group[gp]
  assign(i, Object)
  print(new.group[gp])
}


## creating the list for integartion ####
### define you gene names for the markers ###

# mymarkers <- c("", "") # fill your genes of interest

# immune.combinedall <- integrate(Hall)
# my_Dimplot(immune.combinedall)
# my_dimheat(immune.combinedall,mymarkers)
# my_vlnplot(immune.combinedall,mymarkersnewCluster<-readRDS("IntegratedCD_3_0.rds))
# my_tsne(immune.combinedall)


# save.image(file="MyData.RData")

saveRDS(immune.combinedall, file = "Integratedmodel_3_0.rds")
