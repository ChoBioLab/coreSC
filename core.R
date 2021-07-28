# !/usr/bin/env Rscript

# install.packages("Seurat","dplyr","cowplot")
library(Seurat)
library(tidyverse)
#library(cowplot)
library(patchwork)

params <- read.csv("./config/params.csv", row.names = 1)
samples <- read.csv("./config/samples.csv")

# Convenience functions
#https://support.parsebiosciences.com/hc/en-us/articles/360053078092-Seurat-Tutorial-65k-PBMCs
save_figure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
      width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
      width = width, height = height)
}
print(plots)
dev.off()
}

save_object <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

read_object <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}



# iterate using is the method
samples$dir[1]



# The integrate function. Change the dimensions as per need####
integrate <- function(data) {
  imm.anchs <- FindIntegrationAnchors(object.list = data, dims = 1:13)
  imm.int <- IntegrateData(anchorset = imm.anchs, dims = 1:13)
  DefaultAssay(imm.int) <- "integrated"
  all.genes <- rownames(imm.int)
  imm.int <- ScaleData(imm.int, verbose = FALSE, features = all.genes)
  imm.int <- RunPCA(imm.int, npcs = 13, verbose = FALSE)
  imm.int <- RunUMAP(imm.int, reduction = "pca", dims = 1:13)
  imm.int <- FindNeighbors(imm.int, reduction = "pca", dims = 1:13)
  imm.int <- FindClusters(imm.int, resolution = 0.6)
  return(imm.int)
}

# Creating seurat objects from raw files. Change min # of RNA feature ###
seurat_object <- function(file_name) {
  object <- Read10X(data.dir = samples$dir)
  #colnames(x = object) <- paste(s, colnames(x = object), sep = "_")
  object <- CreateSeuratObject(object, project = samples$project, min.cells = params["min.cell",])
  object@meta.data$object <- file_name
  object <- subset(object, subset = nCount_RNA > params["min_nFeature_RNA",] & nCount_RNA < params["max_nFeature_RNA",])
  object <- NormalizeData(object)
  all.set <- rownames(object)
  object <- ScaleData(object, display.progress = F, features = all.set)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = params["nfeatures",])
  return(object)
}

#my_Dimplot <- function(data) {
#  p1 <- DimPlot(data, reduction = "umap", group.by = "group")
#  p2 <- DimPlot(data, reduction = "umap", label = TRUE)
#  p123 <- plot_grid(p1, p2)
#  png(filename = "Dim_PLot.png")
#  return(plot(p123))
#  dev.off()
#}
#
#my_dimheat <- function(data, markers) {
#  a2 <- DoHeatmap(object = data, features = markers)
#  png(filename = "heatmap.png")
#  return(plot(a2))
#  dev.off()
#}
#
#my_vlnplot <- function(data, markers) {
#  a3 <- VlnPlot(object = data, features = markers)
#  png(filename = "VlnplotCombinedall.png")
#  return(plot(a3))
#  dev.off()
#}
#
#my_tsne <- function(data) {
#  p1 <- TSNEPlot(data, pt.size = 0.5, group.by = "group")
#  p2 <- TSNEPlot(data, pt.size = 0.8)
#  p3 <- plot_grid(p1, p2)
#
#  png(filename = "TSNE_Combinedall.png")
#  return(plot(p3))
#  dev.off()
#}

##############################
# Script, change the PATH
##############################

#files <- list.files("/Volumes/Back_up2_MG/Naiyun/UC/Count")
#p.dir <- "/Volumes/Back_up2_MG/Naiyun/UC/Count"
#
## you can create group label as per the variable names below, this list is manual as of now ###
#new.group <- c("inflamed", "uninflamed", "inflamed", "uninflamed", "inflamed", "uninflamed", "inflamed", "uninflamed", "uninflamed")
#names <- files


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
Hall <- NULL
for (i in names) {
  Object <- get(paste0(i))
  Hall <- append(Object, Hall)
}

### define you gene names for the markers ###

#mymarkers <- c("", "") # fill your genes of interest

#immune.combinedall <- integrate(Hall)
# my_Dimplot(immune.combinedall)
# my_dimheat(immune.combinedall,mymarkers)
# my_vlnplot(immune.combinedall,mymarkersnewCluster<-readRDS("IntegratedCD_3_0.rds))
# my_tsne(immune.combinedall)


# save.image(file="MyData.RData")

saveRDS(immune.combinedall, file = "Integratedmodel_3_0.rds")
