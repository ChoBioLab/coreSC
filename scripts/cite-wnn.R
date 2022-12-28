# !/usr/bin/env Rscript

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(ggplot2)
library(future) # parallelization
library(harmony)
library(limma)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")
load(paste0(out_path, "tmp/preamble_image.RData"))

plan(
  multicore,
  workers = params["future.workers", ]
) # parallelization
options(future.globals.maxSize = params["future.mem", ] * 1024^2)

for (i in 1:nrow(samples)) {
  name <- samples$name[i]
  raw <- Read10X( # pulling data with no filter
    data.dir = samples$dir[i]
  )

  rna <- raw$`Gene Expression`
  adt <- raw$`Antibody Capture`

  x <- CreateSeuratObject(
    counts = rna,
    project = samples$project[i],
    min.cells = params["min.cells", ],
    min.features = params["min.features", ]
  )

  x@meta.data$object <- samples$name[i]
  x@meta.data$group <- samples$group[i]

  x[["ADT"]] <- CreateAssayObject(counts = adt)

  x[["percent.mt"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^MT-"
  )

  x[["percent.rb"]] <- PercentageFeatureSet(
    x,
    pattern = "RPS|RPL"
  )

  x[["percent.hb"]] <- PercentageFeatureSet(
    x,
    pattern = "(?i)^HB[^(P)]"
  )


  # raw qc plots
  p1 <- VlnPlot(
    x,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "percent.rb",
      "percent.hb"
    ),
    ncol = 5
  )

  save_figure(
    p1,
    paste0(samples$name[i], "_rna_unfilt_vln"),
    width = 12,
    height = 6
  )

  p1 <- VlnPlot(
    x,
    features = c(
      "nFeature_ADT",
      "nCount_ADT"
    ),
    ncol = 2
  )

  save_figure(
    p1,
    paste0(samples$name[i], "_adt_unfilt_vln"),
    width = 12,
    height = 6
  )

  p1 <- FeatureScatter(
    x,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )

  p2 <- FeatureScatter(
    x,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

  p3 <- FeatureScatter(
    x,
    feature1 = "nCount_ADT",
    feature2 = "nFeature_ADT"
  )

  save_figure(
    (p1 + p2),
    paste0(samples$name[i], "_unfilt_scatter"),
    width = 12,
    height = 6
  )

  str_section_head("Base Seurat Object") # logging

  DefaultAssay(x) <- "RNA"
  x <- subset(
    x,
    nFeature_RNA > params["min.features", ] &
      nFeature_RNA < params["max.features", ] &
      nCount_RNA > params["min.count.rna", ] &
      nCount_RNA < params["max.count.rna", ] &
      percent.mt < params["max.percent.mt", ] &
      percent.mt > params["min.percent.mt", ]
  )

  # norm, dimred, and clustering
  # clustering is performed on individual samples for RNA QC
  x <- SCTransform(
    x,
    vst.flavor = "v2",
    verbose = FALSE
  ) %>%
    RunPCA(
      npcs = d,
      verbose = FALSE
    ) %>%
    RunUMAP(
      dims = 1:d,
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_"
    )

  # norm, neighborhood formation and clustering
  DefaultAssay(x) <- "ADT"
  VariableFeatures(x) <- rownames(x[["ADT"]])

  x <- NormalizeData(
    x,
    normalization.method = "CLR",
    margin = 2
  ) %>%
    ScaleData() %>%
    RunPCA(reduction.name = "apca") %>%
    FindMultiModalNeighbors(
      reduction.list = list(
        "pca",
        "apca"
      ),
      dims.list = list(
        1:d,
        1:18
      ),
      modality.weight.name = "RNA.weight"
    ) %>%
    RunUMAP(
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_"
    ) %>%
    FindClusters(
      graph.name = "wsnn",
      algorithm = 3,
      resolution = params["res", ],
      verbose = FALSE
    )

  str_section_head("Subset, SCT Normd, Red, Clust") # logging

  # clustering by modality
  x <- RunUMAP(
    x,
    reduction = "pca",
    dims = 1:d,
    assay = "RNA",
    reduction.name = "rna.umap",
    reduction.key = "rnaUMAP_"
  )

  x <- RunUMAP(
    x,
    reduction = "apca",
    dims = 1:d,
    assay = "ADT",
    reduction.name = "adt.umap",
    reduction.key = "adtUMAP_"
  )

  p1 <- DimPlot(
    x,
    reduction = "rna.umap",
    label = TRUE,
    repel = TRUE,
    label.size = 2.5
  ) + NoLegend() +
    ggtitle("RNA")

  p2 <- DimPlot(
    x,
    reduction = "adt.umap",
    label = TRUE,
    repel = TRUE,
    label.size = 2.5
  ) + NoLegend() +
    ggtitle("ADT")

  p3 <- DimPlot(
    x,
    reduction = "wnn.umap",
    label = TRUE,
    repel = TRUE,
    label.size = 2.5
  ) + NoLegend() +
    ggtitle("WNN")

  save_figure(
    p1 + p2 + p3,
    paste0(samples$name[i], "_clustered"),
    width = 18,
    height = 6
  )

  # plotting clusters by specified group
  p1 <- DimPlot(
    x,
    reduction = "wnn.umap",
    group.by = "group",
    label = T
  )

  p2 <- DimPlot(
    x,
    reduction = "wnn.umap",
    label = T
  )

  save_figure(
    (p1 + p2),
    paste0(samples$name[i], "_individual_dimplot"),
    width = 12,
    height = 6
  )

  assign( # giving names to objects
    samples$name[i],
    x
  )
}

if (length(samples$name) == 1) {
  save_object(x, "individual_clustered")
} else { # integrate
  # create and save list of seurat objects
  objects <- list()
  for (i in samples$name) {
    objects <- c(
      objects,
      get(i) # need get() to call object instead of string
    )
  }

  save_object(objects, "individual_clustered")
}

# integration
# https://satijalab.org/signac/1.2.0/articles/integration.html
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#perform-integration-using-pearson-residuals-1
# https://github.com/satijalab/seurat/issues/6089
x <- Reduce(merge, objects)
save_h5(x, "combined")

str_section_noloop("Combined")

DefaultAssay(x) <- "RNA"
x <- SCTransform(
  x,
  vst.flavor = "v2",
  verbose = FALSE
) %>%
  RunPCA(
    npcs = d,
    verbose = FALSE
  ) %>%
  RunUMAP(
    dims = 1:d,
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_"
  )

DefaultAssay(x) <- "ADT"
VariableFeatures(x) <- rownames(x[["ADT"]])

x <- NormalizeData(
  x,
  normalization.method = "CLR",
  margin = 2
) %>%
  ScaleData() %>%
  RunPCA(reduction.name = "apca")

p1 <- DimPlot(
  x,
  group.by = "object"
)

save_figure(
  p1,
  "combined_dimplot"
)

str_section_noloop("Pre-harmony Normalized")

# harmonize
x <- RunHarmony(
  object = x,
  group.by.vars = "object",
  reduction = "apca",
  assay.use = "ADT",
  project.dim = FALSE,
  reduction.save = "harmony_a"
)

x <- RunHarmony(
  object = x,
  group.by.vars = "object",
  reduction = "pca",
  assay.use = "SCT",
  project.dim = FALSE,
  reduction.save = "harmony_r"
)

str_section_noloop("Harmonized")

# integrated multimodal wnn
x <- FindMultiModalNeighbors(
  object = x,
  reduction.list = list(
    "harmony_r",
    "harmony_a"
  ),
  dims.list = list(
    1:d,
    1:18
  ),
  modality.weight.name = "RNA.weight"
) %>%
  RunUMAP(
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    assay = "SCT"
  ) %>%
  FindClusters(
    graph.name = "wsnn",
    algorithm = 3,
    resolution = params["res", ]
  )

str_section_noloop("Integrated")
save_h5(x, "integrated")

p1 <- DimPlot(
  x,
  reduction = "wnn.umap"
)

save_figure(
  p1,
  "integrated_dimplot"
)

DefaultAssay(x) <- "SCT"
x <- PrepSCTFindMarkers(x)

markers <- FindAllMarkers(
  x,
  assay = "SCT",
  verbose = FALSE
)

write.csv(markers, paste0(out_path, "all_markers.csv"))

sessionInfo()

print("End of cite-wnn.R")
