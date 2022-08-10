library(Seurat)
library(Signac)
library(SeuratDisk)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(future)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")
load(paste0(out_path, "tmp/preamble_image.RData"))

plan(multicore) # parallelization
options(future.globals.maxSize = params["future.mem", ] * 1024^2)

# # load H5 reference for cluster mapping
# download.file(
#   url = "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat",
#   destfile = paste0(out_path, "tmp/pbmc_multimodal.h5seurat"),
#   method = "wget",
#   quiet = TRUE
# )

# reference <- LoadH5Seuerat(params["clust.ref", ])

for (i in 1:nrow(samples)) {
  name <- samples$name[i]
  # the 10x hdf5 file contains both data types.
  x <- Read10X_h5(
    paste0(
      samples$dir[i],
      "/filtered_feature_bc_matrix.h5"
    )
  )
  str_section_head("Raw Object")

  # metadata is found in the per_barcode_metrics.csv in cellranger-arc
  metadata <- read.csv(
    file = paste0(
      samples$dir[i],
      "/per_barcode_metrics.csv"
    ),
    header = TRUE,
    row.names = 1
  )

  # extract RNA and ATAC data
  rna_counts <- x$`Gene Expression`
  atac_counts <- x$Peaks

  # Create Seurat object
  x <- CreateSeuratObject(
    counts = rna_counts,
    meta.data = metadata
  )
  str_section_head("Base Seurat Object")

  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")

  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"

  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = paste0(
      samples$dir[i],
      "/atac_fragments.tsv.gz"
    ),
    min.cells = 10,
    annotation = annotations
  )

  x[["ATAC"]] <- chrom_assay
  DefaultAssay(x) <- "ATAC"
  str_section_head("with ATAC")

  p1 <- VlnPlot(x,
    features = c(
      "nCount_ATAC",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3,
    log = TRUE,
    pt.size = 0
  ) + NoLegend()

  save_figure(
    p1,
    paste0(name, "_counts_vln")
  )

  # compute nucleosome signal score per cell
  x <- NucleosomeSignal(object = x)

  # compute TSS enrichment score per cell
  x <- TSSEnrichment(
    object = x,
    fast = FALSE
  )

  # add blacklist ratio and fraction of reads in peaks
  x$pct_reads_in_peaks <- x$atac_peak_region_fragments / x$atac_fragments * 100
  x$blacklist_ratio <- x$blacklist_region_fragments / x$peak_region_fragments
  x$high.tss <- ifelse(
    x$TSS.enrichment > params["tss.score", ],
    "High",
    "Low"
  )
  str_section_head("Nucleosome and TSS Scores")

  p1 <- TSSPlot(
    x,
    group.by = "high.tss"
  ) + NoLegend()

  save_figure(
    p1,
    paste0(name, "_tss")
  )

  x$nucleosome_group <- ifelse(
    x$nucleosome_signal > 4,
    "NS > 4",
    "NS < 4"
  )

  p1 <- FragmentHistogram(
    object = x,
    group.by = "nucleosome_group"
  )

  save_figure(
    p1,
    paste0(name, "_frag_histogram")
  )

  p1 <- VlnPlot(
    object = x,
    features = c(
      "pct_reads_in_peaks",
      "peak_region_fragments",
      "TSS.enrichment",
      "blacklist_ratio",
      "nucleosome_signal"
    ),
    pt.size = 0.1,
    ncol = 5
  )

  save_figure(
    p1,
    paste0(name, "_atac_vln")
  )

  x <- subset(
    x = x,
    subset = nCount_ATAC < params["max.count.atac", ] &
      nCount_ATAC > params["min.count.atac", ] &
      nCount_RNA < params["max.count.rna", ] &
      nCount_RNA > params["min.count.rna", ] &
      pct_reads_in_peaks > params["pct.reads.peaks", ] &
      nucleosome_signal < params["nucleosome", ] &
      TSS.enrichment > params["tss.score", ] &
      percent.mt < params["max.percent.mt", ] &
      blacklist_ratio < 0.05
  )
  str_section_head("Filtered")

  # RNA analysis
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(
    x,
    method = "glmGamPoi",
    verbose = FALSE
  ) %>%
    RunPCA() %>%
    RunUMAP(
      dims = 1:50,
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_"
    )

  # determine anchors between reference and query for mapping
  anchors <- FindTransferAnchors(
    reference = reference,
    query = x,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )

  x <- MapQuery(
    anchorset = anchors,
    query = x,
    reference = reference,
    refdata = list(
      celltype = "celltype"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  str_section_head("WNN Mapped")

  # merge reference and query
  reference$id <- "reference"
  x$id <- "query"
  refquery <- merge(
    reference,
    x
  )

  refquery[["spca"]] <- merge(
    reference[["spca"]],
    x[["ref.spca"]]
  )

  refquery <- RunUMAP(
    refquery,
    reduction = "spca",
    dims = 1:50
  )

  p1 <- DimPlot(
    refquery,
    group.by = "id",
    shuffle = TRUE
  )

  save_figure(
    p1,
    paste0(name, "_mapping_dim")
  )

  save_H5object(refquery, "refquery_object")

  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(x) <- "ATAC"
  x <- RunTFIDF(x)
  x <- FindTopFeatures(
    x,
    min.cutoff = "q0"
  )

  x <- RunSVD(x)
  p1 <- DepthCor(x)

  save_figure(
    p1,
    paste0(name, "_depth")
  )

  x <- RunUMAP(
    x,
    reduction = "lsi",
    dims = 2:50,
    reduction.name = "umap.atac",
    reduction.key = "atacUMAP_"
  )

  x <- FindMultiModalNeighbors(
    x,
    reduction.list = list("pca", "lsi"),
    dims.list = list(1:50, 2:50)
  )

  x <- RunUMAP(
    x,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_"
  )

  x <- FindClusters(
    x,
    graph.name = "wsnn",
    algorithm = 3,
    verbose = FALSE
  )

  #  clust_idents <- na.omit(clusters[, i])
  #  names(clust_idents) <- levels(x)
  #  x <- RenameIdents(x, clust_idents)
  #  x$celltype <- Idents(x)

  p1 <- DimPlot(
    x,
    reduction = "umap.rna",
    group.by = "celltype",
    label = TRUE,
    label.size = 2.5,
    repel = TRUE
  ) +
    ggtitle("RNA") +
    NoLegend()

  p2 <- DimPlot(
    x,
    reduction = "umap.atac",
    group.by = "celltype",
    label = TRUE,
    label.size = 2.5,
    repel = TRUE
  ) +
    ggtitle("ATAC") +
    NoLegend()

  p3 <- DimPlot(
    x,
    reduction = "wnn.umap",
    group.by = "celltype",
    label = TRUE,
    label.size = 2.5,
    repel = TRUE
  ) +
    ggtitle("WNN") +
    NoLegend()

  save_figure(
    p1 + p2 + p3,
    paste0(name, "_clustered"),
    width = 18,
    height = 6
  )
  str_section_head("Clustered")

  ## to make the visualization easier, subset T cell clusters
  celltype.names <- levels(x)
  tcell.names <- grep(
    "CD4|CD8|Treg",
    celltype.names,
    value = TRUE
  )

  tcells <- subset(
    x,
    idents = tcell.names
  )

  p1 <- CoveragePlot(
    tcells,
    region = "CD8A",
    features = "CD8A",
    assay = "ATAC",
    expression.assay = "SCT",
    peaks = FALSE
  )

  save_figure(
    p1,
    paste0(name, "_coverage")
  )

  # Get a list of motif position frequency matrices from the JASPAR database
  pwm_set <- getMatrixSet(
    x = JASPAR2020,
    opts = list(
      species = 9606,
      all_versions = FALSE
    )
  )

  motif.matrix <- CreateMotifMatrix(
    features = granges(x),
    pwm = pwm_set,
    genome = "hg38",
    use.counts = FALSE
  )

  motif.object <- CreateMotifObject(
    data = motif.matrix,
    pwm = pwm_set
  )

  x <- SetAssayData(
    x,
    assay = "ATAC",
    slot = "motifs",
    new.data = motif.object
  )

  # Note that this step can take 30-60 minutes
  x <- RunChromVAR(
    object = x,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  str_section_head("Motif Scored")

  warnings()

  x@meta.data$object <- name
  x@meta.data$group <- samples$group[i]

  assign( # giving names to objects
    name,
    x
  )
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
  save_H5object(objects, "individual")
}

## integration
## https://satijalab.org/signac/articles/integrate_atac.html
# combined <- Reduce(merge, objects)
#
# combined <- FindTopFeatures(
#  combined,
#  min.cutoff = 10
# ) %>%
#  RunTFIDF() %>%
#  RunSVD() %>%
#  RunUMAP(
#    .,
#    reduction = "lsi",
#    dims = 2:30
#  )
#
# anchors <- FindIntegrationAnchors(
#  object.list = objects,
#  anchor.features = rownames(objects[1]),
#  reduction = "rlsi",
#  dims = 2:30
# )
#
## integrate LSI embeddings
# integrated <- IntegrateEmbeddings(
#  anchorset = anchors,
#  reductions = combined[["lsi"]],
#  new.reduction.name = "integrated_lsi",
#  dims.to.integrate = 1:30
# )
#
## create a new UMAP using the integrated embeddings
# integrated <- RunUMAP(
#  integrated,
#  reduction = "integrated_lsi",
#  dims = 2:30
# )
#
# save_object(integrated, "integrated")
# save_H5object(integrated, "integrated")
#

sessionInfo()
