library(Seurat)
library(Signac)
library(SeuratDisk)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(future)

plan(multicore) # parallelization
options(future.globals.maxSize = 3000 * 1024^2)

load("./tmp/preamble_image.RData")

# load H5 reference for cluster mapping
download.file(
  url = "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat",
  destfile = "./tmp/pbmc_multimodal.h5seurat",
  method = "wget"
)

reference <- LoadH5Seurat("./tmp/pbmc_multimodal.h5seurat")

for (i in 1:nrow(samples)) {
  print(paste0(samples$dir[i], "/filtered_feature_bc_matrix.h5"))
  # the 10x hdf5 file contains both data types.
  inputdata.10x <- Read10X_h5(
    paste0(
      samples$dir[i],
      "/filtered_feature_bc_matrix.h5"
    )
  )

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
  rna_counts <- inputdata.10x$`Gene Expression`
  atac_counts <- inputdata.10x$Peaks

  # Create Seurat object
  x <- CreateSeuratObject(
    counts = rna_counts,
    meta.data = metadata
  )

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
    paste0(samples$name[i], "_counts_vln"),
    width = 12,
    height = 6
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
  # x$blacklist_ratio <- x$blacklist_region_fragments / x$peak_region_fragments
  x$high.tss <- ifelse(
    x$TSS.enrichment > 2,
    "High",
    "Low"
  )

  p1 <- TSSPlot(
    x,
    group.by = "high.tss"
  ) + NoLegend()

  save_figure(
    p1,
    paste0(samples$name[i], "_tss"),
    width = 12,
    height = 6
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
    paste0(samples$name[i], "_frag_histogram"),
    width = 12,
    height = 6
  )

  p1 <- VlnPlot(
    object = x,
    features = c(
      "pct_reads_in_peaks",
      "peak_region_fragments",
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    pt.size = 0.1,
    ncol = 5
  )

  save_figure(
    p1,
    paste0(samples$name[i], "_atac_vln"),
    width = 12,
    height = 6
  )

  x <- subset(
    x = x,
    subset = nCount_ATAC < 7e4 &
      nCount_ATAC > 5e3 &
      nCount_RNA < 25000 &
      nCount_RNA > 1000 &
      pct_reads_in_peaks > 15 &
      nucleosome_signal < 2 &
      TSS.enrichment > 1 &
      percent.mt < 20
  )

  # RNA analysis
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(
    x,
    verbose = FALSE
  ) %>%
    RunPCA() %>%
    RunUMAP(
      dims = 1:50,
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_"
    )

  # determine anchors between reference and query
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
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )

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
    paste0(samples$name[i], "_mapping_dim"),
    width = 12,
    height = 6
  )

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
    paste0(samples$name[i], "_depth"),
    width = 12,
    height = 6
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

  p1 <- DimPlot(
    x,
    reduction = "umap.rna",
    group.by = "predicted.celltype.l2",
    label = TRUE,
    label.size = 2.5,
    repel = TRUE
  ) + ggtitle("RNA")

  p2 <- DimPlot(
    x,
    reduction = "umap.atac",
    group.by = "predicted.celltype.l2",
    label = TRUE,
    label.size = 2.5,
    repel = TRUE
  ) + ggtitle("ATAC")

  p3 <- DimPlot(
    x,
    reduction = "wnn.umap",
    group.by = "predicted.celltype.l2",
    label = TRUE,
    label.size = 2.5,
    repel = TRUE
  ) + ggtitle("WNN")

  save_figure(
    p1 + p2 + p3,
    paste0(samples$name[i], "_clustered"),
    width = 12,
    height = 6
  )

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
    paste0(samples$name[i], "_coverage"),
    width = 12,
    height = 6
  )


  assign( # giving names to objects
    samples$name[i],
    x
  )
}

if (length(samples$name) == 1) {
  save_object(x, file = "individual")
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
