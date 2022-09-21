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

## load H5 reference for cluster mapping
# download.file(
#  url = "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat",
#  destfile = paste0(out_path, "tmp/pbmc_multimodal.h5seurat"),
#  method = "wget",
#  quiet = TRUE
# )
#
# reference <- LoadH5Seurat(
#  paste0(
#    out_path,
#    "tmp/pbmc_multimodal.h5seurat"
#  )
# )
#
# reference <- ScaleData(
#  reference,
#  assay = "SCT"
# )
#
# reference <- RunSPCA(
#  reference,
#  assay = "SCT",
#  graph = "wsnn"
# )
#
# reference <- FindNeighbors(
#  object = reference,
#  reduction = "spca",
#  dims = 1:50,
#  graph.name = "spca.annoy.neighbors",
#  k.param = 50,
#  cache.index = TRUE,
#  return.neighbor = TRUE,
#  l2.norm = TRUE
# )
#
# SaveAnnoyIndex(
#  object = reference[["spca.annoy.neighbors"]],
#  file = paste0(
#    out_path,
#    "tmp/preamble_image.RData"
#  )
# )
#
# reference[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(
#  object = reference[["spca.annoy.neighbors"]],
#  file = paste0(
#    out_path,
#    "tmp/preamble_image.RData"
#  )
# )
#
# anchors <- list()

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
  cellranger_peaks <- x$Peaks

  # create Seurat object
  x <- CreateSeuratObject(
    counts = x$`Gene Expression`,
    meta.data = metadata
  )
  str_section_head("Base Seurat Object")

  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")

  # add in the ATAC-seq data
  grange.counts <- StringToGRanges(
    rownames(cellranger_peaks),
    sep = c(":", "-")
  )
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  cellranger_peaks <- cellranger_peaks[as.vector(grange.use), ]
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"

  x[["ATAC"]] <- CreateChromatinAssay(
    counts = cellranger_peaks,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = paste0(
      samples$dir[i],
      "/atac_fragments.tsv.gz"
    ),
    min.cells = 10,
    annotation = annotation
  )

  DefaultAssay(x) <- "ATAC"
  str_section_head("with ATAC")

  p1 <- VlnPlot(x,
    features = c(
      "nCount_ATAC",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3,
    log = TRUE
  ) + NoLegend()

  save_figure(
    p1,
    paste0(name, "_counts_vln")
  )

  # more accurate peak calling using MACS2
  peaks <- CallPeaks(
    x,
    macs2.path = "/home/cho_lab/chris/miniconda3/envs/signac/bin/macs2"
  )

  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(
    peaks,
    pruning.mode = "coarse"
  )

  peaks <- subsetByOverlaps(
    x = peaks,
    ranges = blacklist_hg38_unified,
    invert = TRUE
  )

  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(x),
    features = peaks,
    cells = colnames(x)
  )

  # create a new assay using the MACS2 peak set and add it to the Seurat object
  x[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = paste0(
      samples$dir[i],
      "/atac_fragments.tsv.gz"
    ),
    annotation = annotation
  )

  DefaultAssay(x) <- "peaks"

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
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    ncol = 3
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
      percent.mt < params["max.percent.mt", ]
  )
  str_section_head("Filtered")

  # RNA analysis
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(
    x,
    vst.flavor = "v2",
    verbose = FALSE
  ) %>%
    RunPCA() %>%
    RunUMAP(
      dims = 1:d,
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_"
    )

  #  # determine anchors between reference and query for mapping
  #  anchors[[i]] <- FindTransferAnchors(
  #    reference = reference,
  #    query = x,
  #    normalization.method = "SCT",
  #    reference.reduction = "spca",
  #    reference.neighbors = "spca.annoy.neighbors",
  #    dims = 1:50
  #  )
  #
  #  x <- MapQuery(
  #    anchorset = anchors[[i]],
  #    query = x,
  #    reference = reference,
  #    refdata = list(
  #      celltype = "celltype.l2",
  #      predicted_ADT = "ADT"
  #    ),
  #    reference.reduction = "spca",
  #    reduction.model = "wnn.umap"
  #  )
  #  str_section_head("WNN Mapped")
  #
  #  # merge reference and query
  #  reference$id <- "reference"
  #  x$id <- "query"
  #  refquery <- merge(
  #    reference,
  #    x
  #  )
  #
  #  refquery[["spca"]] <- merge(
  #    reference[["spca"]],
  #    x[["ref.spca"]]
  #  )
  #
  #  refquery <- RunUMAP(
  #    refquery,
  #    reduction = "spca",
  #    dims = 1:50
  #  )
  #
  #  p1 <- DimPlot(
  #    refquery,
  #    group.by = "id",
  #    shuffle = TRUE
  #  )
  #
  #  save_figure(
  #    p1,
  #    paste0(name, "_mapping_dim")
  #  )

  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(x) <- "peaks"
  x <- RunTFIDF(x) %>%
    FindTopFeatures(
      min.cutoff = "q0"
    ) %>%
    RunSVD()

  x <- RegionStats(
    x,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )

  p1 <- DepthCor(x)

  save_figure(
    p1,
    paste0(name, "_macs2_depth")
  )

  x <- RunUMAP(
    x,
    reduction = "lsi",
    dims = 2:d,
    reduction.name = "umap.atac",
    reduction.key = "atacUMAP_"
  ) %>%
    FindMultiModalNeighbors(
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:d, 2:d)
    ) %>%
    RunUMAP(
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_"
    ) %>%
    FindClusters(
      graph.name = "wsnn",
      algorithm = 3,
    )

  #  clust_idents <- na.omit(clusters[, i])
  #  names(clust_idents) <- levels(x)
  #  x <- RenameIdents(x, clust_idents)
  #  x$celltype <- Idents(x)

  #  p1 <- DimPlot(
  #    x,
  #    reduction = "umap.rna",
  #    group.by = "predicted.celltype",
  #    label = TRUE,
  #    label.size = 2.5,
  #    repel = TRUE
  #  ) +
  #    ggtitle("RNA") +
  #    NoLegend()
  #
  #  p2 <- DimPlot(
  #    x,
  #    reduction = "umap.atac",
  #    group.by = "predicted.celltype",
  #    label = TRUE,
  #    label.size = 2.5,
  #    repel = TRUE
  #  ) +
  #    ggtitle("ATAC") +
  #    NoLegend()
  #
  #  p3 <- DimPlot(
  #    x,
  #    reduction = "wnn.umap",
  #    group.by = "predicted.celltype",
  #    label = TRUE,
  #    label.size = 2.5,
  #    repel = TRUE
  #  ) +
  #    ggtitle("WNN") +
  #    NoLegend()
  #
  #  save_figure(
  #    p1 + p2 + p3,
  #    paste0(name, "_clustered"),
  #    width = 18,
  #    height = 6
  #  )
  str_section_head("Clustered")

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
    assay = "peaks",
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

  # finding co-accessible networks with Cicero

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
}

# integration
# https://satijalab.org/signac/1.2.0/articles/integration.html
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#perform-integration-using-pearson-residuals-1
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

DefaultAssay(x) <- "peaks"
x <- RunTFIDF(x) %>%
  FindTopFeatures(
    min.cutoff = "q0"
  ) %>%
  RunSVD() %>%
  RunUMAP(
    reduction = "lsi",
    dims = 2:d,
    reduction.name = "umap.atac"
  )

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
  reduction = "lsi",
  assay.use = "peaks",
  project.dim = FALSE,
  reduction.save = "harmony_p"
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
    "harmony_p"
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

print("End of atac-multi.wnn.R")
