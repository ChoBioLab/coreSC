# Core libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)

# Data manipulation and visualization
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(future)

# Set up parallel processing
plan("multicore", workers = 24)
options(future.globals.maxSize = 32000 * 1024^2) # 50GB

# Define parameters (from your params.csv)
params <- list(
  min.cells = 3,
  min.count.rna = 500,
  max.count.rna = 60000,
  min.count.atac = 1000,
  max.count.atac = 70000,
  min.features = 200,
  max.features = 2000,
  max.percent.mt = 30,
  min.percent.mt = 0,
  pct.reads.peaks = 20,
  nucleosome = 2,
  tss.score = 1,
  dims = 17,
  res = 0.8,
  future.mem = 20000,
  future.workers = 20
)

print("Parameters loaded:")
print(params)

setwd('/mnt/data1/jzhang/projects/perianal-cd/analysis/coreSC')

# Create timestamped output directory
timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
out_path <- paste0('./output/output_', timestamp, '/')

# Create the directory if it doesn't exist
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

cat("Output directory created:", out_path, "\n")

sample_file <- "./config/sample.txt"
samples <- read.table(sample_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

print("Sample information loaded:")
print(paste("Number of samples:", nrow(samples)))

for (i in seq_len(nrow(samples))){
  name <- samples$sample_id[i]
  # the 10x hdf5 file contains both data types.
  x <- Read10X_h5(
    paste0(
      samples$atac_path[i],
      "/filtered_feature_bc_matrix.h5"
    )
  )
  str_section_head("Raw Object")
  # metadata is found in the per_barcode_metrics.csv in cellranger-arc
  metadata <- read.csv(
    file = paste0(
      samples$atac_path[i],
      "/per_barcode_metrics.csv"
    ),
    header = TRUE,
    row.names = 1
  )
  # extract RNA and ATAC data
  rna_counts <- x$`Gene Expression`
  cellranger_peaks <- x$Peaks
  # NEW STEP: Replace RNA counts with CellBender output if available
  if ("rna_path" %in% colnames(samples) && !is.na(samples$rna_path[i]) && samples$rna_path[i] != "") {
    cellbender_path <- samples$rna_path[i]    
    if (file.exists(cellbender_path)) {
      cat("Replacing RNA counts with CellBender output:", basename(cellbender_path), "\n")
      cellbender_counts <- NULL
      # Try standard Read10X_h5 first
    tryCatch({
      cellbender_counts <- Read10X_h5(cellbender_path)
    }, error = function(e) {
      cat("Read10X_h5 failed, trying direct HDF5 approach:", e$message, "\n")      
      # Try direct HDF5 reading
      tryCatch({
        library(rhdf5)
        # Matrix components
        data    <- h5read(cellbender_path, "/matrix/data")
        indices <- h5read(cellbender_path, "/matrix/indices")
        indptr  <- h5read(cellbender_path, "/matrix/indptr")
        shape   <- h5read(cellbender_path, "/matrix/shape")  # [n_genes, n_barcodes]
        cellbender_counts <- new("dgCMatrix",
          Dim = as.integer(shape),
          x = as.numeric(data),
          i = as.integer(indices),
          p = as.integer(indptr)
        )
        barcodes <- h5read(cellbender_path, "/matrix/barcodes")
        genes    <- h5read(cellbender_path, "/matrix/features/name")     # Gene names
        gene_ids <- h5read(cellbender_path, "/matrix/features/id")       # Gene IDs        
        rownames(cellbender_counts) <- genes
        colnames(cellbender_counts) <- barcodes
        H5close()
      }, error = function(e2) {
        warning(paste("Both Read10X_h5 and direct HDF5 failed:", e2$message))
      })
    })
    if (!is.null(cellbender_counts)) {    
        # Align barcodes between original multiome and CellBender
        original_barcodes <- colnames(rna_counts)
        cellbender_barcodes <- colnames(cellbender_counts)
        cat("Original RNA cells:", length(original_barcodes), "\n")
        cat("CellBender RNA cells:", length(cellbender_barcodes), "\n")
        # Find common barcodes
        common_barcodes <- intersect(original_barcodes, cellbender_barcodes)
        # Try barcode format adjustment if no common barcodes found
        if (length(common_barcodes) == 0) {
          if (all(grepl("-1$", original_barcodes)) && !all(grepl("-1$", cellbender_barcodes))) {
            # Original has -1 suffix, CellBender doesn't
            cellbender_barcodes_adj <- paste0(cellbender_barcodes, "-1")
            colnames(cellbender_counts) <- cellbender_barcodes_adj
            common_barcodes <- intersect(original_barcodes, cellbender_barcodes_adj)
            cat("Added -1 suffix to CellBender barcodes\n")
          } else if (!all(grepl("-1$", original_barcodes)) && all(grepl("-1$", cellbender_barcodes))) {
            # CellBender has -1 suffix, original doesn't
            original_barcodes_adj <- paste0(original_barcodes, "-1")
            # Update all relevant objects with new barcodes
            colnames(rna_counts) <- original_barcodes_adj
            colnames(cellranger_peaks) <- original_barcodes_adj
            rownames(metadata) <- original_barcodes_adj
            common_barcodes <- intersect(original_barcodes_adj, cellbender_barcodes)
            cat("Added -1 suffix to original barcodes\n")
          }
        }
        cat("Common barcodes found:", length(common_barcodes), "\n")
        if (length(common_barcodes) > 0) {
          # Subset all data to common barcodes
          rna_counts <- cellbender_counts[, common_barcodes]
          cellranger_peaks <- cellranger_peaks[, common_barcodes]
          metadata <- metadata[common_barcodes, ]
          
          cat("Successfully replaced RNA counts with CellBender data\n")
          cat("Final RNA features:", nrow(rna_counts), "cells:", ncol(rna_counts), "\n")
          cat("Final ATAC features:", nrow(cellranger_peaks), "cells:", ncol(cellranger_peaks), "\n")
        } else {
          warning("No common barcodes found between multiome and CellBender data, keeping original RNA counts")
        }
      } else {
        warning(paste("Failed to load CellBender data: keeping original RNA counts"))
      }
    } else {
      cat("CellBender file not found:", cellbender_path, "- keeping original RNA counts\n")
    }
  } else {
    cat("No CellBender path specified - keeping original RNA counts\n")
  }
  # create Seurat object
  x <- CreateSeuratObject(
    counts = rna_counts,  # Now potentially CellBender-processed
    meta.data = metadata
  )
  str_section_head("Base Seurat Object")
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
  # add in the ATAC-seq data
  grange_counts <- StringToGRanges(
    rownames(cellranger_peaks),
    sep = c(":", "-")
  )
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  cellranger_peaks <- cellranger_peaks[as.vector(grange_use), ]
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  x[["ATAC"]] <- CreateChromatinAssay(
    counts = cellranger_peaks,
    sep = c(":", "-"),
    genome = "hg38",
    fragments = paste0(
      samples$atac_path[i],
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
      "percent.mt",
      "percent.rb",
      "percent.hb"
    ),
    ncol = 5,
    log = TRUE
  ) + NoLegend()
  save_figure(
    p1,
    paste0(name, "_counts_vln")
  )
  # more accurate peak calling using MACS2
  peaks <- CallPeaks(x)
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
      samples$atac_path[i],
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
    x$TSS.enrichment > params$tss.score,
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
    subset = nCount_ATAC < params$max.count.atac &
      nCount_ATAC > params$min.count.atac &
      nCount_RNA < params$max.count.rna &
      nCount_RNA > params$min.count.rna &
      pct_reads_in_peaks > params$pct.reads.peaks &
      nucleosome_signal < params$nucleosome &
      TSS.enrichment > params$tss.score &
      percent.mt < params$max.percent.mt
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
  # ATAC analysis
  # We exclude 1st dim as this is typically correlated with sequencing depth
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
  str_section_head("Clustered")
  # Get a list of motif position frequency matrices from the JASPAR database
  pwm_set <- getMatrixSet(
    x = JASPAR2020,
    opts = list(
      species = 9606,
      all_versions = FALSE
    )
  )
  motif_matrix <- CreateMotifMatrix(
    features = granges(x),
    pwm = pwm_set,
    genome = "hg38",
    use.counts = FALSE
  )
  motif_object <- CreateMotifObject(
    data = motif_matrix,
    pwm = pwm_set
  )
  x <- SetAssayData(
    x,
    assay = "peaks",
    slot = "motifs",
    new.data = motif_object
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
  x@meta.data$inflammation <- samples$diagnosis[i]
  assign( # giving names to objects
    name,
    x
  )
}

if (length(samples$sample_id) == 1) {
  save_object(x, "individual")
} else { # integrate
  # create and save list of seurat objects
  objects <- list()
  for (i in samples$sample_id) {
    objects <- c(
      objects,
      get(i) # need get() to call object instead of string
    )
  }
  save_object(objects, "individual")
}

# integration
# https://satijalab.org/signac/1.2.0/articles/integration.html
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
# Save the processed Seurat object
output_file <- paste0(out_path, "pre-harmony_normalized.rds")
saveRDS(x, output_file)
str_section_noloop("Pre-harmony Normalized")

# try to fix error: element-wise pow(): incompatible matrix dimensions 
# none of below worked
# installing rcpp 1.0.8
#setting the variable to batch correct as factor
#reinstalling harmony, devtools::install_github("immunogenomics/harmony") or devtools::install_github("eddelbuettel/harmony",force = TRUE)
# manually upgrade harmony from 0.1 to 1.2.3
x@meta.data$orig.ident <- as.character(x@meta.data$object)
x@meta.data$object <- as.factor(x@meta.data$object)

# Temporarily disable future processing for integration
plan("sequential")

# Method 1A: Standard Seurat Integration (CCA-based)
DefaultAssay(x) <- "RNA"

# Split by batch/sample for integration
x.list <- SplitObject(x, split.by = "object")

# Find integration anchors
x.anchors <- FindIntegrationAnchors(
  object.list = x.list,
  dims = 1:d
)

# Integrate data
x.integrated <- IntegrateData(
  anchorset = x.anchors,
  dims = 1:d
)

# Re-enable multicore processing after integration
plan("multicore", workers = 24)

# Set integrated assay as default and run analysis
DefaultAssay(x.integrated) <- "integrated"
x.integrated <- ScaleData(x.integrated, verbose = FALSE)
x.integrated <- RunPCA(x.integrated, npcs = d, verbose = FALSE)
x.integrated <- RunUMAP(x.integrated, reduction = "pca", dims = 1:d, 
                       reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

# Process ATAC separately
DefaultAssay(x.integrated) <- "peaks"

x.integrated <- RunTFIDF(x.integrated) %>%
  FindTopFeatures(min.cutoff = "q0") %>%
  RunSVD() %>%
  RunUMAP(reduction = "lsi", dims = 2:d, 
          reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Find multimodal neighbors and run WNN
x.integrated <- FindMultiModalNeighbors(
  x.integrated,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:d, 2:d),
  modality.weight.name = "RNA.weight"
) %>%
  RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap", assay = "SCT",
          reduction.key = "wnnUMAP_") %>%
  FindClusters(graph.name = "wsnn", algorithm = 3, resolution = params$res)

# Replace x with integrated object
x <- x.integrated

if (!inherits(x[["peaks"]], "ChromatinAssay")) {
  cat("Converting peaks assay to ChromatinAssay to store motifs\n")
  
  # Get current data
  peak_counts <- GetAssayData(x, assay = "peaks", slot = "counts")
  
  # Create granges
  peak_granges <- StringToGRanges(rownames(peak_counts), sep = c(":", "-"))
  
  # Get annotation
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  genome(annotation) <- "hg38"
  
  # Recreate as ChromatinAssay
  x[["peaks"]] <- CreateChromatinAssay(
    counts = peak_counts,
    sep = c(":", "-"),
    genome = "hg38",
    annotation = annotation,
    motifs = motif_object
  )
}

# Plot by default clusters
p1 <- DimPlot(
  x,
  reduction = "wnn.umap",
  label = TRUE,
  label.size = 4,
  pt.size = 0.5
) + ggtitle("Colored by Clusters") + NoLegend()

# Plot by object (sample)
p2 <- DimPlot(
  x,
  reduction = "wnn.umap",
  group.by = "object",
  pt.size = 0.5
) + ggtitle("Colored by Sample")

# Combine plots side by side
combined_plot <- p1 | p2

save_figure(
  combined_plot,
  "_integrated_dimplot"
)

str_section_head("Integrated Clustered")

# Get a list of motif position frequency matrices from the JASPAR database
pwm_set <- getMatrixSet(
  x = JASPAR2020,
  opts = list(
    species = 9606,
    all_versions = FALSE
  )
)

motif_matrix <- CreateMotifMatrix(
  features = granges(x),
  pwm = pwm_set,
  genome = "hg38",
  use.counts = FALSE
)

motif_object <- CreateMotifObject(
  data = motif_matrix,
  pwm = pwm_set
)

x <- SetAssayData(
  x,
  assay = "peaks",
  slot = "motifs",
  new.data = motif_object
)

# Note that this step can take 30-60 minutes
x <- RunChromVAR(
  object = x,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
str_section_head("Integrated Motif Scored")

warnings()

DefaultAssay(x) <- "SCT"
x <- PrepSCTFindMarkers(x)

markers <- FindAllMarkers(
  x,
  assay = "SCT",
  verbose = FALSE
)

save_object(x, "integrated")
write.csv(markers, paste0(out_path, "all_markers.csv"))

sessionInfo()

print("End of atac-multi.wnn.R")

