library(Seurat)
library(Signac)
library(SeuratDisk)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(future)

plan(multicore) # parallelization
options(future.globals.maxSize = 2000 * 1024^2)

# load H5 reference for cluster mapping
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

# the 10x hdf5 file contains both data types.
inputdata.10x <- Read10X_h5("/mnt/cho_lab/disk2/ctastad/projects/perianal-cd/analysis/cellranger/cr_arc_count_2021-12-23_1559/AA5_Non/outs/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
x <- CreateSeuratObject(counts = rna_counts)
x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

frag.file <- "/mnt/cho_lab/disk2/ctastad/projects/perianal-cd/analysis/cellranger/cr_arc_count_2021-12-23_1559/AA5_Non/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = frag.file,
  min.cells = 10,
  min.features = 200
  annotation = annotations
)

x[["ATAC"]] <- chrom_assay

VlnPlot(x,
  features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0
) + NoLegend()

# compute nucleosome signal score per cell
x <- NucleosomeSignal(object = x)

# compute TSS enrichment score per cell
x <- TSSEnrichment(object = x, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
x$pct_reads_in_peaks <- x$peak_region_fragments / x$passed_filters * 100
x$blacklist_ratio <- x$blacklist_region_fragments / x$peak_region_fragments
x$high.tss <- ifelse(x$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(x, group.by = 'high.tss') + NoLegend()
x$nucleosome_group <- ifelse(x$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = x, group.by = 'nucleosome_group')

VlnPlot(
  object = x,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

x <- subset(
  x = x,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# RNA analysis
DefaultAssay(x) <- "RNA"
x <- SCTransform(x, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

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
refquery <- merge(reference, x)
refquery[["spca"]] <- merge(reference[["spca"]], x[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = "spca", dims = 1:50)
DimPlot(refquery, group.by = "id", shuffle = TRUE)

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(x) <- "ATAC"
x <- RunTFIDF(x)
x <- FindTopFeatures(x, min.cutoff = "q0")
x <- RunSVD(x)

DepthCor(x)

x <- RunUMAP(x, reduction = "lsi", dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

x <- FindMultiModalNeighbors(x, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
x <- RunUMAP(x, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
x <- FindClusters(x, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# perform sub-clustering on cluster 6 to find additional structure
x <- FindSubCluster(x, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(x) <- "sub.cluster"

p1 <- DimPlot(x, reduction = "umap.rna", group.by = "predicted.celltype.l2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(x, reduction = "umap.atac", group.by = "predicted.celltype.l2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(x, reduction = "wnn.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

gene.activities <- GeneActivity(x)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
x[['RNA']] <- CreateAssayObject(counts = gene.activities)
x <- NormalizeData(
  object = x,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(x$nCount_RNA)
)

## to make the visualization easier, subset T cell clusters
#celltype.names <- levels(x)
celltype.names <- unique(x@meta.data$predicted.celltype.l2)
tcell.names <- grep("CD4|CD8|Treg", celltype.names, value = TRUE)
tcells <- subset(x, idents = tcell.names)
CoveragePlot(tcells, region = "CD8A", features = "CD8A", assay = "ATAC", expression.assay = "SCT", peaks = FALSE)

saveRDS(x, file = "cluster.RDS")

