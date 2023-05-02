library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(future)

args <- commandArgs(trailingOnly = TRUE)
out_path <- paste0(args[1], "/")
load(paste0(out_path, "tmp/preamble_image.RData"))
obj <- read_object("individual_clustered")

plan(
  multicore,
  workers = params["future.workers", ]
) # parallelization
options(future.globals.maxSize = params["future.mem", ] * 1024^2)

name <- unique(obj)

# Calculate perturbation signature (PRTB).
x <- CalcPerturbSig(
  object = obj,
  assay = "RNA",
  slot = "data",
  gd.class = "gene",
  nt.cell.class = "NT",
  reduction = "pca",
  ndims = 40,
  num.neighbors = 20,
  split.by = "replicate",
  new.assay.name = "PRTB"
)

# Prepare PRTB assay for dimensionality reduction:
# Normalize data, find variable features and center data.
DefaultAssay(object = x) <- "PRTB"

# Use variable features from RNA assay.
VariableFeatures(object = x) <- VariableFeatures(object = x[["RNA"]])
x <- ScaleData(
  object = x,
  do.scale = FALSE,
  do.center = TRUE
)

# Run PCA to reduce the dimensionality of the data.
x <- RunPCA(
  object = x,
  reduction.key = "prtbpca",
  reduction.name = "prtbpca"
)

# Run UMAP to visualize clustering in 2-D.
x <- RunUMAP(
  object = x,
  dims = 1:40,
  reduction = "prtbpca",
  reduction.key = "prtbumap",
  reduction.name = "prtbumap"
)

# Generate plots to check if clustering is driven by biological replicate ID,
# cell cycle phase or target gene class.
p1 <- DimPlot(
  object = x,
  group.by = "replicate",
  reduction = "prtbumap",
  pt.size = 0.2,
  cols = "Dark2",
  label = FALSE,
  repel = TRUE
) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Biological Replicate") +
  ylab("UMAP 2") +
  xlab("UMAP 1")

p2 <- DimPlot(
  object = x,
  group.by = "Phase",
  reduction = "prtbumap",
  pt.size = 0.2, label = FALSE, repel = TRUE
) +
  ggtitle("Cell Cycle Phase") +
  ylab("UMAP 2") +
  xlab("UMAP 1")

p3 <- DimPlot(
  object = x,
  group.by = "crispr",
  reduction = "prtbumap",
  split.by = "crispr",
  ncol = 1,
  pt.size = 0.2,
  cols = c("grey39", "goldenrod3")
) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1")

# Visualize plots.
save_figure(
  (p1 / p2 + plot_layout(guides = "auto") | p3),
  paste0(name, "_pertub_stat_dimplot"),
  width = 12,
  height = 6
)

# Run mixscape.
x <- RunMixscape(
  object = x,
  assay = "PRTB",
  slot = "scale.data",
  labels = "gene",
  nt.class.name = "NT",
  min.de.genes = 5,
  iter.num = 10,
  de.assay = "RNA",
  verbose = FALSE,
  prtb.type = "KO"
)

# Calculate percentage of KO cells for all target gene classes.
y <- prop.table(table(x$mixscape_class.global, x$NT), 2)
y <- reshape2::melt(y)
y$Var2 <- as.character(y$Var2)
test <- y[which(y$Var1 == "KO"), ]
test <- test[order(test$value, decreasing = TRUE), ]
new_levels <- test$Var2

y$Var2 <- factor(
  y$Var2,
  levels = new_levels
)

y$Var1 <- factor(
  y$Var1,
  levels = c("NT", "NP", "KO")
)

y$gene <- sapply(
  as.character(y$Var2),
  function(x) strsplit(x, split = "g")[[1]][1]
)

y$guide_number <- sapply(
  as.character(y$Var2),
  function(x) strsplit(x, split = "g")[[1]][2]
)

y <- y[-c(which(y$gene == "NT")), ]

p1 <- ggplot(y, aes(x = guide_number, y = value * 100, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
  ylab("% of cells") +
  xlab("sgRNA") +
  theme(
    axis.text.x = element_text(size = 18, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  facet_wrap(vars(gene), ncol = 5, scales = "free") +
  labs(fill = "mixscape class") +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

save_figure(
  p1,
  paste0(name, "_by_class"),
  width = 12,
  height = 6
)

# Remove non-perturbed cells and run LDA dimensionality reduction
Idents(x) <- "mixscape_class.global"
sub <- subset(x, idents = c("KO", "NT"))

# Run LDA.
sub <- MixscapeLDA(
  object = sub,
  assay = "RNA",
  pc.assay = "PRTB",
  labels = "gene",
  nt.label = "NT",
  npcs = 10,
  logfc.threshold = 0.25,
  verbose = FALSE
)

# Use LDA results to run UMAP and visualize cells on 2-D.
# note the number of the dimensions to be used is equal to the number of
# labels minus one (to account for NT cells).
sub <- RunUMAP(
  object = sub,
  dims = 1:11,
  reduction = "lda",
  reduction.key = "ldaumap",
  reduction.name = "ldaumap"
)

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)

# Set colors for each perturbation.
col <- setNames(object = hue_pal()(12), nm = levels(sub$mixscape_class))
names(col) <- c(names(col)[1:7], "NT", names(col)[9:12])
col[8] <- "grey39"

p1 <- DimPlot(
  object = sub,
  reduction = "ldaumap",
  repel = TRUE,
  label.size = 5,
  label = TRUE,
  cols = col
) + NoLegend()

p2 <- p1 +
  scale_color_manual(values = col, drop = FALSE) +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme

save_figure(
  p2,
  paste0(name, "_LDA_clusters"),
  width = 12,
  height = 6
)

print("End of mixscape.R")
