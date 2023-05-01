library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
out_path <- paste0(args[1], "/")
load(paste0(out_path, "tmp/preamble_image.RData"))

plan(
  multicore,
  workers = params["future.workers", ]
) # parallelization
options(future.globals.maxSize = params["future.mem", ] * 1024^2)

for (i in seq_len(nrow(samples))) {
}

print("End of mixscape.R")
