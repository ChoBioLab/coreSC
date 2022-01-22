
# !/usr/bin/env Rscript

library(Seurat)
library(dplyr)
# library(future) # parallelization
library(parallel)

# plan(multicore) # parallelization
# options(future.globals.maxSize = 2000 * 1024^2)

args <- commandArgs(trailingOnly = T)
out_path <- paste0(args[1], "/")

load(paste0(out_path, "tmp/preamble_image.RData"))

numCores <- detectCores()

output <- mclapply(
  mc.cores = 4,
  X = files,
  FUN = process_file
)
