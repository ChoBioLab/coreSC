# !/usr/bin/env Rscript

library(tidyverse)

# capture bash vars
args <- commandArgs(trailingOnly = T)

out_path <- paste0(args[1], "/")
params <- read.csv("./config/params.csv", row.names = 1)
samples <- read.csv("./config/samples.csv")
# clusters <- read_csv("./config/clust_ann.csv") %>%
#    remove_rownames %>%
#    column_to_rownames(var = "clust_num")

# package install check and load
packages <- c(
  "Seurat",
  "tidyverse",
  "cowplot",
  "patchwork",
  "future"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")

# convenience functions
save_figure <- function(plots, name, type = "pdf", width, height) {
  if (type == "png") {
    png(paste0(out_path, name, ".", type),
      width = 12, height = 6, units = "in", res = NA
    )
  } else {
    pdf(paste0(out_path, name, ".", type),
      width = 12, height = 6
    )
  }
  print(plots)
  dev.off()
}

save_object <- function(object, name) {
  saveRDS(object, paste0(out_path, name, ".RDS"))
}

save_h5 <- function(object, name) {
  SaveH5Seurat(object, paste0(out_path, name))
}

read_object <- function(name) {
  readRDS(paste0(out_path, name, ".RDS"))
}

str_section_head <- function(title) {
  message("")
  message(strrep("=", 80))
  message(strrep("=", 80))
  message("")
  message(paste(name, title))
  message("")
  message(strrep("=", 80))
  message(strrep("=", 80))
  message("")
  str(x)
}

str_section_noloop <- function(title) {
  message("")
  message(strrep("=", 80))
  message(strrep("=", 80))
  message("")
  message(paste(title))
  message("")
  message(strrep("=", 80))
  message(strrep("=", 80))
  message("")
  str(x)
}

d <- params["dims", ]

save.image(paste0(out_path, "tmp/preamble_image.RData"))

print("End of preamble.R")
