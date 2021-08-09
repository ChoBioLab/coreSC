# !/usr/bin/env Rscript

# capture bash vars
args <- commandArgs(trailingOnly = T)

out_path <- paste0(args[1], "/")
params <- read.csv("./config/params.csv", row.names = 1)
samples <- read.csv("./config/samples.csv")

# package install check and load
packages <- c(
  "Seurat",
  "tidyverse",
  "cowplot",
  "patchwork"
)

package_check <- lapply(
  packages,
  function(x) {
    if (!require(x, character.only = T)) {
      install.packages(x, dependencies = T)
    }
  }
)

# convenience functions
save_figure <- function(plots, name, type = "png", width, height, res) {
  if (type == "png") {
    png(paste0(out_path, name, ".", type),
      width = width, height = height, units = "in", res = 200
    )
  } else {
    pdf(paste0(out_path, name, ".", type),
      width = width, height = height
    )
  }
  print(plots)
  dev.off()
}

save_object <- function(object, name) {
  saveRDS(object, paste0(out_path, name, ".RDS"))
}

read_object <- function(name) {
  readRDS(paste0(out_path, name, ".RDS"))
}

str_section_head <- function(title) {
  message("")
  message(strrep("=", 80))
  message(strrep("=", 80))
  message("")
  message(paste(samples$name[i], title))
  message("")
  message(strrep("=", 80))
  message(strrep("=", 80))
  message("")
  str(x)
}

save.image("./tmp/base_image.RData")

print("End of preamble.R")
