# !/usr/bin/env Rscript

install.packages("devtools")
#install.packages("dplyr")
#install.packages("cowplot")
#install.packages("patchwork")

require(devtools)

install_version("Seurat", version = "4.0.2")

