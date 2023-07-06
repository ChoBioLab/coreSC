# coreSC

This pipeline serves as the master process to execute multi-sample preprocessing, normalization, dimensional reduction, integration, and clustering for 10X single cell RNAseq data using the Seurat framework. This documentation assumes a working understanding of Seurat.

## Features

## Dependencies
- coreSC is run with the [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) container execution framework.

## Setup

### First Time
```sh

# clone the repo
git clone https://github.com/ChoBioLab/coreSC.git --recurse-submodules

# put config templates in-place
cd coreSC/config && cp templates/* .
```

### Config
1. Required configuration
    * `env`
    * `params.csv`
    * `samples.csv`
    * `ann`
1. Confirm parallel memory use with future
    * The `future.mem` param gives `RAM / thread`. Each individual task needs an adequate threshold of RAM to complete its work. Also `future.mem * future.workers` gives the total memory allocation. This should live under the available system RAM for the job as a whole. If either of these considerations are not met, the run will fail!

## Usage
- Execution can be carried out with the `run` script.

`./run`
- `-v` *version* [REQUIRED]
- `-h` *harmonize* [NULL]
- `-a` *atac-multi subroutine* [NULL]
- `-c` *cite-seq subroutine* [NULL]

- Executing `run` alone will perform a standard, single-modal regularization and integration if relevant.
- Specifying a Seurat container version to use is required. Versions can be found at https://gallery.ecr.aws/chobiolab

```sh
# examples

./run -v v4-r2           # regularize and integrate with SCTransform method

./run -v v4-r2 -h TRUE   # substitute integration method for harmony

./run -v v4-r2 -a TRUE   # alternatively, pass samples through atac-multiome routine
```

### Output
```sh
coreSC/output/output_2022-12-21_18.59.26
├── all_markers.csv                 # output of FindAllMarkers()
├── combined_dimplot_red.pdf
├── individual_clustered.RDS        # list object of all individual objects (regularized & clustered)
├── integrated.RDS                  # single, integrated object (regularized & clustered)
├── log.txt
├── params.csv                      # input config params.csv
├── sample1_individual_dimplot.pdf
├── sample1_unfilt_scatter.pdf
├── sample1_unfilt_vln.pdf
├── sample1_var_features.pdf
├── sample2_individual_dimplot.pdf
├── sample2_unfilt_scatter.pdf
├── sample2_unfilt_vln.pdf
├── sample2_var_features.pdf
└── samples.csv                     # input config samples.csv

coreSC/output/output_2023-07-06_19.32.36
├── all_markers.csv                     # output of FindAllMarkers() unannotated
├── annotation                          # CELLTYPIST OUTPUT
│   ├── annd_all_markers.csv            # output of FindAllMarkers() annotated
│   ├── celltypist-log.txt
│   ├── decision_matrix.csv             # celltypist output
│   ├── integrated-annd_2023-07-06.RDS  # ANNOTATED OBJECT
│   ├── predicted_labels.csv            # celltypist output
│   ├── probability_matrix.csv          # celltypist output
│   └── qc.csv                          # qc metrics for annotation
├── combined_dimplot_red.pdf
├── individual_clustered.RDS            # list of individual, norm'd, clustered objects
├── integrated.h5seurat                 # integrated hdf5 object unannotated
├── integrated.RDS                      # integrated seurat object unannotated
├── log.txt
├── params.csv                          # copy of process params
├── sample1_individual_dimplot.pdf
├── sample1_unfilt_scatter.pdf
├── sample1_unfilt_vln.pdf
├── sample1_var_features.pdf
└── samples.csv                         # copy of sample sheet

1 directory, 19 files
```

## Reference

### File Tree
```sh
coreSC/
├── annotation
│   ├── README.md
│   ├── run                     # annotation execution script
│   └── src
│       ├── apply-ann.R         # process to apply annotation back to input
│       ├── getopts             # annotation run arguments
│       ├── matrix-convert.R    # seurat to sparse matrix conversion
│       └── qc.R                # qc metrics generation process
├── config                      # COPY CONFIG FILES HERE
│   └── templates
│       ├── ann                 # env vars for annotation
│       ├── env                 # general env vars
│       ├── params.csv          # pipeline parameters
│       └── samples.csv         # sample details
├── LICENSE
├── main                        # orchestration script
├── output
├── README.md
├── run                         # EXECUTION SCRIPT
└── scripts
    ├── atac-multi-wnn.R        # full atac-multiome routine
    ├── cite-wnn.R              # full cite-seq multiome routine
    ├── create-multi-norm.R     # create object routine for multiomic samples
    ├── create-object-norm.R    # create object routine with classical normalization
    ├── create-object-sct.R     # create object routine with SCTransform normalization
    ├── getopts                 # run script arguments
    ├── harmonize.R             # harmony integration subroutine
    ├── integrate.R             # standard seurat integration subroutine
    ├── main.R                  # legacy
    ├── mixscape.R              # method for crispr preprocessing
    └── preamble.R              # place setting subroutine

6 directories, 26 files
```

### params.csv

| Var | Description |
| ----- | ----- |
| min.cells | minimum number of cells |
| min.count.rna | minimum number of rna molecules per cell |
| max.count.rna | max number of rna molecules per cell |
| min.count.atac | minimum number of atac fragments per cell |
| max.count.atac | max number of atac fragments per cell |
| min.features | minimum number of genes per cell |
| max.features | max number of genes per cell |
| max.percent.mt | max threshold for mitochondrial content |
| min.percent.mt | minimum threshold for mitochondrial content |
| pct.reads.peaks | proportion of reads found in peaks |
| nucleosome | nucleosome score |
| tss.score | transcription start site score |
| dims | number of dimensions to include |
| res | clustering resolution |
| future.mem | memory allocation per thread in GB |
| future.workers | number of parallel threads |

### samples.csv

| Var | Description |
| ----- | ----- |
| name | sample name (must be unique, cannot be integer) |
| dir | path to directory containing matrices |
| project | project variable name (applied as object metadata) |
| group | group label (applied as object metadata) |

