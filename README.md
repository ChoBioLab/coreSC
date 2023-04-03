# coreSC

This pipeline serves as the master process to execute multi-sample preprocessing, normalization, dimensional reduction, integration, and clustering for 10X single cell RNAseq data.

## Features

## Dependencies
- coreSC is run with the [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) container execution framework.

## Setup
```sh
git clone https://github.com/ChoBioLab/coreSC.git
```

## Usage

total mem use is future.mem * future.workers

### Template Reference

##### params.csv

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
| future.mem | memory allocation per thread in MB |
| future.workers | number of parallel threads |

##### samples.csv

| Var | Description |
| ----- | ----- |
| name | sample name (must be unique, cannot be integer) |
| dir | path to directory containing matrices |
| project | project variable name |
| group | grouping metadata parameter |

## Output
