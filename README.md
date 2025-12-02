# arbalistIO

**arbalistIO** provides streamlined functions for importing single-cell ATAC-seq and Multiome (ATAC + RNA) data into an [MultiAssayExperiment](https://bioconductor.org/packages/MultiAssayExperiment). It works efficiently with standard outputs from 10x Cell Ranger outputs, using C++ backends to parse fragment files and generate HDF5-backed sparse matrices.

## Installation

You can install the development version:

```r
BiocManager::install("Jfortin1/arbalistIO")
```

## Usage

`createArbalistMAE` takes paths to your fragment files and (optionally) filtered feature matrices to construct the experiment object.

```r
library(arbalistIO)
library(GenomicRanges)

# Paths to your 10x Cell Ranger outputs
sample_name <- "Sample1"
fragment_file <- "path/to/atac_fragments.tsv.gz"
feature_matrix_h5 <- "path/to/filtered_feature_bc_matrix.h5"

# Define genomic regions of interest (e.g., gene promoters)
gene_regions <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = c(1000, 5000), width = 500),
  name = c("GeneA", "GeneB")
)

mae <- createArbalistMAE(
  sample.names = sample_name,
  fragment.files = fragment_file,
  filtered.feature.matrix.files = feature_matrix_h5,
  gene.grs = gene_regions,
  tile.size = 500,
  seq.lengths = c(chr1 = 10000)
)

print(mae)

# Access specific experiments
tiles <- experiments(mae)[["TileMatrix500"]]
rna <- experiments(mae)[["GeneExpressionMatrix"]]

head(colData(mae))
```

The package also includes functions to process fragment files and generate appropriate tile or region matrices, e.g.

```r
# Mock a fragment file
tmp_frag <- tempfile(fileext = ".fragments.gz")
arbalistIO::mockFragmentFile(
  output.file = tmp_frag,
  seq.lengths = c(chr1 = 10000, chr2 = 5000),
  num.fragments = 1000,
  cell.names = c("CellA", "CellB", "CellC")
)

# Generate a tile matrix from the mock file
out_h5 <- tempfile(fileext = ".h5")
matrix_res <- arbalistIO::saveTileMatrix(
  fragment.file = tmp_frag,
  output.file = out_h5,
  output.name = "tiles",
  tile.size = 500
)

print(matrix_res$counts)
```

