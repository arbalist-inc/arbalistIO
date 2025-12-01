library(methods)
library(Matrix)
library(rhdf5)

# Helper to create a mock Cell Ranger HDF5 file
mockCellRangerH5 <- function(
  filepath,
  n_genes = 100,
  n_cells = 50,
  sample_name = "Sample1"
) {
  if (file.exists(filepath)) {
    file.remove(filepath)
  }

  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("rhdf5 package is required for testing")
  }

  rhdf5::h5createFile(filepath)
  rhdf5::h5createGroup(filepath, "matrix")

  # Create sparse matrix data
  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = 1),
    nrow = n_genes,
    ncol = n_cells
  )
  sparse_counts <- as(counts, "dgCMatrix")

  # Write matrix components in CSR
  rhdf5::h5write(sparse_counts@x, filepath, "matrix/data")
  rhdf5::h5write(sparse_counts@i, filepath, "matrix/indices")
  rhdf5::h5write(sparse_counts@p, filepath, "matrix/indptr")
  rhdf5::h5write(dim(sparse_counts), filepath, "matrix/shape")

  # Write barcodes
  barcodes <- paste0(seq_len(n_cells), "-1")
  rhdf5::h5write(barcodes, filepath, "matrix/barcodes")

  # Write features
  rhdf5::h5createGroup(filepath, "matrix/features")
  feature_ids <- paste0("GENE", seq_len(n_genes))
  feature_names <- paste0("Gene", seq_len(n_genes))
  feature_types <- rep("Gene Expression", n_genes)

  # Add intervals for half, some NA's to test filter.features.without.intervals
  intervals <- c(rep("chr1:100-200", n_genes / 2), rep("NA", n_genes / 2))

  rhdf5::h5write(feature_ids, filepath, "matrix/features/id")
  rhdf5::h5write(feature_names, filepath, "matrix/features/name")
  rhdf5::h5write(feature_types, filepath, "matrix/features/feature_type")
  rhdf5::h5write(intervals, filepath, "matrix/features/interval")
  rhdf5::h5write(rep("Genome1", n_genes), filepath, "matrix/features/genome")
}
