#' Create a mock fragment file
#'
#' Mock up a fragment file for examples and testing.
#'
#' @param output.file String containing a path to an output file.
#' @param seq.lengths Named integer vector containing the lengths of the
#'   reference sequences used for alignment. Vector names should correspond to
#'   the names of the sequences.
#' @param num.fragments Integer scalar, the average number of fragments per
#'   cell.
#' @param cell.names Character vector containing the cell names. The length of
#'   this vector is used as the total number of cells.
#' @param width.range Integer vector of length 2, containing the range of
#'   possible fragment widths in base pairs.
#' @param read.range Integer vector of length 2, containing the range of the
#'   read count supporting each fragment.
#' @param comments Character vector of comments to be added to the start of the
#'   file.
#' @param compressed Logical scalar indicating whether the output file should be
#'   compressed.
#'
#' @return A fragment file is created at \code{output.file}. \code{NULL} is
#'   invisibly returned.
#'
#' @examples
#' temp <- tempfile(fileext=".fragments.gz")
#' mockFragmentFile(temp, c(chrA=1000, chrB=200000, chrC=200),
#'     num.fragments=100, cell.names=LETTERS)
#'
#' X <- read.table(temp)
#' head(X)
#'
#' @author Aaron Lun
#'
#' @importFrom stats runif
#' @importFrom utils write.table
#' @export
mockFragmentFile <- function(output.file,
                             seq.lengths,
                             num.fragments,
                             cell.names,
                             width.range = c(10, 1000),
                             read.range = c(1, 10),
                             comments = NULL,
                             compressed = TRUE) {
  number <- num.fragments * length(cell.names)
  seq <- sample(names(seq.lengths), number, replace = TRUE)
  limits <- (seq.lengths - 1L)[seq] # avoid overlapping the end position.
  starts <- floor(runif(number) * limits) 
  ends <- pmin(starts + floor(runif(number, width.range[1], width.range[2])), limits+1L) # since ends are exclusive, ok to overlap the end position
  
  o <- order(factor(seq, names(seq.lengths)), starts)
  df <- data.frame(
    seq,
    as.integer(starts),
    as.integer(ends),
    name = sample(cell.names, number, replace = TRUE),
    count = as.integer(floor(runif(number, read.range[1], read.range[2])))
  )
  df <- df[o, ]
  
  handle <- if (compressed) {
    gzfile(output.file, open = "wb")
  } else {
    file(output.file, open = "wb")
  }
  on.exit(close(handle))
  
  if (length(comments)) {
    writeLines(con = handle, paste0("# ", comments), sep = "\n")
  }
  write.table(
    file = handle,
    df,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE,
    eol = "\n"
  )
  
  invisible(NULL)
}


#' Create a mock Cell Ranger h5 file
#'
#' Mock up a Cell Ranger h5 file for examples and testing.
#'
#' @param filepath String containing a path to the output file.
#' @param n.genes Integer indicating the number of genes 
#' @param n.cells Integer indicating the number of cells
#' @param cell.names Character vector containing the cell names. The length of
#' this vector is used as the total number of cells.
#' 
#'
#' @return A h5 file mimicking a typical Cell Ranger output is created at \code{filepath}. 
#'
#' @examples
#' temp <- tempfile(fileext = ".h5")
#' mockCellRangerH5(temp, cell.names = LETTERS)
#' rhdf5::h5ls(temp)
#'
#' @author Jayaram Kancherla
#'
#' @importFrom rhdf5 h5createFile h5createGroup h5write
#' @export
mockCellRangerH5 <- function(
    filepath,
    n.genes = 100,
    n.cells = NULL,
    cell.names = NULL
) {
  
  if (file.exists(filepath)) {
    file.remove(filepath)
  }
  
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("rhdf5 package is required for testing")
  }
  
  if(is.null(cell.names) && is.null(n.cells)) {
    stop("either n.cells or cell.names must be provided")
  } else if(!is.null(cell.names) && is.null(n.cells)) {
    n.cells <- length((cell.names))
  } else if(is.null(cell.names) && !is.null(n.cells)) {
    cell.names <- paste0("Cell", seq_len(n.cells))
  } else if(n.cells != length(cell.names)) {
    stop("n.cells and cell.names do not match")
  }
  
  h5createFile(filepath)
  h5createGroup(filepath, "matrix")
  
  # Create sparse matrix data
  counts <- matrix(
    stats::rpois(n.genes * n.cells, lambda = 1),
    nrow = n.genes,
    ncol = n.cells
  )
  sparse_counts <- as(counts, "dgCMatrix")
  
  # Write matrix components in CSR
  h5write(sparse_counts@x, filepath, "matrix/data")
  h5write(sparse_counts@i, filepath, "matrix/indices")
  h5write(sparse_counts@p, filepath, "matrix/indptr")
  h5write(dim(sparse_counts), filepath, "matrix/shape")
  
  # Write barcodes
  barcodes <- cell.names
  h5write(barcodes, filepath, "matrix/barcodes")
  
  # Write features
  h5createGroup(filepath, "matrix/features")
  feature_ids <- paste0("GENE", seq_len(n.genes))
  feature_names <- paste0("Gene", seq_len(n.genes))
  feature_types <- rep("Gene Expression", n.genes)
  
  # Add intervals for half, some NA's to test filter.features.without.intervals
  n_with_intervals <- ceiling(n.genes / 2)
  n_na <- n.genes - n_with_intervals
  intervals <- c(rep("chr1:100-200", n_with_intervals), rep("NA", n_na))
  
  h5write(feature_ids, filepath, "matrix/features/id")
  h5write(feature_names, filepath, "matrix/features/name")
  h5write(feature_types, filepath, "matrix/features/feature_type")
  h5write(intervals, filepath, "matrix/features/interval")
  h5write(rep("Genome1", n.genes), filepath, "matrix/features/genome")
}
