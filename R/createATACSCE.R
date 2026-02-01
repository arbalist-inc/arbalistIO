#' Create a SingleCellExperiment from fragment files
#'
#' Create a SingleCellExperiment from fragment files, storing the tile matrix.
#'
#' @param fragment.files Vector of strings specifying fragment files. Vector
#'   names need to be sample names.
#' @param output.dir String containing the directory where files should be
#'   output.
#' @param tile.size Integer scalar specifying the size of the tiles in base
#'   pairs.
#' @param seq.lengths Named integer vector containing the lengths of the
#'   reference sequences used for alignment.
#' @param barcodes.list A List with samples as names and the values a vector of
#'   barcodes for that sample.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix
#'   creation should be parallelized.
#'
#' @return A \linkS4class{SingleCellExperiment} containing the tile matrix.
#'
#' @author Natalie Fox
#' @importFrom BiocParallel bpparam
#' @export
createATACSCE <- function(fragment.files,
                          output.dir = tempdir(),
                          tile.size = 500,
                          seq.lengths = NULL,
                          barcodes.list = NULL,
                          BPPARAM = bpparam()) {
  if (is.null(names(fragment.files))) {
    stop('Please add sample names as the names on fragment.files')
  }
  
  # check that fragment headers contained required information
  if (is.null(seq.lengths)) {
    info <- .processFragmentHeader(fragment.files[1])
    if (!'reference_path' %in% names(info)) {
      stop(
        'the fragment file header does not have reference_path information so please specify the seq.lengths argument'
      )
    }
  }
  
  createSCEFromFragments(
    fragment.files = fragment.files,
    output.dir = output.dir,
    matrix.name = paste0('TileMatrix', tile.size),
    worker.fun = .saveTileMatrixCall,
    BPPARAM = BPPARAM,
    tile.size = tile.size,
    seq.lengths = seq.lengths,
    barcodes.list = barcodes.list
  )
}
