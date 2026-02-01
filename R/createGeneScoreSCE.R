#' Create a SingleCellExperiment for Gene Scores from fragment files
#'
#' Create a SingleCellExperiment from fragment files, storing the gene score matrix.
#'
#' @param fragment.files Vector of strings specifying fragment files. Vector
#'   names need to be sample names.
#' @param gene.grs Genomic Ranges specifying gene coordinates for creating the
#'   gene score matrix.
#' @param output.dir String containing the directory where files should be
#'   output.
#' @param barcodes.list A List with samples as names and the values a vector of
#'   barcodes for that sample.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating how matrix
#'   creation should be parallelized.
#'
#' @return A \linkS4class{SingleCellExperiment} containing the gene score matrix.
#'
#' @author Natalie Fox
#' @importFrom BiocParallel bpparam
#' @export
createGeneScoreSCE <- function(fragment.files,
                               gene.grs,
                               output.dir = tempdir(),
                               barcodes.list = NULL,
                               BPPARAM = bpparam()) {
  if (is.null(names(fragment.files))) {
    stop('Please add sample names as the names on fragment.files')
  }
  
  createSCEFromFragments(
    fragment.files = fragment.files,
    output.dir = output.dir,
    matrix.name = 'GeneAccessibilityMatrix',
    worker.fun = .saveRegionMatrixCall,
    BPPARAM = BPPARAM,
    regions = gene.grs,
    barcodes.list = barcodes.list
  )
}
