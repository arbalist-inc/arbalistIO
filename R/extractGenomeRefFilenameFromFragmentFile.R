#' Retrieves genome reference file name
#'
#' Helper function to extract genome reference file name from fragment file
#' header.
#'
#' @param fragment.file String specifying fragment file name.
#' @return Character vector containing the file path to the genome reference
#'   file.
#' @author Natalie Fox
#'
#' @examples
#' # Mock a fragment file
#' f <- tempfile(fileext=".tsv.gz")
#' mockFragmentFile(f, c(chr1=1000), 10, LETTERS[1:5],
#'                  comments=c("reference_path=/path/to/ref"))
#'
#' extractGenomeRefFilenameFromFragmentFile(f)
#'
#' @export
extractGenomeRefFilenameFromFragmentFile <- function(fragment.file) {
  info <- .processFragmentHeader(fragment.file)
  
  if (!"reference_path" %in% names(info))
    return(NULL)
  
  fai <- file.path(info$reference_path, "fasta", "genome.fa.fai")
  return(fai)
}
