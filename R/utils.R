.processFragmentHeader <- function(file) {
  handle <- gzfile(file, open = "rb")
  on.exit(close(handle))
  all.headers <- character(0)
  
  chunk <- 100
  repeat {
    lines <- readLines(handle, n = chunk)
    header <- startsWith(lines, "#")
    all.headers <- c(all.headers, sub("^# ", "", lines[header]))
    if (length(lines) < chunk || !all(header)) {
      break
    }
  }
  
  field <- sub("=.*", "", all.headers)
  value <- sub("[^=]+=", "", all.headers)
  split(value, field)
}

#' Retrieves genome reference file name
#'
#' Helper function to extract genome reference file name from fragment file header.
#'
#' @param fragment.file String specifying the file name of the fragment file
#' @return Character vector containing the file path to the genome reference
#'   file.
#' @author Natalie Fox
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


#' Filter duplicate features
#'
#' Keeps one feature from each set of duplicated features based on specified
#' criteria, such as the feature with the highest values. Other duplicate
#' features are removed.
#'
#' @param se A \link[SummarizedExperiment]{SummarizedExperiment} object containing
#'    RNA data
#' @param mcol.name String specifying the column name in the rowData 
#'     corresponding to feature names
#' @param assay.name String specifying the name of the assay
#' @param summary.stat Function to summarize each feature (row) of the
#'   experiment. Default is \code{sum}
#' @param selection.metric Function to select the row to keep when there are
#'   duplicate rows with the same \code{mcol.name}. Default is \code{max}
#'
#' @importFrom SummarizedExperiment mcols assay
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} with duplicate features resolved.
#' @author Natalie Fox, Jayaram Kancherla
#' @examples
#' library(SummarizedExperiment)
#'
#' # Create a mock SummarizedExperiment with duplicates
#' counts <- matrix(1:10, ncol=2)
#' rdata <- DataFrame(name = c("GeneA", "GeneA", "GeneB", "GeneC", "GeneC"))
#' se <- SummarizedExperiment(assays=list(counts=counts), rowData=rdata)
#'
#' # Filter duplicates, keeping the row with the max sum
#' se_filtered <- filterDuplicateFeatures(se, mcol.name="name", selection.metric=max)
#'
#' @export
filterDuplicateFeatures <- function(se,
                                    mcol.name = "name",
                                    assay.name = 1,
                                    summary.stat = sum,
                                    selection.metric = max) {
  ids <- mcols(se)[[mcol.name]]
  dup_flag <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
  
  # if no duplicates, return
  if (!any(dup_flag)) {
    return(se)
  } else {
    
    # precompute summary statistic for all rows
    row_summary <- apply(assay(se, assay.name), 1, summary.stat)
    
    dup_indices <- which(dup_flag)
    dup_ids <- ids[dup_flag]
    
    # Split duplicate indices by ID, pick one row per group
    idx_list <- split(dup_indices, dup_ids)
    
    selected_dup_rows <- vapply(
      idx_list,
      FUN.VALUE = integer(1L),
      FUN = function(ix) {
        vals <- row_summary[ix]
        ix[which(vals == selection.metric(vals))[1]]
      }
    )
    
    non_dup_rows <- which(!dup_flag)
    keep <- sort(c(non_dup_rows, selected_dup_rows))
    se[keep, ]
  }
}
