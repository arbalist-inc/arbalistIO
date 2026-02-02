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
#' @param fragment.file String specifying fragment file name
#' @return Named string vector
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
#' @param se \linkS4class{SummarizedExperiment}
#' @param mcol.name String specifying the colname for the experiment rowData
#' @param summary.stat Function to summarize each feature (row) of the
#'   experiment
#' @param selection.metric Function to select the row to keep when there are
#'   duplicate rows with the same mcol.name
#'
#' @importFrom SummarizedExperiment mcols assay
#' @return A \linkS4class{SummarizedExperiment} with duplicate features resolved.
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
                                    mcol.name = 'name',
                                    summary.stat = sum,
                                    selection.metric = max) {
  duplicate.values <- names(which(table(mcols(se)[, mcol.name]) > 1))
  if (length(duplicate.values) == 0) {
    return(se)
  }
  non.duplicate.rows <- which(!mcols(se)[, mcol.name] %in% duplicate.values)
  duplicate.rows <- which(mcols(se)[, mcol.name] %in% duplicate.values)
  selected.duplicate.rows <- sapply(duplicate.values, function(i) {
    duplicate.rows <- which(mcols(se)[, mcol.name] %in% i)
    row.summary.stats <- apply(assay(se)[duplicate.rows, ], 1, summary.stat)
    return(duplicate.rows[which(row.summary.stats == selection.metric(row.summary.stats))[1]])
  })
  
  se <- se[sort(c(non.duplicate.rows, selected.duplicate.rows)), ]
  
  return(se)
}
