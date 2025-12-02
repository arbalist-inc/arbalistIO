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
