#' Filter duplicate features
#'
#' Keeps one feature from each set of duplicated features based on specified
#' criteria, such as the feature with the highest values. Other duplicate
#' features are removed.
#'
#' @param se \linkS4class{SummarizedExperiment}.
#' @param mcol.name String specifying the column name containing
#'   feature id's in the experiment's rowData.
#' @param summary.stat Function to summarize each feature (row) of the
#'   experiment.
#' @param selection.metric Function to select the row to keep when there are
#'   duplicate rows with the same \code{mcol.name}.
#'
#' @importFrom SummarizedExperiment mcols assay
#' @return A \linkS4class{SummarizedExperiment} with duplicate features
#'   resolved.
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
                                    summary.stat = sum,
                                    selection.metric = max) {
  ids <- mcols(se)[[mcol.name]]
  dup_flag <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
  
  # if no duplicates, return
  if (!any(dup_flag)) {
    return(se)
  }
  
  # precompute summary statistic for all rows
  row_summary <- apply(assay(se), 1, summary.stat)
  
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
