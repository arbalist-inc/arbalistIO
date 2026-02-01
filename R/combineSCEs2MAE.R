#' Combine SingleCellExperiments into a MultiAssayExperiment
#'
#' Combine ATAC and RNA SingleCellExperiments into a MultiAssayExperiment.
#' 
#' @param atac.sce A \linkS4class{SingleCellExperiment} containing ATAC data.
#' @param rna.sce A \linkS4class{SingleCellExperiment} containing RNA data.
#' @param sample.annotation A data.frame containing sample annotation.
#' @param use.alt.exp Logical for selecting the MultiAssayExperiment structure.
#' @param main.exp.name String containing the name of the experiment that will
#'   be the main experiment when use.alt.exp is TRUE.
#'
#' @return A \linkS4class{MultiAssayExperiment}
#'
#' @author Natalie Fox
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList listToMap
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment altExp<- mainExpName
#' @export
combineSCEs2MAE <- function(atac.sce = NULL,
                            rna.sce = NULL,
                            sample.annotation = NULL,
                            use.alt.exp = FALSE,
                            main.exp.name = 'TileMatrix500') {
  
  if (is.null(atac.sce) && is.null(rna.sce)) {
    stop("At least one of atac.sce or rna.sce must be provided.")
  }

  exp.list <- list()
  if (!is.null(atac.sce)) {
      # Use the mainExpName if available, otherwise default to 'atac'
      name <- mainExpName(atac.sce)
      if (is.null(name)) name <- "atac"
      exp.list[[name]] <- atac.sce
  }
  if (!is.null(rna.sce)) {
      name <- mainExpName(rna.sce)
      if (is.null(name)) name <- "GeneExpressionMatrix"
      exp.list[[name]] <- rna.sce
  }

  # If asked for, use alternative experiments for the MAE structure
  if (use.alt.exp) {
    if (!main.exp.name %in% names(exp.list) & length(exp.list) > 1) {
      cell.union <- colnames(exp.list[[1]])
      cell.intersect <- colnames(exp.list[[1]])
      for (i in 2:length(exp.list)) {
        cell.union <- union(cell.union, colnames(exp.list[[i]]))
        cell.intersect <- intersect(cell.intersect, colnames(exp.list[[i]]))
      }
      if (length(cell.union) == length(cell.intersect)) {
        # Convert the list of SCEs to main/alternative experiments in one SCE
        main.sce <- exp.list[[main.exp.name]]
        for (matrix.name in setdiff(names(exp.list), main.exp.name)) {
          altExp(main.sce, matrix.name) <- exp.list[[matrix.name]]
        }
        mainExpName(main.sce) <- main.exp.name
        if ('GeneExpressionMatrix' %in% names(exp.list)) {
          exp.list <- list(multiome = main.sce)
        } else {
          exp.list <- list(atac = main.sce)
        }
      }
    }
  }
  
  # Create the sample map for the MultiAssayExperiment
  el <- ExperimentList(exp.list)
  maplist <- lapply(exp.list, function(se) {
    data.frame(
      primary = se$Sample,
      colname = colnames(se),
      stringsAsFactors = FALSE
    )
  })
  sampMap <- listToMap(maplist)
  
  # Create and annotate the MultiAssayExperiment
  mae <- MultiAssayExperiment(el,
                              sampleMap = sampMap,
                              colData = DataFrame(row.names = unique(sampMap$primary)))

  if (!is.null(sample.annotation)) {
       colData(mae) <- cbind(colData(mae), sample.annotation[rownames(colData(mae)), , drop=FALSE])
  }
  
  return(mae)
}
