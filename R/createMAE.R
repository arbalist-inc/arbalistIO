#' Combine SingleCellExperiments into a MultiAssayExperiment
#'
#' Combine ATAC and RNA SingleCellExperiments into a MultiAssayExperiment.
#'
#' @param atac.sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'    containing ATAC data
#' @param rna.sce A \link[SingleCellExperiment]{SingleCellExperiment} object
#'     containing RNA data
#' @param sample.annotation A data.frame containing sample annotation.
#' @param rna.name String indicating the name of \code{rna.sce}
#' @param atac.name String indicating the name of \code{atac.sce}
#'
#' @return A \link[MultiAssayExperiment]{MultiAssayExperiment}.
#' @examples
#' # create rna sce
#' temp_rna <- tempfile(fileext = ".h5")
#' mockCellRangerH5(temp_rna, cell.names = LETTERS)
#' rna_sce <- createRNASCE(temp_rna, sample.names="sample1")
#'
#' # create atac sce
#' temp_atac <- tempfile(fileext=".fragments.gz")
#' mockFragmentFile(temp_atac, c(chrA=1000, chrB=200000, chrC=200),
#'                  num.fragments=100, cell.names=LETTERS)
#' tile_sce <- createTileSCE(fragment.files=temp_atac,
#'                           sample.names="sample1",
#'                           tile.size=500,
#'                           seq.lengths = c(chrA=1000, chrB=200000, chrC=200))
#' # create combined MAE
#' mae <- createArbalistMAE(tile_sce, rna_sce)
#' file.remove(list.files(tempdir(), pattern=".h5", full.names = TRUE))
#' @author Natalie Fox, Jayaram Kancherla
#' @importFrom MultiAssayExperiment MultiAssayExperiment ExperimentList listToMap
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment altExp<- mainExpName
#' @export
createArbalistMAE <- function(atac.sce,
                              rna.sce,
                              sample.annotation = NULL,
                              rna.name = "GeneExpressionMatrix",
                              atac.name = "TileMatrix500") {
  exp.list <- list()
  exp.list[[rna.name]] <- rna.sce
  exp.list[[atac.name]] <- atac.sce
  
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
    colData(mae) <- cbind(colData(mae), sample.annotation[rownames(colData(mae)), , drop =
                                                            FALSE])
  }
  
  return(mae)
}


#' Import a single cell ATAC-seq MultiAssayExperiment
#'
#' Import results from scATAC-seq or multiome Cell Ranger results directories
#' into a MultiAssayExperiment.
#'
#' @param cellranger.dirs Named character vector specifying a Cell Ranger
#'   scATAC-seq or multiome results directory. Vector names need to be sample
#'   names.
#' @param fragment.file Character string indicating the name of the fragment file
#' @param h5.file Character string indicating the name of feature matrix h5 file
#' @param barcode.file Character string indicating the name of the barcode file
#' @param summary.file Character string indicating the name of the summary file
#' @inheritParams createArbalistMAE
#' @inheritParams createRegionSCE
#' @inheritParams createTileSCE
#' @inheritParams createRNASCE
#'
#' @return A \link[MultiAssayExperiment]{MultiAssayExperiment}.
#'
#' @author Natalie Fox, Xiaosai Yao
#' @export
#'
createMAEFromCellranger <- function(cellranger.dirs,
                                    sample.names,
                                    output.dir = tempdir(),
                                    tile.size = 500,
                                    seq.lengths = NULL,
                                    regions = NULL,
                                    atac.name = "TileMatrix500",
                                    rna.name = "GeneExpressionMatrix",
                                    fragment.file = "atac_fragments.tsv.gz",
                                    h5.file = "filtered_feature_bc_matrix.h5",
                                    barcode.file = "per_barcode_metrics.csv",
                                    summary.file = "summary.csv",
                                    filter.features.without.intervals = TRUE,
                                    BPPARAM = bpparam()) {
  sample.names <- names(cellranger.dirs)
  
  fragment.files <- .getFilesFromResDirs(cellranger.dirs, fragment.file)
  
  filtered.feature.matrix.files <- .getFilesFromResDirs(cellranger.dirs, h5.file)
  
  barcode.annotation.files <- .getFilesFromResDirs(cellranger.dirs, barcode.file)
  
  sample.annotation.files <- .getFilesFromResDirs(cellranger.dirs, summary.file)
  
  if (is.null(regions)) {
    atac.sce <- createTileSCE(
      fragment.files,
      sample.names,
      output.dir = output.dir,
      tile.size = tile.size,
      seq.lengths = seq.lengths,
      barcodes.list = barcode.annotation.files,
      BPPARAM = BPPARAM
    )
  } else {
    atac.sce <- createRegionSCE(
      fragment.files,
      sample.names,
      regions = regions,
      output.dir = output.dir,
      barcodes.list = barcode.annotation.files,
      BPPARAM = BPPARAM
    )
    
  }
  
  
  
  rna.sce <- createRNASCE (
    filtered.feature.matrix.files,
    sample.names,
    feature.type = 'Gene Expression',
    filter.features.without.intervals = filter.features.without.intervals
  )
  
  # Create the MultiAssayExperiment from the list of Experiments
  mae <- createArbalistMAE(
    atac.sce,
    rna.sce,
    sample.annotation = sample.annotation.files,
    rna.name =  rna.name,
    atac.name = atac.name
  )
  
}



.getFilesFromResDirs <- function(res.dirs, file.name) {
  selected.files <- lapply(res.dirs, function(x, file.name) {
    potential.file <- file.path(x, file.name)
    if (file.exists(potential.file)) {
      return(potential.file)
    } else {
      warning(potential.file, " does not exist")
      potential.file <- NA
    }
  })
  selected.files <- unlist(selected.files)
}
