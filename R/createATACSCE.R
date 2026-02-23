#' Create a SingleCellExperiment from fragment files
#'
#' Create a SingleCellExperiment from fragment files, storing the tile matrix.
#'
#' @param fragment.files Character vector specifying fragment files
#' @param sample.names Character vector specifying sample names corresponding
#'     to fragment.files
#' @param output.dir String containing the directory where files should be
#'   output
#' @param tile.size Integer scalar specifying the size of the tiles in base
#'   pairs
#' @param seq.lengths Named integer vector containing the lengths of the
#'   reference sequences used for alignment. Vector names should correspond to
#'   the names of the sequences, in the same order of occurrence as in the
#'   fragment file. If \code{NULL}, this is obtained from the reference genome
#'   used by Cellranger (itself located by scanning the header of the fragment
#'   file). 
#' @param barcodes.list A named list with samples as names and each list element
#'     as a character vector of barcodes. If \code{NULL}, all barcodes are extracted
#' @param matrix.name Character string indicating the name of the matrix
#' @param BPPARAM A \link[BiocParallel]{BiocParallelParam} object indicating how matrix
#'   creation should be parallelized.
#'
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} containing the tile matrix.
#'
#' @examples
#' temp1 <- tempfile(fileext=".fragments.gz")
#' mockFragmentFile(temp1, c(chrA=1000, chrB=200000, chrC=200),
#'                  num.fragments=100, cell.names=LETTERS)
#' temp2 <- tempfile(fileext=".fragments.gz")
#' mockFragmentFile(temp2, c(chrA=1000, chrB=200000, chrC=200),
#'                  num.fragments=100, cell.names=LETTERS)
#' tile_sce <- createTileSCE(fragment.files=c(temp1,temp2), 
#'                           sample.names=c("sample1", "sample2"),
#'                           tile.size=500, 
#'                           seq.lengths = c(chrA=1000, chrB=200000, chrC=200))
#' file.remove(list.files(tempdir(), pattern=".h5", full.names = TRUE))
#' @author Natalie Fox, Jayaram Kancherla, Xiaosai Yao
#' @importFrom BiocParallel bpparam
#' @importFrom stats na.omit
#' @export
createTileSCE <- function(fragment.files,
                          sample.names, 
                          output.dir = tempdir(),
                          tile.size = 500,
                          seq.lengths = NULL,
                          barcodes.list = NULL,
                          matrix.name = "TileMatrix",
                          BPPARAM = bpparam()) {
  
  names(fragment.files) <- sample.names
  fragment.files <- na.omit(fragment.files)
  
  # check that fragment headers contained required information
  if (is.null(seq.lengths)) {
    info <- .processFragmentHeader(fragment.files[1])
    if (!"reference_path" %in% names(info)) {
      stop(
        "the fragment file header does not have reference_path information \
        so please specify the seq.lengths argument"
      )
    }
  }
  
  createSCEFromFragments(
    fragment.files = fragment.files,
    output.dir = output.dir,
    matrix.name = paste0(matrix.name, tile.size),
    worker.fun = .saveTileMatrixCall,
    BPPARAM = BPPARAM,
    tile.size = tile.size,
    seq.lengths = seq.lengths,
    barcodes.list = barcodes.list
  )
}

#' Create a SingleCellExperiment for Gene Scores from fragment files
#'
#' Create a SingleCellExperiment from fragment files, storing the gene score matrix.
#'
#' @inheritParams createTileSCE 
#' @param regions Genomic Ranges specifying gene coordinates for creating the
#'   gene score matrix.
#' @importFrom BiocParallel bpparam
#' @importFrom GenomicRanges GRanges
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} containing the gene score matrix.
#' @examples
#' temp1 <- tempfile(fileext=".fragments.gz")
#' mockFragmentFile(temp1, c(chrA=1000, chrB=200000, chrC=200),
#'                  num.fragments=100, cell.names=LETTERS)
#' temp2 <- tempfile(fileext=".fragments.gz")
#' mockFragmentFile(temp2, c(chrA=1000, chrB=200000, chrC=200),
#'                  num.fragments=100, cell.names=LETTERS)
#' fragment_files <- c(sample1=temp1, sample2=temp2)
#' region_sce <- createRegionSCE(fragment.files=c(temp1,temp2),
#'                               sample.names=c("sample1", "sample2"),                                
#'                               region=GenomicRanges::GRanges(c("chrA:500-1000", "chrB:1000-2000")))
#' file.remove(list.files(tempdir(), pattern=".h5", full.names = TRUE))
#' @author Natalie Fox, Jayaram Kancherla, Xiaosai Yao
#' @export
#' 
#' 

createRegionSCE <- function(fragment.files,
                            sample.names,
                            output.dir = tempdir(),
                            regions,
                            barcodes.list = NULL,
                            matrix.name = "GeneAccessibilityMatrix",
                            BPPARAM = bpparam()) {
  names(fragment.files) <- sample.names
  
  fragment.files <- na.omit(fragment.files)
  
  createSCEFromFragments(
    fragment.files = fragment.files,
    output.dir = output.dir,
    matrix.name = matrix.name,
    worker.fun = .saveRegionMatrixCall,
    BPPARAM = BPPARAM,
    regions = regions,
    barcodes.list = barcodes.list
  )
}



createSCEFromFragments <- function(fragment.files,
                                   output.dir,
                                   matrix.name,
                                   worker.fun,
                                   BPPARAM = bpparam(),
                                   ...) {
  # check if the hdf5 files already exist
  output.file.names <- file.path(output.dir, 
                                 paste0(matrix.name,'_', names(fragment.files),'.h5'))
  
  if (any(table(output.file.names) > 1)) {
    # if two file names are the same then add a random component to the file name.
    output.file.names <- tempfile(
      pattern = paste0(matrix.name, "_", names(fragment.files), "_"),
      tmpdir = output.dir,
      fileext = ".h5"
    )
  }
  names(output.file.names) <- names(fragment.files)
  if (any(file.exists(output.file.names))) {
    stop(output.file.names[which(file.exists(output.file.names))[1]],
         ' already exists. We do not want to overwrite the file in case it is being used. 
        Either remove the file if you think it is safe to do so or specify a different output.dir.'
    )
  }
  
  # Parallelizing per sample
  res.list <- bptry(
    bplapply(
      seq_along(fragment.files),
      worker.fun,
      fragment.files = fragment.files,
      output.file.names = output.file.names,
      BPPARAM = BPPARAM,
      ...
    )
  )
  
  # Extract counts (H5SparseMatrix) from results
  # Use lapply to preserve list structure and names
  tile.res.list <- lapply(res.list, function(x) {
    x$counts
  })
  names(tile.res.list) <- names(fragment.files)
  
  # Use the first sample's ranges as the usage
  tile.grs <- res.list[[1]]$tiles
  
  # check that the tiles are the same for all samples
  for (i in setdiff(seq_along(res.list), 1)) {
    if (length(tile.grs) != length(res.list[[i]]$tiles) ||
        !all(tile.grs == res.list[[i]]$tiles)) {
      stop("Matrix GRanges do not match")
    }
  }
  
  sce <- getSCEFromH5List(tile.res.list, tile.grs)
  mainExpName(sce) <- matrix.name
  
  return(sce)
}


#' Create a SingleCellExperiment from a list of delayed matrices
#'
#' Create a SingleCellExperiment from a list of delayed matrices using
#' AmalgamatedArray.
#'
#' @param h5.res.list List containing delayed matrices with HDF5 backends that
#'   will be combined using \link[alabaster.matrix]{AmalgamatedArray} into a 
#'   \link[SingleCellExperiment]{SingleCellExperiment} object. List
#'   item names should be the sample name for the delayed matrix.
#' @param grs \link[GenomicRanges]{GRanges} object to be used for the rowRanges 
#'    of the resulting \link[SingleCellExperiment]{SingleCellExperiment} object
#'
#' @return A SingleCellExperiment
#'
#' @author Natalie Fox
#' @importFrom alabaster.matrix AmalgamatedArray
#' @importFrom SingleCellExperiment SingleCellExperiment mainExpName<-
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges<- colData<-
#' @importFrom methods as
#' @importFrom BiocParallel bptry bplapply bpparam

getSCEFromH5List <- function(h5.res.list, grs) {
  # Combined the per sample results into one matrix
  if (length(h5.res.list) == 1) {
    mat <- h5.res.list[[1]]
  } else {
    mat <- AmalgamatedArray(h5.res.list, along = 2)
  }
  
  # Map the cells back to samples and update colnames to include sample names
  cell.to.sample <- unlist(lapply(names(h5.res.list), function(x) {
    rep(x, ncol(h5.res.list[[x]]))
  }))
  
  # Create a SingleCellExperiment
  # Ensure colnames are not NULL to avoid paste0 issues
  cnames <- colnames(mat)
  if (is.null(cnames)) {
    cnames <- character(ncol(mat))
  }
  new.cnames <- paste0(cell.to.sample, "#", cnames)
  
  # Set dimnames on matrix directly
  dimnames(mat) <- list(NULL, new.cnames)
  
  mat.list <- list(counts = mat)
  
  se <- SummarizedExperiment(mat.list, rowRanges = grs)
  colData(se)$Sample <- as.character(cell.to.sample)
  sce <- as(se, "SingleCellExperiment")
  
  return(sce)
}

.saveTileMatrixCall <- function(sample.name,
                                fragment.files,
                                output.file.names,
                                tile.size,
                                seq.lengths,
                                barcodes.list = NULL) {
  
  barcodes <- NULL
  if (!is.null(barcodes.list)) {
    barcodes <- barcodes.list[[sample.name]]
  }
  if (!is.null(barcodes)) {
    barcodes <- as.character(barcodes)
  }
  
  tile.res <- saveTileMatrix(
    as.character(fragment.files[sample.name]),
    output.file = as.character(output.file.names[sample.name]),
    output.name = "tile_matrix",
    tile.size = tile.size,
    seq.lengths = seq.lengths,
    barcodes = barcodes
  )
  return(tile.res)
}

.saveRegionMatrixCall <- function(sample.name,
                                  fragment.files,
                                  output.file.names,
                                  regions,
                                  barcodes.list = NULL) {
  
  barcodes <- NULL
  if (!is.null(barcodes.list)) {
    barcodes <- barcodes.list[[sample.name]]
  }
  if (!is.null(barcodes)) {
    barcodes <- as.character(barcodes)
  }
  
  matrix.res <- saveRegionMatrix(
    as.character(fragment.files[sample.name]),
    output.file = as.character(output.file.names[sample.name]),
    output.name = "gene_matrix",
    regions = regions,
    barcodes = barcodes
  )
  return(list(counts = matrix.res, tiles = regions))
}

