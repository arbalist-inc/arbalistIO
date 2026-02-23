library(MultiAssayExperiment)

test_that("createArbalistMAE works", {
  tmp_dir <- tempdir()
  sample_name <- "TestSample"

  frag_file <- tempfile(
    pattern = "fragments",
    tmpdir = tmp_dir,
    fileext = ".tsv.gz"
  )
  seq_lengths <- c(chr1 = 1000)
  mockFragmentFile(
    frag_file,
    seq_lengths,
    num.fragments = 100,
    cell.names = LETTERS[1:5]
  )

  feat_file <- paste0(tmp_dir, "/", "filtered_feature_bc_matrix.h5")
  mockCellRangerH5(feat_file, n.genes = 10, n.cells = 5, cell.names = LETTERS[1:5])

  gene_grs <- GRanges("chr1", IRanges(10, 100))

  region.sce <- createRegionSCE(fragment.files=frag_file,
                                sample.names=sample_name,
                                region=gene_grs,
                                output.dir = tmp_dir)

  rna.sce <- createRNASCE(feat_file, sample.names=sample_name)
  
  mae <- createArbalistMAE(region.sce,
                           rna.sce,
                           sample.annotation = NULL,
                           rna.name = "GeneExpressionMatrix",
                           atac.name = "GeneAccessibilityMatrix")
  


  expect_s4_class(mae, "MultiAssayExperiment")

  exp_names <- names(experiments(mae))
  expect_true("GeneAccessibilityMatrix" %in% exp_names)
  expect_true("GeneExpressionMatrix" %in% exp_names)

})

test_that("createMAEFromCellranger works", {
  tmp_dir <- tempdir()
  file.remove(list.files(tmp_dir, pattern=".h5", full.names = TRUE))
  temp_rna <- file.path(tmp_dir,"filtered_feature_bc_matrix.h5")
  mockCellRangerH5(temp_rna, cell.names = LETTERS)
  temp_atac <- file.path(tmp_dir,"atac_fragments.tsv.gz")
  mockFragmentFile(temp_atac,
                   c(chrA=1000, chrB=200000, chrC=200),
                   num.fragments=100,
                   cell.names=LETTERS)
  seq.lengths <-  c(chrA=1000, chrB=200000, chrC=200)
  # create combined MAE
  mae <- createMAEFromCellranger(cellranger.dirs = tmp_dir,
                                 sample.names = "mock",
                                 seq.lengths = seq.lengths,
                                 consistent.barcodes = FALSE)
  file.remove(list.files(tmp_dir, pattern=".h5", full.names = TRUE))
  
  expect_s4_class(mae, "MultiAssayExperiment")
  
  exp_names <- names(experiments(mae))
  expect_true("TileMatrix500" %in% exp_names)
  expect_true("GeneExpressionMatrix" %in% exp_names)
})


