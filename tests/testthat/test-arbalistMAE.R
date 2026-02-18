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

  #expect_true("fragment_file" %in% names(colData(mae)))
  #expect_equal(as.vector(colData(mae)[1, "fragment_file"]), frag_file)
})

# test_that("createArbalistMAE validates input lengths", {
#   expect_error(
#     createArbalistMAE(
#       sample.names = c("A", "B"),
#       fragment.files = c("file1"),
#       filtered.feature.matrix.files = c("file1", "file2")
#     ),
#     "fragment.files and sample names need to be the same length"
#   )
# })
