test_that("createRNASCE parses H5 along with interval filtering", {
  tmp_h5 <- tempfile(fileext = ".h5")
  mockCellRangerH5(tmp_h5, n.genes = 10, n.cells = 5, cell.names = LETTERS[1:5])

  sce <- suppressWarnings(
    createRNASCE(
      h5.files = tmp_h5,
      sample.names = "Sample1",
      filter.features.without.intervals = TRUE
    )
  )

  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(nrow(sce), 5)
  expect_true("Sample" %in% names(colData(sce)))
  expect_equal(unique(colData(sce)$Sample), "Sample1")

  sce_full <- suppressWarnings(
    createRNASCE(
      h5.files = tmp_h5,
      sample.names = "Sample1",
      filter.features.without.intervals = FALSE
    )
  )
  expect_equal(nrow(sce_full), 10)
})

test_that("createRNASCE handles multiple samples", {
  tmp1 <- tempfile(fileext = ".h5")
  mockCellRangerH5(tmp1, n.genes = 10, n.cells = 5, cell.names = LETTERS[1:5])
  tmp2 <- tempfile(fileext = ".h5")
  mockCellRangerH5(tmp2, n.genes = 10, n.cells = 10, cell.names = LETTERS[1:10])

  sce <- suppressWarnings(
    createRNASCE(
      h5.files = c(tmp1, tmp2),
      sample.names = c("S1", "S2")
    )
  )

  expect_equal(ncol(sce), 5+10)
  expect_equal(unique(colData(sce)$Sample), c("S1", "S2"))
  
})
