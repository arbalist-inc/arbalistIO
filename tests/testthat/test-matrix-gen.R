test_that("saveTileMatrix creates correct output structure", {
  tmp_frag <- tempfile(fileext = ".fragments.gz")
  tmp_h5 <- tempfile(fileext = ".h5")
  seq.lengths <- c(chr1 = 10000, chr2 = 5000)

  mockFragmentFile(
    tmp_frag,
    seq.lengths,
    num.fragments = 500,
    cell.names = LETTERS[1:5]
  )

  expect_silent(
    res <- saveTileMatrix(
      fragment.file = tmp_frag,
      output.file = tmp_h5,
      output.name = "tiles",
      seq.lengths = seq.lengths,
      tile.size = 500
    )
  )

  expect_true(file.exists(tmp_h5))
  expect_s4_class(res$counts, "DelayedArray")
  expect_s4_class(res$tiles, "GRanges")
  expect_equal(ncol(res$counts), 5)

  # num of tiles: ceiling(10000/500) [for chr1] + ceiling(5000/500) [for chr2]
  expect_equal(length(res$tiles), 30)
})

test_that("saveRegionMatrix handles overlaps", {
  tmp_frag <- tempfile(fileext = ".fragments.gz")
  tmp_h5 <- tempfile(fileext = ".h5")
  seq.lengths <- c(chr1 = 2000)

  mockFragmentFile(
    tmp_frag,
    seq.lengths,
    num.fragments = 100,
    cell.names = c("CellA")
  )
  regions <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(1, 500), end = c(200, 700))
  )

  expect_silent(
    res <- saveRegionMatrix(
      fragment.file = tmp_frag,
      output.file = tmp_h5,
      output.name = "genes",
      regions = regions
    )
  )

  expect_s4_class(res, "DelayedArray")
  expect_equal(nrow(res), 2)
  expect_equal(ncol(res), 1)
})

test_that("saveTileMatrix fails on missing inputs", {
  expect_error(saveTileMatrix(
    fragment.file = "nonexistent.tsv.gz",
    output.file = "out.h5",
    output.name = "test"
  ))
})
