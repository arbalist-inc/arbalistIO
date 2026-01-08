test_that("saveRegionMatrix detects unsorted fragments", {
  # Fragments are not sorted by start position within the same chromosome
  # chr1 100 200
  # chr1 50 150
  
  temp <- tempfile(fileext = ".gz")
  frag.mock <- data.frame(
    chr = c("chr1", "chr1"),
    start = c(100, 50),
    end = c(200, 150),
    cells = c("cell1", "cell2"),
    counts = c(1, 1)
  )
  
  handle <- gzfile(temp, "w")
  write.table(frag.mock, file = handle, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  close(handle)
  
  regions <- GRanges("chr1:1-1000")
  
  expect_error(
    saveRegionMatrix(temp, output.file = tempfile(), output.name = "out", regions = regions),
    "fragment start .* is less than the previous fragment start"
  )
})

test_that("saveRegionMatrix detects unordered sequences", {
  
  temp <- tempfile(fileext = ".gz")
  frag.mock <- data.frame(
    chr = c("chr1", "chr2", "chr1"),
    start = c(100, 100, 200),
    end = c(200, 200, 300),
    cells = c("c1", "c1", "c1"),
    counts = c(1, 1, 1)
  )
  
  handle <- gzfile(temp, "w")
  write.table(frag.mock, file = handle, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  close(handle)
  
  regions <- GRanges(c("chr1:1-1000", "chr2:1-1000"))
  
  expect_error(
    saveRegionMatrix(temp, output.file = tempfile(), output.name = "out", regions = regions),
    "unordered sequence names"
  )
})

test_that("saveTileMatrix detects unsorted fragments", {
  # Fragments are not sorted by start position
  
  temp <- tempfile(fileext = ".gz")
  frag.mock <- data.frame(
    chr = c("chr1", "chr1"),
    start = c(100, 50),
    end = c(200, 150),
    cells = c("cell1", "cell2"),
    counts = c(1, 1)
  )
  
  handle <- gzfile(temp, "w")
  write.table(frag.mock, file = handle, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  close(handle)
  
  seq.lengths <- c(chr1 = 1000)
  
  expect_error(
    saveTileMatrix(temp, output.file = tempfile(), output.name = "out", seq.lengths = seq.lengths),
    "fragment start .* is less than the previous fragment start"
  )
})

test_that("saveTileMatrix detects sequence order mismatch", {
  # This is already checked by cpp code
  
  temp <- tempfile(fileext = ".gz")
  frag.mock <- data.frame(
    chr = c("chr2", "chr1"),
    start = c(100, 100),
    end = c(200, 200),
    cells = c("c1", "c1"),
    counts = c(1, 1)
  )
  
  handle <- gzfile(temp, "w")
  write.table(frag.mock, file = handle, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  close(handle)
  
  seq.lengths <- c(chr1 = 1000, chr2 = 1000)
  
  expect_error(
    saveTileMatrix(temp, output.file = tempfile(), output.name = "out", seq.lengths = seq.lengths),
    "order of sequences in fragment file differs"
  )
})