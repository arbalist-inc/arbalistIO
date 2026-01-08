test_that("Test .sanitizeRegions work correctly",{
  
  # disjoint regions
  regions <- GRanges(c("chr1:1-100","chr1:101-200","chr1:201-300"))
  solo <- .extractNonOverlaps(regions)
  sanitized <- .prepareRegions(solo, regions)
  
  expected_sanitized <- list()
  expected_sanitized[["ids"]][["chr1"]] <- seq_len(length(regions)) - 1L
  expected_sanitized[["starts"]][["chr1"]] <- start(regions) - 1L
  expected_sanitized[["ends"]][["chr1"]] <- end(regions)
  
  expect_identical(sanitized, expected_sanitized)
  
  
  
  # overlapping regions
  regions <- GRanges(c("chr1:1-100","chr1:50-200"))
  solo <- .extractNonOverlaps(regions)
  sanitized <- .prepareRegions(solo, regions)
  
  region1_unique <- setdiff(regions[1],regions[2])
  region2_unique <- setdiff(regions[2],regions[1])
  
  expected_sanitized <- list()
  expected_sanitized[["ids"]][["chr1"]] <- c(rep(1L, length(region1_unique)) - 1L, rep(2L, length(region2_unique))-1L)
  expected_sanitized[["starts"]][["chr1"]] <- c(start(region1_unique) - 1L, start(region2_unique) - 1L)
  expected_sanitized[["ends"]][["chr1"]] <- c(end(region1_unique), end(region2_unique))
  
  expect_identical(sanitized, expected_sanitized)
  
  

  
})


test_that("Test saveRegionMatrix counts fragments correctly", {
  
# fragments files are 0-based, inclusive starts and exclusive ends
# chr1 1 200 cell1 2
# chr1 489 1005 cell2 1


temp <- tempfile(fileext = ".gz")
temp.out <- tempfile(fileext=".h5")
frag.mock <- data.frame(chr=c("chr1", "chr1"),
                        start=c(1, 489),
                        end=c(200, 1005),
                        cells=c("cell1","cell2"),
                        counts=c(2,1))


frag.mock

handle <- gzfile(temp, "w")
write.table(
  file = handle,
  frag.mock,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t",
  quote = FALSE,
  eol = "\n"
)

close(handle)

# regions are 1-based granges with inclusive starts and ends
# test start regions overlap. 

regions <- GRanges(c("chr1:2-50","chr1:490-1000"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example1", regions=regions)
counted


# ranges
#chr1:2-50 
#chr1:490-1000

# expected counts
#   cell1 cell2
#1.   1     0
#2.   0     1

expected_counts <- matrix(as.raw(c(1, 0, 0, 1)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)

# test start regions no overlap
regions <- GRanges(c("chr1:3-50","chr1:491-1000"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example2", regions=regions)
counted

# ranges
#chr1:3-50 
#chr1:491-1000

# expected counts
#   cell1 cell2
#1.   0     0
#2.   0     0

expected_counts <- matrix(as.raw(c(0, 0, 0, 0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)

# test start regions no overlap
regions <- GRanges(c("chr1:10-199","chr1:500-1004"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example3", regions=regions)
counted

# ranges
#chr1:10-199 
#chr1:500-1004

# expected counts
#   cell1 cell2
#1.   0     0
#2.   0     0

expected_counts <- matrix(as.raw(c(0, 0, 0, 0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


# test end regions overlap
regions <- GRanges(c("chr1:10-200","chr1:500-1005"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example4", regions=regions)
counted

# ranges
#chr1:10-200 
#chr1:500-1005

# expected counts
#   cell1 cell2
#1.   1     0
#2.   0     1

expected_counts <- matrix(as.raw(c(1, 0, 0, 1)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)



# test contiguous regions
regions <- GRanges(c("chr1:401-1000", "chr1:1001-1500"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example5", regions=regions)
counted


# ranges
#chr1:401-1000
#chr1:1001-1500

# expected counts
#   cell1 cell2
#1.   0     1
#2.   0     1

# test overlapping regions no overlap
regions <- GRanges(c("chr1:800-1100", "chr1:900-1500"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example6", regions=regions)
counted

# ranges
#chr1:800-1100
#chr1:900-1500

# expected counts
#   cell1 cell2
#1.   0     0
#2.   0     0

expected_counts <- matrix(as.raw(c(0, 0, 0, 0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


# test overlapping regions overlap
regions <- GRanges(c("chr1:800-950", "chr1:900-1500"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example7", regions=regions)
counted

# ranges
#chr1:800-950
#chr1:900-1500

# expected counts
#   cell1 cell2
#1.   0     0
#2.   0     1

expected_counts <- matrix(as.raw(c(0, 0, 0, 1)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


# test nested regions no overlap
regions <- GRanges(c("chr1:800-1500", "chr1:900-1500"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example8", regions=regions)
counted

# ranges
#chr1:800-1500
#chr1:900-1500

# expected counts
#   cell1 cell2
#1.   0     0
#2.   0     0

expected_counts <- matrix(as.raw(c(0, 0, 0, 0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


# test nested regions overlap
regions <- GRanges(c("chr1:800-1500", "chr1:1100-1500"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example9", regions=regions)
counted

# ranges
#chr1:800-1500
#chr1:1100-1500

# expected counts
#   cell1 cell2
#1.   0     1
#2.   0     0

expected_counts <- matrix(as.raw(c(0, 0, 1, 0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


# test sorting of regions
regions <- GRanges(c("chr1:1100-1500", "chr1:800-1500"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example10", regions=regions)
counted

# ranges
#chr1:1100-1500
#chr1:800-1500

# expected counts
#   cell1 cell2
#1.   0     0     # since all of it overlapped with 800-1500
#2.   0     1

expected_counts <- matrix(as.raw(c(0, 0, 0, 1)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


# GRangesList
regions <- GRangesList(GRanges(c("chr1:1-100","chr1:101-200","chr1:201-300")), GRanges(c("chr1:400-500")))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example11", regions=regions)
counted

unlisted <- unlist(reduce(regions))

expected_counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example12", regions=unlisted)
expect_identical(as.matrix(counted), as.matrix(expected_counted))

})