test_that("saveRegionMatrix counts fragments correctly", {
  
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
#chr1 2-50 
#chr1 490-1000

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
#chr1 3-50 
#chr1 490-1000

# expected counts
#   cell1 cell2
#1.   0     0
#2.   0     0

expected_counts <- matrix(as.raw(c(0, 0, 0, 0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)

# test start regions no overlap
regions <- GRanges(c("chr1:10-198","chr1:500-1003"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example3", regions=regions)
counted

# ranges
#chr1 10-198 
#chr1 500-1003

# expected counts
#   cell1 cell2
#1.   0     0
#2.   0     0

expected_counts <- matrix(as.raw(c(0, 0, 0, 0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


# test end regions overlap
regions <- GRanges(c("chr1:10-199","chr1:500-1004"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example4", regions=regions)
counted

# ranges
#chr1 10-199 
#chr1 500-1004

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
#chr1 chr1:401-1000
#chr1 chr1:1001-1500

# expected counts
#   cell1 cell2
#1.   0     1
#2.   0     1

# test overlapping regions no overlap
regions <- GRanges(c("chr1:800-1100", "chr1:900-1500"))
counted <- saveRegionMatrix(temp, output.file=temp.out, output.name="example6", regions=regions)
counted

# ranges
#chr1 chr1:1000-1100
#chr1 chr1:900-1500

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
#chr1 chr1:890-950
#chr1 chr1:900-1500

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
#chr1 chr1:800-1500
#chr1 chr1:900-1500

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
#chr1 chr1:800-1500
#chr1 chr1:1100-1500

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
#chr1 chr1:1100-1500
#chr1 chr1:800-1500

# expected counts
#   cell1 cell2
#1.   0     0     # since all of it overlapped with 800-1500
#2.   0     1

expected_counts <- matrix(as.raw(c(0, 0, 0, 1)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")
expect_identical(as.matrix(counted), expected = expected_counts)


})