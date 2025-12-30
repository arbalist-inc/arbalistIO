test_that("saveTileMatrix counts fragments correctly", {
# example 1: check for multiple fragments, multiple cells and multiple counts

# fragments files are 0-based, inclusive starts and exclusive ends
# chr1 0 100 cell1 1
# chr1 489 1005 cell1 1
# chr1 1 200 cell2 2

temp <- tempfile(fileext = ".gz")
temp.out <- tempfile(fileext=".h5")
frag.mock <- data.frame(chr=c("chr1", "chr1", "chr1"),
                        start=c(0,489,1),
                        end=c(100,1005,200),
                        cells=c("cell1","cell1","cell2"),
                        counts=c(1,1,2))

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

seq.lengths <- c(chr1=1500)
counted <- saveTileMatrix(temp, output.file=temp.out, output.name="example2", seq.lengths=seq.lengths)
counted


# expected

#ranges
#chr1 1-500 
#chr1 501-1000
#chr1 1001-1500

#counts
#   cell1 cell2
#1.   2     1
#2.   0     0 
#3.   1     0



expected_counts <- matrix(as.raw(c(2,0,1,1,0,0)), ncol=2)
colnames(expected_counts) <- c("cell1","cell2")

expect_identical(as.matrix(counted$counts), expected = expected_counts)



# example 2: check for boundary cases

# fragments
# chr1 0 499 cell1 1
# chr1 0 500 cell2 1
# chr1 0 501 cell3 1
# chr1 0 502 cell4 1
# chr1 998 1005 cell5 1
# chr1 999 1005 cell6 1
# chr1 1000 1005 cell7 1

temp <- tempfile(fileext = ".gz")
temp.out <- tempfile(fileext=".h5")
frag.mock <- data.frame(chr=rep("chr1",7),
                        start=c(0,0,0,0,998,999,1000),
                        end=c(499,500,501,502,1005,1005,1005),
                        cells=c("cell1","cell2","cell3", "cell4", "cell5","cell6","cell7"),
                        counts=rep(1,7))

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

seq.lengths <- c(chr1=1500)
counted <- saveTileMatrix(temp, output.file=temp.out, output.name="example2", seq.lengths=seq.lengths)
counted


# expected

#ranges
#chr1 1-500 
#chr1 501-1000
#chr1 1001-1500

#counts
#   cell1 cell2 cell3 cell4 cell5 cell6 cell7
#1.   1     1    1     1      0    0     0
#2.   0     0    1     1      1    1     0
#3.   0     0    0     0      1    1     1



expected_counts <- matrix(as.raw(c(1,0,0,
                                   1,0,0,
                                   1,1,0,
                                   1,1,0,
                                   0,1,1,
                                   0,1,1,
                                   0,0,1)), 
                          ncol=7)
colnames(expected_counts) <- paste0("cell",1:7)

expect_identical(as.matrix(counted$counts), expected = expected_counts)
})