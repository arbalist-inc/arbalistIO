test_that("filterDuplicateFeatures filters features correctly",{
  counts <- matrix(1:10, ncol=2)
  rdata <- DataFrame(name = c("GeneA", "GeneA", "GeneB", "GeneC", "GeneC"))
  se <- SummarizedExperiment(assays=list(counts=counts), rowData=rdata)
  
  # Filter duplicates, keeping the row with the max sum
  se_filtered <- filterDuplicateFeatures(se, mcol.name="name", selection.metric=max)
  filtered_matrix <- assay(se_filtered)
  rownames(filtered_matrix) <- mcols(se_filtered)$name
  colnames(filtered_matrix) <- c("V1","V2")
  
  expected_matrix <- aggregate(counts, by = list(rdata$name), max)
  rownames(expected_matrix) <- expected_matrix$Group.1
  expected_matrix <- as.matrix(expected_matrix[,-1])
  
  expect_equal(filtered_matrix, expected_matrix)
  
})

test_that("extractGenomeRefFilenameFromFragmentFile extracts genome file",{
  f <- tempfile(fileext=".tsv.gz")
  mockFragmentFile(f, c(chr1=1000), 10, LETTERS[1:5],
                   comments=c("reference_path=/path/to/ref"))
  
  out <- extractGenomeRefFilenameFromFragmentFile(f)
  expected <- "/path/to/ref/fasta/genome.fa.fai"
  expect_equal(out,expected)
  
})
