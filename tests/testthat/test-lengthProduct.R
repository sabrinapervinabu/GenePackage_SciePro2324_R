# tests/testthat/test-lengthProduct.R
library(testthat)
library(GenePackage)  
library(GenomicRanges)
test_that("lengthProduct function works correctly", {
  lncrna_gene <- new("LongNonCodingRNAGene",
                     id = "ENSG00000212350",
                     hugoSymbol = "XIST",
                     name = "XIST",
                     lncRNAID = "lncXIST",
                     rnaSequence = "AUGCUAGC",
                     description = "X-inactive specific transcript",
                     structure = GRanges(seqnames = "chrX", ranges = IRanges(start = 5000, end = 20000)))
  
  expect_equal(lengthProduct(lncrna_gene), 8L)  # Expected length of "AUGCUAGC" is 8
})
