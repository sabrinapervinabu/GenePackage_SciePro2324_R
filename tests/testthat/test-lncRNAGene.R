# tests/testthat/test-lncRNAGene.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("LongNonCodingRNAGene class works correctly", {
  gr <- GRanges(seqnames = "chr11", ranges = IRanges(start = 100000, end = 110000))
  lncrna_gene <- new("LongNonCodingRNAGene",
                     id = "ENSG00000224207",
                     hugoSymbol = "HOTAIR",
                     name = "HOTAIR",
                     description = "HOX Transcript Antisense RNA",
                     structure = gr,
                     lncRNAID = "lncHOTAIR",
                     rnaSequence = "AUGCUCGAA...")
  
  expect_s4_class(lncrna_gene, "LongNonCodingRNAGene")
  expect_equal(lncrna_gene@hugoSymbol, "HOTAIR")
  expect_equal(lncrna_gene@lncRNAID, "lncHOTAIR")
})
