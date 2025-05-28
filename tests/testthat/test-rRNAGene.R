# tests/testthat/test-rRNAGene.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("RibosomalRNAGene class works correctly", {
  gr <- GRanges(seqnames = "chr13", ranges = IRanges(start = 100000, end = 101000))
  rrna_gene <- new("RibosomalRNAGene",
                   id = "ENSG00000212346",
                   hugoSymbol = "18S",
                   name = "18S ribosomal RNA",
                   rRNAID = "RRNA18S",
                   rnaSequence = "AGCUAGCU...",
                   description = "18S rRNA is a component of the small ribosomal subunit.",
                   structure = gr)
  
  expect_s4_class(rrna_gene, "RibosomalRNAGene")
  expect_equal(rrna_gene@hugoSymbol, "18S")
  expect_equal(rrna_gene@rRNAID, "RRNA18S")
})
