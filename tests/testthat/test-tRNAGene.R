# tests/testthat/test-tRNAGene.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("TransferRNAGene class works correctly", {
  gr <- GRanges(seqnames = "chr6", ranges = IRanges(start = 50000, end = 50500))
  trna_gene <- new("TransferRNAGene",
                   id = "ENSG00000212347",
                   hugoSymbol = "TRNA",
                   name = "Transfer RNA for Alanine",
                   tRNAID = "TRNAAla",
                   aminoAcid = "Alanine",
                   rnaSequence = "GGGUUACCC...",
                   description = "This tRNA carries alanine to the ribosome during translation.",
                   structure = gr)
  
  expect_s4_class(trna_gene, "TransferRNAGene")
  expect_equal(trna_gene@hugoSymbol, "TRNA")
  expect_equal(trna_gene@aminoAcid, "Alanine")
})
