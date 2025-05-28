# tests/testthat/test-miRNAGene.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("MicroRNAGene class works correctly", {
  gr <- GRanges(seqnames = "chr7", ranges = IRanges(start = 30000, end = 31000))
  mirna_gene <- new("MicroRNAGene",
                    id = "ENSG00000212351",
                    hugoSymbol = "MIR21",
                    name = "miR-21",
                    description = "MicroRNA 21, involved in cancer progression",
                    structure = gr,
                    miRNAID = "miR-21",
                    rnaSequence = "AAGCUUUC...",
                    seedSequence = "AAGCUU")
  
  expect_s4_class(mirna_gene, "MicroRNAGene")
  expect_equal(mirna_gene@hugoSymbol, "MIR21")
  expect_equal(mirna_gene@miRNAID, "miR-21")
})
