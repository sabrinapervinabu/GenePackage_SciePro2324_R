# tests/testthat/test-Pseudogene.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("Pseudogene class works correctly", {
  gr <- GRanges(seqnames = "chr11", ranges = IRanges(start = 300000, end = 305000))
  pseudo_gene <- new("Pseudogene",
                     id = "ENSG00000212348",
                     hugoSymbol = "ΨBRCA1",
                     name = "Pseudogene of BRCA1",
                     parentGeneID = "ENSG00000012048",
                     functionalStatus = "non-functional",
                     description = "This is a pseudogene derived from the BRCA1 gene.",
                     structure = gr)
  
  expect_s4_class(pseudo_gene, "Pseudogene")
  expect_equal(pseudo_gene@hugoSymbol, "ΨBRCA1")
  expect_equal(pseudo_gene@parentGeneID, "ENSG00000012048")
})
