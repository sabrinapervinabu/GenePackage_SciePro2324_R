# tests/testthat/test-hugoSymbol.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("hugoSymbol function works correctly", {
  gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1000, end = 5000))
  protein_gene <- new("ProteinCodingGene",
                      id = "ENSG00000139618",
                      hugoSymbol = "BRCA2",
                      name = "BRCA2",
                      description = "Breast cancer type 2 susceptibility protein",
                      structure = gr,
                      proteinID = "NP_000050.2",
                      proteinSequence = "MSSQLQALFQS...",
                      exonCount = 27L)
  
  expect_equal(hugoSymbol(protein_gene), "BRCA2")
  
  hugoSymbol(protein_gene) <- "NEW_SYMBOL"
  expect_equal(hugoSymbol(protein_gene), "NEW_SYMBOL")
})
