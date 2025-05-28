# tests/testthat/test-ProteinCodingGene.R
library(testthat)
library(GenePackage) 
library(GenomicRanges)
test_that("ProteinCodingGene class works correctly", {
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
  
  expect_s4_class(protein_gene, "ProteinCodingGene")  
  expect_equal(protein_gene@hugoSymbol, "BRCA2")      
  expect_equal(protein_gene@exonCount, 27L)            
})
