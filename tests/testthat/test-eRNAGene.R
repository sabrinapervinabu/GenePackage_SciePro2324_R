# tests/testthat/test-eRNAGene.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("EnhancerRNAGene class works correctly", {
  gr <- GRanges(seqnames = "chr8", ranges = IRanges(start = 150000, end = 151000))
  erna_gene <- new("EnhancerRNAGene",
                   id = "ENSG00000212349",
                   hugoSymbol = "eRNA",
                   name = "Enhancer RNA near MYC",
                   targetGeneID = "ENSG00000136997",
                   rnaSequence = "AUCGUAUC...",
                   description = "This eRNA is associated with the regulation of the MYC gene.",
                   structure = gr)
  
  expect_s4_class(erna_gene, "EnhancerRNAGene")
  expect_equal(erna_gene@hugoSymbol, "eRNA")
  expect_equal(erna_gene@targetGeneID, "ENSG00000136997")
})
