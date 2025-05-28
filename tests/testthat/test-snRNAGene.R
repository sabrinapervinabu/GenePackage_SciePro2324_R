# tests/testthat/test-snRNAGene.R
library(testthat)
library(GenePackage)
library(GenomicRanges)
test_that("SmallNuclearRNAGene class works correctly", {
  gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 2000, end = 3000))
  snrna_gene <- new("SmallNuclearRNAGene",
                    id = "ENSG00000212345",
                    hugoSymbol = "U1",
                    name = "U1 small nuclear RNA",
                    snRNAID = "SNRNP70",
                    rnaSequence = "AUGCUAGC...",
                    description = "U1 snRNA is involved in the splicing of pre-mRNA.",
                    structure = gr)
  
  expect_s4_class(snrna_gene, "SmallNuclearRNAGene")
  expect_equal(snrna_gene@hugoSymbol, "U1")
  expect_equal(snrna_gene@snRNAID, "SNRNP70")
})
