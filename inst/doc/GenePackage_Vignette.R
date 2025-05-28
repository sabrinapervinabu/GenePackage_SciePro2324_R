## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(GenePackage)

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}
library(GenomicRanges)

## ----example_ProteinCodingGene------------------------------------------------
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1000, end = 5000))
protein_gene <- new("ProteinCodingGene",
                    id = "ENSG00000139618",
                    hugoSymbol = "BRCA2",
                    name = "BRCA2",
                    description = "Breast cancer type 2 susceptibility protein",
                    structure = gr,
                    proteinID = "NP_000050.2",
                    proteinSequence = "MSSQLQALFQS...",
                    exonCount = as.integer(27))
protein_gene

## ----example_lncRNA-----------------------------------------------------------
gr <- GRanges(seqnames = "chr11", ranges = IRanges(start = 100000, end = 110000))
lncrna_gene <- new("LongNonCodingRNAGene",
                   id = "ENSG00000224207",
                   hugoSymbol = "HOTAIR",
                   name = "HOTAIR",
                   description = "HOX Transcript Antisense RNA",
                   structure = gr,
                   lncRNAID = "lncHOTAIR",
                   rnaSequence = "AUGCUCGAA...")
lncrna_gene

## ----example_miRNA------------------------------------------------------------
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
mirna_gene

## ----example_snRNA------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 2000, end = 3000))
snrna_gene <- SmallNuclearRNAGene(
  id = "ENSG00000212345",
  hugoSymbol = "U1",
  name = "U1 small nuclear RNA",
  snRNAID = "SNRNP70",
  rnaSequence = "AUGCUAGC...",
  description = "U1 snRNA is involved in the splicing of pre-mRNA.",
  structure = gr
)
snrna_gene

## ----example_rRNA-------------------------------------------------------------
gr <- GRanges(seqnames = "chr13", ranges = IRanges(start = 100000, end = 101000))
rrna_gene <- RibosomalRNAGene(
  id = "ENSG00000212346",
  hugoSymbol = "18S",
  name = "18S ribosomal RNA",
  rRNAID = "RRNA18S",
  rnaSequence = "AGCUAGCU...",
  description = "18S rRNA is a component of the small ribosomal subunit.",
  structure = gr
)
rrna_gene

## ----example_tRNA-------------------------------------------------------------
gr <- GRanges(seqnames = "chr6", ranges = IRanges(start = 50000, end = 50500))
trna_gene <- TransferRNAGene(
  id = "ENSG00000212347",
  hugoSymbol = "TRNA",
  name = "Transfer RNA for Alanine",
  tRNAID = "TRNAAla",
  aminoAcid = "Alanine",
  rnaSequence = "GGGUUACCC...",
  description = "This tRNA carries alanine to the ribosome during translation.",
  structure = gr
)
trna_gene

## ----example_pseudogenes------------------------------------------------------
gr <- GRanges(seqnames = "chr11", ranges = IRanges(start = 300000, end = 305000))
pseudo_gene <- Pseudogene(
  id = "ENSG00000212348",
  hugoSymbol = "ΨBRCA1",
  name = "Pseudogene of BRCA1",
  parentGeneID = "ENSG00000012048",
  functionalStatus = "non-functional",
  description = "This is a pseudogene derived from the BRCA1 gene.",
  structure = gr
)
pseudo_gene

## ----example_eRNA-------------------------------------------------------------
gr <- GRanges(seqnames = "chr8", ranges = IRanges(start = 150000, end = 151000))
erna_gene <- EnhancerRNAGene(
  id = "ENSG00000212349",
  hugoSymbol = "eRNA",
  name = "Enhancer RNA near MYC",
  targetGeneID = "ENSG00000136997",
  rnaSequence = "AUCGUAUC...",
  description = "This eRNA is associated with the regulation of the MYC gene.",
  structure = gr
)
erna_gene

## ----example_getter-----------------------------------------------------------
hugoSymbol(protein_gene)

## ----example_setter-----------------------------------------------------------
hugoSymbol(protein_gene) <- "NEW_SYMBOL"
hugoSymbol(protein_gene)

## ----example_lengthProduct----------------------------------------------------
lncrna_gene <- LongNonCodingRNAGene(
  id = "ENSG00000212350",
  hugoSymbol = "XIST",
  name = "XIST",
  lncRNAID = "lncXIST",
  rnaSequence = "AUGCUAGC...",
  description = "X-inactive specific transcript",
  structure = GRanges(seqnames = "chrX", ranges = IRanges(start = 5000, end = 20000))
)
lengthProduct(lncrna_gene)

## ----session-info, echo=FALSE-------------------------------------------------
sessionInfo()

