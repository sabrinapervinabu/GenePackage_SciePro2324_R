# GenePackage

**GenePackage** is an R package that defines an extensible set of S4 classes for representing and manipulating different types of genes (protein-coding, non-coding RNA, pseudogenes, etc.), along with their biological attributes and genomic structures.

---

## Installation
To install from a local `.tar` source:
```r
install.packages("GenePackage_1.0.0.tar", repos = NULL, type = "source")
# Load the package
library(GenePackage)
```
## Overview
The core functionality is based on a virtual base class Gene, from which specialized gene types inherit. These include:
- ProteinCodingGene
- LongNonCodingRNAGene
- MicroRNAGene
- SmallNuclearRNAGene
- RibosomalRNAGene
- TransferRNAGene
- EnhancerRNAGene
- Pseudogene
Each class contains specific slots to hold biologically relevant information (e.g., sequences, IDs, structure using GRanges).

### Usage Example
#### ProteinCodingGene
```r
library(GenePackage)
library(GenomicRanges)

# Create a simple ProteinCodingGene object
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 100, end = 500), strand = "+")
gene <- ProteinCodingGene(
  id = "ENSG000001", 
  hugoSymbol = "TP53", 
  name = "Tumor Protein P53",
  proteinID = "P04637",
  proteinSequence = "MENSDQE...",
  exonCount = 11L,
  description = "Tumor suppressor gene",
  structure = gr
)

# Access gene attributes
hugoSymbol(gene)
lengthProduct(gene)
```
#### LongNonCodingRNAGene
```r
# Create a LongNonCodingRNAGene object
lnc_gene <- LongNonCodingRNAGene(
  id = "ENSG000002", 
  hugoSymbol = "MALAT1", 
  name = "Metastasis Associated Lung Adenocarcinoma Transcript 1",
  lncRNAID = "lncMALAT1",
  rnaSequence = paste(rep("AUGCUU", 50), collapse = ""),  # synthetic sequence
  description = "A well-known lncRNA involved in metastasis",
  structure = gr
)

# Access HUGO symbol
hugoSymbol(lnc_gene)

# Calculate RNA sequence length
lengthProduct(lnc_gene)
```

### Requirements
This package depends on:
- GenomicRanges
- methods
Install them with:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
```

