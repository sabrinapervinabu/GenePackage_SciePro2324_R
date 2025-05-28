#' @importFrom GenomicRanges GRanges
#' @importFrom methods new


library(GenomicRanges) 

#### gene class ####
#' Class Gene: a class to represent a generic gene with its attributes.
#' @slot id Gene ID (character)
#' @slot hugoSymbol HUGO Symbol (character)
#' @slot name Gene name (character)
#' @slot description Description of the gene (character)
#' @slot structure Gene structure, typically as a `GRanges` object
#' @exportClass Gene
setClass("Gene",
         slots = c(id = "character",
                   hugoSymbol = "character", 
                   name = "character",         
                   description = "character",  
                   structure = "GRanges"),
         prototype = list(description = NA_character_,
                          structure = GRanges()),
         contains = "VIRTUAL"
)

#### set of S4 classes for different gene types ####
#' Class "ProteinCodingGene"
#' @slot proteinID A unique identifier for the protein product encoded by the gene (character).
#' @slot proteinSequence The amino acid sequence of the protein encoded by the gene (character).
#' @slot exonCount The number of exons in the gene (integer).
#' @exportClass ProteinCodingGene
#' @family Gene classes
setClass("ProteinCodingGene",
         contains = "Gene", 
         slots = c(proteinID = "character", 
                   proteinSequence = "character", 
                   exonCount = "integer")
)

#' Class "LongNonCodingRNAGene"
#' @slot lncRNAID A unique identifier for the long non-coding RNA (lncRNA) product (character).
#' @slot rnaSequence The nucleotide sequence of the RNA encoded by the lncRNA gene (character).
#' @exportClass LongNonCodingRNAGene
#' @family Gene classes
setClass("LongNonCodingRNAGene",
         contains = "Gene",
         slots = c(lncRNAID = "character",
                   rnaSequence = "character")
)

#' Class "MicroRNAGene"
#' @slot miRNAID A unique identifier for the microRNA (miRNA) product (character).
#' @slot rnaSequence The nucleotide sequence of the RNA encoded by the miRNA gene (character).
#' @slot seedSequence The seed sequence of the miRNA, typically a crucial part of miRNA-target recognition (character).
#' @exportClass MicroRNAGene
#' @family Gene classes
setClass("MicroRNAGene",
         contains = "Gene",
         slots = c(miRNAID = "character", 
                   rnaSequence = "character", 
                   seedSequence = "character")                             
)

#' Class "SmallNuclearRNAGene"
#' @slot snRNAID A unique identifier for the small nuclear RNA (snRNA) product (character).
#' @slot rnaSequence The nucleotide sequence of the RNA encoded by the snRNA gene (character).
#' @exportClass SmallNuclearRNAGene
#' @family Gene classes
setClass("SmallNuclearRNAGene",
         contains = "Gene",
         slots = c(snRNAID = "character",
                   rnaSequence = "character")
)

#' Class "RibosomalRNAGene"
#' @slot rRNAID A unique identifier for the ribosomal RNA (rRNA) product (character).
#' @slot rnaSequence The nucleotide sequence of the RNA encoded by the rRNA gene (character).
#' @exportClass RibosomalRNAGene
#' @family Gene classes
setClass("RibosomalRNAGene",
         contains = "Gene",
         slots = c(rRNAID = "character",
                   rnaSequence = "character")
)

#' Class "TransferRNAGene"
#' @slot tRNAID A unique identifier for the transfer RNA (tRNA) product (character).
#' @slot aminoAcid The amino acid carried by the tRNA (character).
#' @slot rnaSequence The nucleotide sequence of the RNA encoded by the tRNA gene (character).
#' @exportClass TransferRNAGene
#' @family Gene classes
setClass("TransferRNAGene",
         contains = "Gene",
         slots = c(tRNAID = "character",
                   aminoAcid = "character",
                   rnaSequence = "character")
)

#' Class "Pseudogene"
#' @slot parentGeneID The unique identifier for the parent gene from which the pseudogene originated (character).
#' @slot functionalStatus The functional status of the pseudogene, indicating whether it is transcribed, translated, or non-functional (character).
#' @exportClass Pseudogene
#' @family Gene classes
setClass("Pseudogene",
         contains = "Gene",
         slots = c(parentGeneID = "character",
                   functionalStatus = "character")
)

#' Class "EnhancerRNAGene"
#' @slot targetGeneID The unique identifier for the gene targeted by the enhancer RNA (eRNA) (character).
#' @slot rnaSequence The nucleotide sequence of the RNA encoded by the enhancer RNA gene (character).
#' @exportClass EnhancerRNAGene
#' @family Gene classes
setClass("EnhancerRNAGene",
         contains = "Gene",
         slots = c(targetGeneID = "character",
                   rnaSequence = "character")
)

#### constructor functions ####
#' Create a Gene Object
#'
#' This function creates an object of class `Gene`, representing a generic gene with its attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param symbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param description A brief description of the gene (character).
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object indicating chromosome, start, end, and strand.
#' @param product A list containing gene products, such as protein sequences or RNA sequences.
#' @return An object of class `Gene`.
#' @export
Gene <- function(id, symbol, name, description, structure, product) {
  new("Gene", id = id, symbol = symbol, name = name, description = description, structure = structure, product = product)
}

#' Create a ProteinCodingGene Object
#'
#' This function creates an object of class `ProteinCodingGene`, representing a protein-coding gene with additional attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param proteinID A unique identifier for the protein product encoded by the gene (character).
#' @param proteinSequence The amino acid sequence of the protein encoded by the gene (character).
#' @param exonCount The number of exons in the gene (integer).
#' @param description A brief description of the gene (character). Defaults to NA.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object. Defaults to an empty `GRanges` object.
#' @return An object of class `ProteinCodingGene`.
#' @export
ProteinCodingGene <- function(id, hugoSymbol, name, proteinID, proteinSequence, exonCount, description = NA_character_, structure = GRanges()) {
  new("ProteinCodingGene", id = id, hugoSymbol = hugoSymbol, name = name, proteinID = proteinID, proteinSequence = proteinSequence, exonCount = exonCount, description = description, structure = structure)
}

#' Create a LongNonCodingRNAGene Object
#'
#' This function creates an object of class `LongNonCodingRNAGene`, representing a long non-coding RNA gene with additional attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param lncRNAID A unique identifier for the long non-coding RNA (lncRNA) product (character).
#' @param rnaSequence The nucleotide sequence of the RNA encoded by the lncRNA gene (character).
#' @param description A brief description of the gene (character). Defaults to NA.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object. Defaults to an empty `GRanges` object.
#' @return An object of class `LongNonCodingRNAGene`.
#' @export
LongNonCodingRNAGene <- function(id, hugoSymbol, name, lncRNAID, rnaSequence, description = NA_character_, structure = GRanges()) {
  new("LongNonCodingRNAGene", id = id, hugoSymbol = hugoSymbol, name = name, lncRNAID = lncRNAID, rnaSequence = rnaSequence, description = description, structure = structure)
}

#' Create a MicroRNAGene Object
#'
#' This function creates an object of class `MicroRNAGene`, representing a microRNA gene with additional attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param miRNAID A unique identifier for the microRNA (miRNA) product (character).
#' @param seedSequence The seed sequence of the miRNA, typically a crucial part of miRNA-target recognition (character).
#' @param description A brief description of the gene (character). Defaults to NA.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object. Defaults to an empty `GRanges` object.
#' @return An object of class `MicroRNAGene`.
#' @export
MicroRNAGene <- function(id, hugoSymbol, name, miRNAID, seedSequence, description = NA_character_, structure = GRanges()) {
  new("MicroRNAGene", id = id, hugoSymbol = hugoSymbol, name = name, miRNAID = miRNAID, seedSequence = seedSequence, description = description, structure = structure)
}

#' Create a SmallNuclearRNAGene Object
#'
#' This function creates an object of class `SmallNuclearRNAGene`, representing a small nuclear RNA gene with additional attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param snRNAID A unique identifier for the small nuclear RNA (snRNA) product (character).
#' @param rnaSequence The nucleotide sequence of the RNA encoded by the snRNA gene (character).
#' @param description A brief description of the gene (character). Defaults to NA.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object. Defaults to an empty `GRanges` object.
#' @return An object of class `SmallNuclearRNAGene`.
#' @export
SmallNuclearRNAGene <- function(id, hugoSymbol, name, snRNAID, rnaSequence, description = NA_character_, structure = GRanges()) {
  new("SmallNuclearRNAGene", id = id, hugoSymbol = hugoSymbol, name = name, snRNAID = snRNAID, rnaSequence = rnaSequence, description = description, structure = structure)
}

#' Create a RibosomalRNAGene Object
#'
#' This function creates an object of class `RibosomalRNAGene`, representing a ribosomal RNA gene with additional attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param rRNAID A unique identifier for the ribosomal RNA (rRNA) product (character).
#' @param rnaSequence The nucleotide sequence of the RNA encoded by the rRNA gene (character).
#' @param description A brief description of the gene (character). Defaults to NA.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object. Defaults to an empty `GRanges` object.
#' @return An object of class `RibosomalRNAGene`.
#' @export
RibosomalRNAGene <- function(id, hugoSymbol, name, rRNAID, rnaSequence, description = NA_character_, structure = GRanges()) {
  new("RibosomalRNAGene", id = id, hugoSymbol = hugoSymbol, name = name, rRNAID = rRNAID, rnaSequence = rnaSequence, description = description, structure = structure)
}

#' Create a TransferRNAGene Object
#'
#' This function creates an object of class `TransferRNAGene`, representing a transfer RNA (tRNA) gene with its attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param tRNAID A unique identifier for the transfer RNA (tRNA) product (character).
#' @param aminoAcid The amino acid carried by the tRNA (character).
#' @param rnaSequence The nucleotide sequence of the RNA encoded by the tRNA gene (character).
#' @param description A brief description of the gene (character). Defaults to `NA`.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object indicating chromosome, start, end, and strand. Defaults to an empty `GRanges` object.
#' @return An object of class `TransferRNAGene`.
#' @export
TransferRNAGene <- function(id, hugoSymbol, name, tRNAID, aminoAcid, rnaSequence, description = NA_character_, structure = GRanges()) {
  new("TransferRNAGene", id = id, hugoSymbol = hugoSymbol, name = name, tRNAID = tRNAID, aminoAcid = aminoAcid, rnaSequence = rnaSequence, description = description, structure = structure)
}

#' Create a Pseudogene Object
#'
#' This function creates an object of class `Pseudogene`, representing a pseudogene with its attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param parentGeneID The unique identifier for the parent gene from which the pseudogene originated (character).
#' @param functionalStatus The functional status of the pseudogene, indicating whether it is transcribed, translated, or non-functional (character).
#' @param description A brief description of the gene (character). Defaults to `NA`.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object indicating chromosome, start, end, and strand. Defaults to an empty `GRanges` object.
#' @return An object of class `Pseudogene`.
#' @export
Pseudogene <- function(id, hugoSymbol, name, parentGeneID, functionalStatus, description = NA_character_, structure = GRanges()) {
  new("Pseudogene", id = id, hugoSymbol = hugoSymbol, name = name, parentGeneID = parentGeneID, functionalStatus = functionalStatus, description = description, structure = structure)
}

#' Create an EnhancerRNAGene Object
#'
#' This function creates an object of class `EnhancerRNAGene`, representing an enhancer RNA (eRNA) gene with its attributes.
#'
#' @param id A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).
#' @param hugoSymbol The HUGO (Human Genome Organisation) symbol for the gene (character).
#' @param name The full name of the gene (character).
#' @param targetGeneID The unique identifier for the gene targeted by the enhancer RNA (eRNA) (character).
#' @param rnaSequence The nucleotide sequence of the RNA encoded by the eRNA gene (character).
#' @param description A brief description of the gene (character). Defaults to `NA`.
#' @param structure The genomic structure of the gene, typically represented as a `GRanges` object indicating chromosome, start, end, and strand. Defaults to an empty `GRanges` object.
#' @return An object of class `EnhancerRNAGene`.
#' @export
EnhancerRNAGene <- function(id, hugoSymbol, name, targetGeneID, rnaSequence, description = NA_character_, structure = GRanges()) {
  new("EnhancerRNAGene", id = id, hugoSymbol = hugoSymbol, name = name, targetGeneID = targetGeneID, rnaSequence = rnaSequence, description = description, structure = structure)
}





#### accessor function ####
#getter 
#' Get the HUGO Symbol of a Gene Object
#' This function retrieves the HUGO Symbol from an object of class `Gene`.
#' @param gene An object of class `Gene`.
#' @return The HUGO Symbol associated with the gene (character).
#' @export
#' @docType methods
#' @rdname hugoSymbol
setGeneric("hugoSymbol", function(gene) standardGeneric("hugoSymbol"))

#' @rdname hugoSymbol
#' @export
#' @aliases hugoSymbol,Gene-method
setMethod("hugoSymbol", "Gene", function(gene) {
  return(gene@hugoSymbol)
})

#setter
#' Set the HUGO Symbol of a Gene Object
#' This function sets the HUGO Symbol for an object of class `Gene`.
#' @param gene An object of class `Gene`.
#' @param value The HUGO Symbol to assign to the gene (character).
#' @return The updated gene object with the new HUGO Symbol.
#' @export
#' @docType methods
#' @rdname hugoSymbol
setGeneric("hugoSymbol<-", function(gene, value) standardGeneric("hugoSymbol<-"))

#' @rdname hugoSymbol
#' @export
#' @aliases hugoSymbol<-,Gene-method
setMethod("hugoSymbol<-", "Gene", function(gene, value) {
  gene@hugoSymbol <- value
  return(gene)
})



#### lengthProduct #### 
#' Calculate the Length of the Gene Product
#'
#' This generic function calculates the length of the gene product (e.g., protein or RNA sequence) for a given gene object.
#'
#' @param gene An object of a subclass of `Gene` that represents a specific gene type.
#' @return The length of the gene product (integer).
#' @export
setGeneric("lengthProduct", function(gene) standardGeneric("lengthProduct"))


#' @describeIn lengthProduct Calculate the length of the protein sequence for a ProteinCodingGene.
#' @export
setMethod("lengthProduct", "ProteinCodingGene", function(gene) {
  nchar(gene@proteinSequence)
})


#' @describeIn lengthProduct Calculate the length of the RNA sequence for a LongNonCodingRNAGene.
#' @export
setMethod("lengthProduct", "LongNonCodingRNAGene", function(gene) {
  nchar(gene@rnaSequence)
})


#' @describeIn lengthProduct Calculate the length of the RNA sequence for a MicroRNAGene.
#' @export
setMethod("lengthProduct", "MicroRNAGene", function(gene) {
  nchar(gene@rnaSequence)
})


#' @describeIn lengthProduct Calculate the length of the RNA sequence for a SmallNuclearRNAGene.
#' @export
setMethod("lengthProduct", "SmallNuclearRNAGene", function(gene) {
  nchar(gene@rnaSequence)
})


#' @describeIn lengthProduct Calculate the length of the RNA sequence for a RibosomalRNAGene.
#' @export
setMethod("lengthProduct", "RibosomalRNAGene", function(gene) {
  nchar(gene@rnaSequence)
})


#' @describeIn lengthProduct Calculate the length of the RNA sequence for a TransferRNAGene.
#' @export
setMethod("lengthProduct", "TransferRNAGene", function(gene) {
  nchar(gene@rnaSequence)
})


#' @describeIn lengthProduct Calculate the length of the RNA sequence for an EnhancerRNAGene.
#' @export
setMethod("lengthProduct", "EnhancerRNAGene", function(gene) {
  nchar(gene@rnaSequence)
})








