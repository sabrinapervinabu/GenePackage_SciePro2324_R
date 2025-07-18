% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneClasses.R
\name{RibosomalRNAGene}
\alias{RibosomalRNAGene}
\title{Create a RibosomalRNAGene Object}
\usage{
RibosomalRNAGene(
  id,
  hugoSymbol,
  name,
  rRNAID,
  rnaSequence,
  description = NA_character_,
  structure = GRanges()
)
}
\arguments{
\item{id}{A unique identifier for the gene, such as Ensembl or NCBI gene ID (character).}

\item{hugoSymbol}{The HUGO (Human Genome Organisation) symbol for the gene (character).}

\item{name}{The full name of the gene (character).}

\item{rRNAID}{A unique identifier for the ribosomal RNA (rRNA) product (character).}

\item{rnaSequence}{The nucleotide sequence of the RNA encoded by the rRNA gene (character).}

\item{description}{A brief description of the gene (character). Defaults to NA.}

\item{structure}{The genomic structure of the gene, typically represented as a \code{GRanges} object. Defaults to an empty \code{GRanges} object.}
}
\value{
An object of class \code{RibosomalRNAGene}.
}
\description{
This function creates an object of class \code{RibosomalRNAGene}, representing a ribosomal RNA gene with additional attributes.
}
