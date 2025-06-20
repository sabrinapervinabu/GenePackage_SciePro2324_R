% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneClasses.R
\name{lengthProduct}
\alias{lengthProduct}
\alias{lengthProduct,ProteinCodingGene-method}
\alias{lengthProduct,LongNonCodingRNAGene-method}
\alias{lengthProduct,MicroRNAGene-method}
\alias{lengthProduct,SmallNuclearRNAGene-method}
\alias{lengthProduct,RibosomalRNAGene-method}
\alias{lengthProduct,TransferRNAGene-method}
\alias{lengthProduct,EnhancerRNAGene-method}
\title{Calculate the Length of the Gene Product}
\usage{
lengthProduct(gene)

\S4method{lengthProduct}{ProteinCodingGene}(gene)

\S4method{lengthProduct}{LongNonCodingRNAGene}(gene)

\S4method{lengthProduct}{MicroRNAGene}(gene)

\S4method{lengthProduct}{SmallNuclearRNAGene}(gene)

\S4method{lengthProduct}{RibosomalRNAGene}(gene)

\S4method{lengthProduct}{TransferRNAGene}(gene)

\S4method{lengthProduct}{EnhancerRNAGene}(gene)
}
\arguments{
\item{gene}{An object of a subclass of \code{Gene} that represents a specific gene type.}
}
\value{
The length of the gene product (integer).
}
\description{
This generic function calculates the length of the gene product (e.g., protein or RNA sequence) for a given gene object.
}
\section{Methods (by class)}{
\itemize{
\item \code{lengthProduct(ProteinCodingGene)}: Calculate the length of the protein sequence for a ProteinCodingGene.

\item \code{lengthProduct(LongNonCodingRNAGene)}: Calculate the length of the RNA sequence for a LongNonCodingRNAGene.

\item \code{lengthProduct(MicroRNAGene)}: Calculate the length of the RNA sequence for a MicroRNAGene.

\item \code{lengthProduct(SmallNuclearRNAGene)}: Calculate the length of the RNA sequence for a SmallNuclearRNAGene.

\item \code{lengthProduct(RibosomalRNAGene)}: Calculate the length of the RNA sequence for a RibosomalRNAGene.

\item \code{lengthProduct(TransferRNAGene)}: Calculate the length of the RNA sequence for a TransferRNAGene.

\item \code{lengthProduct(EnhancerRNAGene)}: Calculate the length of the RNA sequence for an EnhancerRNAGene.

}}
