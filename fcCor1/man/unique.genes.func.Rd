% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FC_corr_Data.formatting&analysis.functions.R
\name{unique.genes.func}
\alias{unique.genes.func}
\title{The function identofies probes that code for identical genes and takes the mean gene expression value across probes}
\usage{

  \method{unique}{genes.func}(probe.vs.gene = expr.formatted.gene.annotation,
  expr.mean = NULL)
}
\arguments{
\item{probe.vs.gene}{two column data frame containing the gene name assigned to each probe number}

\item{expr.mean}{data frame containing gene expression data for each strain of interests for a specific sex}
}
\value{
a newly formatted data frame containing gene xpression data for each strain of interest

ind.test.4 with only probes for a unique gene (no duplicate genes)
}
\description{
The first probe in the list is kept an contains reformatted, mean gene expression values
}
