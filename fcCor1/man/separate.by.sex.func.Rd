% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FC_corr_Data.formatting&analysis.functions.R
\name{separate.by.sex.func}
\alias{separate.by.sex.func}
\title{The Function seperate.by.sex.func uses finds indices corresponding to Male and Female data and then seperates annotation,
expression, and phenotype data by sex}
\usage{
separate.by.sex.func(samp.annotation = NULL, expr = NULL,
  pheno = NULL)
}
\arguments{
\item{samp.annotation}{data frame containing the strain and sex for every mouse number}

\item{expr}{data frame containing gene expression data for each mouse number for each probe number}

\item{pheno}{data frame containing the phenotypic data for all recorded phenotypes across all mouse numbers}
}
\value{
a list containing annotation, expression and phenotype data seperated by sex
}
\description{
The Function seperate.by.sex.func uses finds indices corresponding to Male and Female data and then seperates annotation,
expression, and phenotype data by sex
}
