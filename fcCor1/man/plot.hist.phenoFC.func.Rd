% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FC_corr_graphs.function.R
\name{plot.hist.phenoFC.func}
\alias{plot.hist.phenoFC.func}
\title{Histogram of the strain specific FC of a specific phenotype}
\usage{
plot.hist.phenoFC.func(dPhno.interest = NULL, xlab = "")
}
\arguments{
\item{dPhno.interest}{data frame with all the phenotype FC for each strain of a specific phenotype}

\item{xlab}{x axis name for histogram}

\item{dPhno.new.name}{Name that you want to assign to the reorganized data frame for the pheno of interest}
}
\value{
one histogram and a list called 'resu' with the reorganized data frame 'DPhno.interest'
histogram of the Fc of a desired phenotype
reorganized data frame for a specific phenotype containing pheno value for each strain (NA's excluded)
}
\description{
The function returns a histogram showing the distribution of the strain specific Fold-Change of a specific phenotype
}
