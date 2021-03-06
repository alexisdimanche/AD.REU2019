% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FC_corr_graphs.function.R
\name{deg.ma.plot}
\alias{deg.ma.plot}
\title{function for plotting DEG analysis results}
\usage{
deg.ma.plot(data = NULL, fname = NULL, fc.cut = 2, fdr.cut = 0.05)
}
\arguments{
\item{data}{is the "output" file from the function deg.analysis containing the FC for every gene and the AVG expression
across MAle and female (and other information of statistical significance)}

\item{fname}{file name of the pdf}

\item{fc.cut}{minimum absolute value for the FC threshold}

\item{fdr.cut}{maximum p value for color annotation}
}
\value{
PDF containing a color annotated scatter plot of the Avg expression vs log2FC for every gene
}
\description{
The function outputs a pdF of a scatter plot of the average expression vs the log2 Fold_change for every gene
Only the genes satisfying threshold values for log2FC and pvalue are color annoted: red = positive FC, blue = negative FC
}
