% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetchCbioData.r
\name{makeExpressionTable}
\alias{makeExpressionTable}
\title{Make a table of RNA expression values.}
\usage{
makeExpressionTable(studyList, genes, molecularProfile)
}
\arguments{
\item{studyList}{List of studies}

\item{genes}{Character vector of genes to evaluate}

\item{molecularProfile}{The type of molecular profile to query. For most
patient cohorts it is "rna_seq_v2_mrna", for cell lines "rna_seq_mrna".}
}
\value{
A tibble of RNA expression values, grouped by cancer origin,
one row per sample.
}
\description{
Make a table of RNA expression values.
}
