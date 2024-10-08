% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_methods.R
\name{get_label_level}
\alias{get_label_level}
\title{Get labels at a specified taxonomy level}
\usage{
get_label_level(data, level_label, tax_level_labels = chemont_tax_levels)
}
\arguments{
\item{data}{A data.frame with data that has been classified by some taxonomy.}

\item{level_label}{A string indicating a taxonomy level of the classified
data.}

\item{tax_level_labels}{A vector of all of the taxonomy levels, in order.
Default is \code{\link{chemont_tax_levels}}, the levels of the ChemOnt
taxonomy.}
}
\value{
The unique labels corresponding to the given level and classified
data.
}
\description{
This function takes in a data.table of chemicals with classification data and
a taxonomy level, and returns all the unique labels for the given taxonomy
level.
}
