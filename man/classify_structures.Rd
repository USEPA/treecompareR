% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_classyfireR.R
\name{classify_structures}
\alias{classify_structures}
\title{Query ClassyFire by structure}
\usage{
classify_structures(input = NULL, tax_level_labels = chemont_tax_levels, ...)
}
\arguments{
\item{input}{A character vector of structural identifiers: SMILES strings or
InChi strings. May optionally be named. If so, the names will be returned as a
column named  \code{identifier} in the output data.frame. If not named, the
vector itself will be returned as a column named  \code{identifier} in the
output data.frame}

\item{tax_level_labels}{By default, the list of taxonomy levels for
ClassyFire: \code{kingdom, superclass, class, subclass, level5, ...
  level11}.}

\item{...}{Other arguments as for \code{\link{query_classyfire}}.}
}
\value{
A data frame with ClassyFire classifications for each input
structural identifier. Will contain columns \code{identifier},
\code{smiles}, \code{inchikey}, one column for each taxonomy level (defined
in argument \code{tax_level_labels}), and \code{report}. The final
\code{report} column explains what happened if a classification could not
be obtained.
}
\description{
This function takes a vector of structural identifiers (SMILES strings or
InChi strings) and queries the ClassyFire API to get classifications for each
one.
}
\seealso{
\code{\link{classify_inchikeys}}
}
