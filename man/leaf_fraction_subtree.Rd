% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_diagrams.R
\name{leaf_fraction_subtree}
\alias{leaf_fraction_subtree}
\title{Leaf fraction subtree}
\usage{
leaf_fraction_subtree(
  data_1,
  data_2,
  name_1 = "data_1",
  name_2 = "data_2",
  show_labels = FALSE,
  tax_level_labels = chemont_tax_levels,
  tree = NULL
)
}
\arguments{
\item{data_1}{A data.table of chemicals, classifications, and column
\code{terminal_label}.}

\item{data_2}{A data.table of chemicals, classifications, and column
\code{terminal_label}}

\item{name_1}{Alternate parameter for name of first data.table.}

\item{name_2}{Alternate parameter for name of second data.table.}

\item{show_labels}{Alternate parameter indicating whether to show tip labels.}

\item{tax_level_labels}{An alternate parameter giving the taxonomy levels if
not using ClassyFire taxonomy.}

\item{tree}{An alternate parameter giving a taxonomy if not using ChemOnt.}
}
\value{
A ggtree plot.
}
\description{
This function takes in two data.tables, plots the subtree induced by the
first data.table and colors the tips based on the proportion of chemicals
from the second data.table that make up the chemicals from the first, grouped
by tip (or terminal) label.
}
\seealso{
\code{\link{add_terminal_label}}
}
