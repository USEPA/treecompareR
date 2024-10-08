% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_methods.R
\name{resnik_similarity}
\alias{resnik_similarity}
\title{Resnik Similarity}
\usage{
resnik_similarity(
  tree = NULL,
  label_1 = NULL,
  label_2 = NULL,
  node_1 = NULL,
  node_2 = NULL
)
}
\arguments{
\item{tree}{The underlying tree being examined.}

\item{label_1}{A label for a node within the \code{tree} object.}

\item{label_2}{A label for a node within the \code{tree} object.}

\item{node_1}{An alternate parameter for a node within the \code{tree} object.}

\item{node_2}{An alternate parameter for a node within the \code{tree} object.}
}
\value{
The Resnik similarity value between the two input nodes.
}
\description{
This function takes in a tree and two input nodes (either labels or node
numbers) and returns the Resnik similarity values of the nodes based on the
tree structure. The function uses the formulation as described in
\href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
Furthermore, this function is a wrapper for a RCPP function. For rapid
calculation of similarity values of several pairs of nodes, consider using
\code{\link{similarity_matrix}}.
}
\references{
\insertRef{lin1998information}{treecompareR}

\insertRef{resnik1995using}{treecompareR}
}
\seealso{
\code{\link{jaccard_similarity}}, \code{\link{lin_similarity}},
\code{\link{jiang_conrath_similarity}}, \code{\link{similarity_matrix}}
}
