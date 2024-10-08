% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_methods.R
\name{general_JiangConrath_similarity}
\alias{general_JiangConrath_similarity}
\title{Jiang and Conrath similarity}
\usage{
general_JiangConrath_similarity(
  tree,
  label_A = NULL,
  label_B = NULL,
  node_A = NULL,
  node_B = NULL
)
}
\arguments{
\item{tree}{A phylo object representing a rooted tree, with an information
content attribute IC.}

\item{label_A}{The first node label.}

\item{label_B}{The second node label.}

\item{node_A}{Alternate parameter, the first node number.}

\item{node_B}{Alternate parameter, the second node number.}
}
\value{
The Jiang and Conrath similarity in the given tree of the pair of
nodes.
}
\description{
This determines the Jiang and Conrath similarity for two input nodes in a
given tree. This uses the formulation as described in
\href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
}
\examples{

tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
}
\references{
\insertRef{seco2004intrinsic}{treecompareR}

\insertRef{jiang1997semantic}{treecompareR}
}
