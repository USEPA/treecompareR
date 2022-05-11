#' Sample US EPA CompTox Dashboard results for a data set
#'
#' A data set containing a list of 726 chemicals with related information
#'
#' @format A csv representing 726 rows and 18 columns
#' \describe{
#' \item{ï..DTXSID}
#' \item{PREFERRED_NAME}
#' \item{CASRN}
#' \item{INCHIKEY}
#' \item{IUPAC.NAME}
#' \item{SMILES}
#' \item{INCHI.STRING}
#' \item{MOLECULAR.FORMULA}
#' \item{AVERAGE.MASS}
#' \item{MONOISOTOPIC.MASS}
#' \item{X..of.SOURCES}
#' \item{X..OF.PUBMED.ARTICLES}
#' \item{PUBCHEM.DATA.SOURCES}
#' \item{CPDAT.COUNT}
#' \item{QC.Level}
#' \item{X..ToxCast.Active}
#' \item{X..ToxCast.Active.1}
#' \item{Total.Assays"}
#' }
'Chemical List BIOSOLIDS-2022-05-10.csv'




#' Parent and child data for chemont tree
#'
#' A data set containing the information required to build the chemont tree.
#'
#' @format A data frame with 4825 rows and 3 columns:
#' \describe{
#' \item{Name}{Name of the chemont label}
#' \item{ID}{ID of chemont label}
#' \item{Parent_ID}{ID of parent chemont label}
#' }
'chemont_parent_child'


#' Jaccard similarity matrix for chemont tree
#'
#' A matrix containing the Jaccard similarity values for all pairs of nodes in
#' the chemont tree.
#'
#' @format A symmetric matrix with 4824 rows and 4824 columns:
#' \describe{
#' \item{rownames}{tip and node labels of chemont tree}
#' \item{colnames}{tip and node labels of chemont tree}
#' }
"chemont_jaccard"

#' Resnik similarity matrix for chemont tree
#'
#' A matrix containing the Resnik similarity values for all pairs of nodes in
#' the chemont tree. This is derived using information content as defined in "An
#' intrinsic information content metric for semantic similarity in WordNet",
#' Seco et al 2004.
#'
#' @format A symmetric matrix with 4824 rows and 4824 columns:
#' \describe{
#' \item{rownames}{tip and node labels of chemont tree}
#' \item{colnames}{tip and node labels of chemont tree}
#' }
"chemont_resnik_IC_SVH"

#' Lin similarity matrix for chemont tree
#'
#' A matrix containing the Lin similarity values for all pairs of nodes in
#' the chemont tree. This is derived using information content as defined in "An
#' intrinsic information content metric for semantic similarity in WordNet",
#' Seco et al 2004.
#'
#' @format A symmetric matrix with 4824 rows and 4824 columns:
#' \describe{
#' \item{rownames}{tip and node labels of chemont tree}
#' \item{colnames}{tip and node labels of chemont tree}
#' }
"chemont_lin_IC_SVH"

#' Jiang and Conrath similarity matrix for chemont tree
#'
#' A matrix containing the Jiang and Conrath similarity values for all pairs of nodes in
#' the chemont tree. This is derived using information content as defined in "An
#' intrinsic information content metric for semantic similarity in WordNet",
#' Seco et al 2004.
#'
#' @format A symmetric matrix with 4824 rows and 4824 columns:
#' \describe{
#' \item{rownames}{tip and node labels of chemont tree}
#' \item{colnames}{tip and node labels of chemont tree}
#' }
"chemont_jiangconrath_IC_SVH"
