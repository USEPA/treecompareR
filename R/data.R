#' BIOSOLIDS2021 chemical list from the US EPA CompTox Chemicals Dashboard
#'
#' A data set containing a list of 726 chemicals with related information
#'
#' \url{https://comptox.epa.gov/dashboard/chemical-lists/BIOSOLIDS2021}
#'
#' @format A csv representing 726 rows and 18 columns
#' \describe{
#' \item{DTXSID}{A chemical identifier}
#' \item{PREFERRED_NAME}{A character string of the name}
#' \item{CASRN}{A character string of the CASRN identifier}
#' \item{INCHIKEY}{A character of the InChIKey identifier}
#' \item{SMILES}{A character string of the SMILES identifier}
#' }
#' @seealso \code{\link{BIOSOLIDS2021_class}}
'BIOSOLIDS2021'

#' USGSWATER chemical list from the US EPA CompTox Chemicals Dashboard
#'
#' A data set containing a list of 707 chemicals with related information
#'
#' \url{https://comptox.epa.gov/dashboard/chemical-lists/USGSWATER}
#'
#' @format A csv representing 707 rows and 18 columns
#' \describe{
#' \item{DTXSID}{A chemical identifier}
#' \item{PREFERRED_NAME}{A character string of the name}
#' \item{CASRN}{A character string of the CASRN identifier}
#' \item{INCHIKEY}{A character of the InChIKey identifier}
#' \item{SMILES}{A character string of the SMILES identifier}
#' }
#' @seealso \code{\link{USGSWATER_class}}
'USGSWATER'

#' BIOSOLIDS2021 chemical list with ClassyFire classifications
#'
#' A data set containing a list of 726 chemicals with related information
#'
#' Original data source:
#' \url{https://comptox.epa.gov/dashboard/chemical-lists/BIOSOLIDS2021}.
#' ClassyFire classifications obtained via ClassyFire API.
#'
#' @format A csv representing 726 rows and 18 columns \describe{ \item{DTXSID}{A
#'   chemical identifier} \item{PREFERRED_NAME}{A character string of the name}
#'   \item{CASRN}{A character string of the CASRN identifier} \item{INCHIKEY}{A
#'   character of the InChIKey identifier} \item{SMILES}{A character string of
#'   the SMILES identifier} \item{kingdom}{ClassyFire kingdom label}
#'   \item{superclass}{ClassyFire superclass label} \item{class}{ClassyFire
#'   class label} \item{subclass}{ClassyFire subclass label} \item{level
#'   5}{ClassyFire level 5 label} \item{level 6}{ClassyFire level 6 label}
#'   \item{level 7}{ClassyFire level 7 label} \item{level 8}{ClassyFire level 8
#'   label} \item{level 9}{ClassyFire level 9 label} \item{level 10}{ClassyFire
#'   level 10 label} \item{level 11}{ClassyFire level 11 label} }
#' @seealso \code{\link{BIOSOLIDS2021}}
'BIOSOLIDS2021_class'

#' USGSWATER chemical list with ClassyFire classifications
#'
#' A data set containing a list of 707 chemicals with related information
#'
#' Original data source:
#' \url{https://comptox.epa.gov/dashboard/chemical-lists/USGSWATER}. ClassyFire
#' classifications obtained via ClassyFire API.
#'
#' @format A csv representing 707 rows and 18 columns \describe{ \item{DTXSID}{A
#'   chemical identifier} \item{PREFERRED_NAME}{A character string of the name}
#'   \item{CASRN}{A character string of the CASRN identifier} \item{INCHIKEY}{A
#'   character of the InChIKey identifier} \item{SMILES}{A character string of
#'   the SMILES identifier} }
#' @seealso \code{\link{USGSWATER}}
'USGSWATER_class'

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
