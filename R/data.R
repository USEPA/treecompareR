#' Sample US EPA CompTox Dashboard results for a data set
#'
#' A data set containing a list of 726 chemicals with related information
#'
#' @format A csv representing 726 rows and 18 columns
#' \describe{
#' \item{DTXSID}{A chemical identifier}
#' \item{PREFERRED_NAME}{A character string of the name}
#' \item{CASRN}{A character string of the CASRN identifier}
#' \item{INCHIKEY}{A character of the InChIKey identifier}
#' \item{IUPAC_NAME}{A character string of the IUPAC identifier}
#' \item{SMILES}{A character string of the SMILES identifier}
#' \item{INCHI_STRING}{A character string of the InCh String identifer}
#' \item{MOLECULAR_FORMULA}{A character string of the molecular formula}
#' \item{AVERAGE_MASS}{The molecular mass}
#' \item{MONOISOTOPIC_MASS}{The monoisotopic mass}
#' \item{NUMBER_OF_SOURCES}{Number of references}
#' \item{NUMBER_OF_PUBMED_ARTICLES}{Number of articles in which the chemical is referenced}
#' \item{PUBCHEM_DATA_SOURCES}{Number of references in the PubChem database}
#' \item{CPDAT_COUNT}{Number of references in the CPDAT}
#' \item{QC_Level}{Character giving quality control level of data}
#' \item{NUMBER_ToxCast_Active}{Number of taxcast active}
#' \item{PERCENT_ToxCast_Active}{Percent of toxcast active}
#' \item{Total_Assays}{Total assays referencing the chemical}
#' }
'chemical_list_biosolids_2022_05_10'




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
