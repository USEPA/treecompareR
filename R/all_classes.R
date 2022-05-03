
#' HazardComparisonDashboard class
#'
#' @description An S4 class to store the results from the Hazard Comparison Dashboard API
#'
#' @slot meta A list of the DTXSID, CASRN, name, SMILES, and InChIKey
#'
#' @export

setClass(
  Class = 'HazardComparisonDashboard',
  representation = representation(
    meta = 'list'
  )
)
