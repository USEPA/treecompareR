#' ctxR chemical identifiers
#'
#' @param input_list A list or vector of CASRNS or DTXSIDS
#' @param verbose A Boolean indicating whether some 'progress report' is given.
#'
#' @return A data.table with columns `INPUT`, `DTXSID`, `CASRN`,
#' `PREFFERED_NAME`, `SMILES`, `INCHIKEY` and each row corresponding to a
#' searched term.
#' @export
#' @import data.table
#' @examplesIf FALSE
#' # Grab identifiers
#' input_list <- c('DTXSID7020182', '67-64-1')
#' identifiers <- chemical_identifiers(input_list = input_list)
chemical_identifiers <- function(input_list,
                                 verbose = FALSE){

  n <- length(input_list)
  INPUT <- NULL
  DTXSID <- NULL
  dtxsid <- NULL
  CASRN <- NULL
  casrn <- NULL
  PREFERRED_NAME <- NULL
  preferredName <- NULL
  SMILES <- NULL
  smiles <- NULL
  INCHIKEY <- NULL
  inchikey <- NULL
  searchValue <- NULL
  . <- NULL



  new_table <- data.table::data.table('INPUT' = character(),
                                      'DTXSID' = character(),
                                      'CASRN' = character(),
                                      'PREFERRED_NAME' = character(),
                                      'SMILES' = character(),
                                      'INCHIKEY' = character())



  pages <- ceiling(n/200)

  for (i in 1:pages){
    start_index <- 1 + (i-1)*200
    end_index <- min(n, 200 + (i-1)*200)

    temp_chemicals <- ctxR::chemical_equal_batch(word_list = input_list[start_index:end_index])
    temp_chemicals <- data.table::data.table(temp_chemicals)
    temp_chemicals[, INPUT := searchValue]
    temp_chemicals[, DTXSID := dtxsid]
    temp_chemicals[, CASRN := casrn]
    temp_chemicals[, PREFERRED_NAME := preferredName]
    temp_chemicals[, SMILES := smiles]

    temp_inchikeys <- ctxR::get_chemical_details_batch(DTXSID = temp_chemicals$DTXSID,
                                                       Projection = 'chemicalidentifier')
    temp_inchikeys <- data.table::data.table(temp_inchikeys)
    temp_inchikeys[, INCHIKEY := inchikey]
    temp_inchikeys[, DTXSID := dtxsid]


    chemical <- data.table::merge.data.table(x = temp_chemicals,
                                             y = temp_inchikeys[, .(INCHIKEY,
                                                                    DTXSID)],
                                             by.x = 'DTXSID',
                                             by.y = 'DTXSID',
                                             all = TRUE)

    chemical <- chemical[, .(INPUT, DTXSID, CASRN, PREFERRED_NAME, SMILES, INCHIKEY)]
    new_table <- data.table::rbindlist(list(new_table, chemical))

    if (verbose){
      print(paste('Finished batch', i, 'of', pages))
    }
  }

  return(new_table)
}
