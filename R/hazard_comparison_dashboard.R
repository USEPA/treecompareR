#' Get chemical identifiers
#'
#' This function retrieves relevant chemical information for input chemical
#' identifiers. This uses the
#' \href{https://hazard.sciencedataexperts.com/#/}{Hazard Comparison Dashboard}
#' API to retrieve relevant chemical information.
#'
#' @param input A string representing one of DTXSID, CAS, SMILES, NAME,
#'   InChIKey.
#' @param type A string specifying which type of input format the `input`
#'   parameter takes.
#' @return An object of class 'HazardComparisonDashboard'.
#' @export
#' @import httr
#' @import methods
#'
get_chemical_identifiers <- function(input, type = c('AnyId', 'DTXSID', 'CAS')){
  type <- match.arg(type)
  #print(type)
  entity_url <- 'http://hazard-dev.sciencedataexperts.com/api/resolver/lookup?query='


  entity_query <- paste0(entity_url, input, '&idType=', type, '&fuzzy=Not&mol=false')
  response <- httr::RETRY(
    verb = 'get',
    url = entity_query,
    times = 10,
    terminate_on = c(404),
    quiet = T
  )

  if (response$status_code == 429){
    stop('Request rate limit exceeded!')
  }

  if (response$status_code == 404){
    cat(crayon::red(clisymbols::symbol$cross, input), '\n')
  }

  if (response$status_code == 200){
    content <- httr::content(response)

    if (length(content) == 0){
      cat(crayon::red(clisymbols::symbol$cross, input), '\n')
      return(invisible(NULL))
    } else {
      cat(crayon::green(clisymbols::symbol$tick, input), '\n')
    }



    dtxsid <- content[[1]]$sid
    casrn <- content[[1]]$casrn
    name <- content[[1]]$name
    smiles <- content[[1]]$smiles
    InChIKey <- content[[1]]$inchiKey

    object <- methods::new('HazardComparisonDashboard')

    object@meta <-
      list(
        dtxsid = dtxsid,
        casrn = casrn,
        name = name,
        smiles = smiles,
        inchikey = InChIKey
      )


    return(object)
  }


}

#' Batch chemical identifiers.
#'
#' This function collects chemical identifiers and store them in a data.table.
#' It uses the \code{\link{get_chemical_identifiers}} function to access the
#' \href{https://hazard.sciencedataexperts.com/#/}{Hazard Comparison Dashboard}
#' API to retrieve relevant chemical information.
#'
#' @param input_list A list of input values valid for
#'   `get_chemical_identifiers()`
#' @return A data.table with the input list as a column, and chemical
#'   identifiers as subsequent columns
#' @export
#' @import data.table
#'
#' @seealso \code{\link{get_chemical_identifiers}}
batch_chemical_identifiers <- function(input_list){
  INPUT <- NULL
  new_table <- data.table('INPUT' = input_list,
                          'DTXSID' = character(length(input_list)),
                          'CASRN' = character(length(input_list)),
                          'PREFERRED_NAME' = character(length(input_list)),
                          'SMILES' = character(length(input_list)),
                          'INCHIKEY' = character(length(input_list)))

  list_of_identifiers <- purrr::map(input_list, get_chemical_identifiers)

  for (i in seq_along(list_of_identifiers)){
    if (is.null(list_of_identifiers[[i]])){
      new_row <- rep('', 5)
    } else {
      current_item <- list_of_identifiers[[i]]@meta
      new_row <- list(current_item$dtxsid, current_item$casrn, current_item$name, current_item$smiles, current_item$inchikey)
      null_indices <- which(sapply(new_row, is.null))
      new_row[null_indices] <- ''
    }
    new_table[INPUT == input_list[[i]], c('DTXSID', 'CASRN', 'PREFERRED_NAME', 'SMILES', 'INCHIKEY') := as.list(new_row)]
  }

  return(new_table)
}



#' Convert SMILES to hex
#'
#' This is a helper function for converting SMILES string to API acceptable
#' format. This is used by \code{\link{get_toxprints}} for making API calls to
#' \href{https://hazard.sciencedataexperts.com/#/}{Hazard Comparison Dashboard}
#' to retrieve relevant chemical information.
#'
#' @param smiles A string representing a chemical in SMILES format
#' @return A string converted for use in API calls of the
#'   \href{https://hazard.sciencedataexperts.com/#/}{Hazard Comparison
#'   Dashboard}.
convert_smiles_to_hex <- function(smiles = NULL){
  if (is.null(smiles)){
    smiles <- readline('Input a SMILES string:')
  }
    hex_replacement <- c('[(]' = '%28',
                            '[)]' = '%29',
                            '\\[' = '%5b',
                            '[]]' = '%5d',
                            '[\\.]' = '%2e',
                            '[-]' = '%2d',
                            '[=]' = '%3d',
                            '[#]' = '%23',
                            '[$]' = '%24',
                            '[:]' = '%3a',
                            '[/]' = '%2f',
                            '[\\\\]' = '%5c',
                            '[@]' = '%40')
    new_smiles <- stringr::str_replace_all(smiles, hex_replacement)
    return(new_smiles)
}


#' Get toxprints
#'
#' This is a function for retrieving toxprint data for a given chemical from the
#' \href{https://hazard.sciencedataexperts.com/#/}{Hazard Comparison Dashboard}
#' API to retrieve relevant chemical information.
#'
#' @param smiles A string representing a chemical in SMILES format
#' @return A list of chemical descriptors, including ToxPrints descriptors
#' @export
#' @import httr
#'
get_toxprints <- function(smiles) {
  entity_url <- 'https://hazard.sciencedataexperts.com/api/descriptors?smiles='

  converted_smiles <- convert_smiles_to_hex(smiles)

  entity_query <- paste0(entity_url, converted_smiles, '&type=toxprints&headers=false')
  #print(entity_query)
  response <- httr::RETRY(
    verb = 'get',
    url = entity_query,
    times = 10,
    terminate_on = c(404),
    quiet = T
  )

  #return(response)

  if (response$status_code == 429){
    stop('Request rate limit exceeded!')
  }

  if (response$status_code == 404){
    cat(crayon::red(clisymbols::symbol$cross, smiles), '\n')
  }

  if (response$status_code == 200){
    content <- httr::content(response)

    if (length(content) == 0){
      cat(crayon::red(clisymbols::symbol$cross, smiles), '\n')
      return(invisible(NULL))
    } else {
      cat(crayon::green(clisymbols::symbol$tick, smiles), '\n')
    }

    information <- content$chemicals[[1]]
    return(information)
  }

}
