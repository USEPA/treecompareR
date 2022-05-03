#' This function retrieves relevant chemical information for input chemical identifiers
#'
#' @param input A string representing one of DTXSID, CAS, SMILES, NAME, InChIKey
#' @param type A string specifying which type of input format the `input` parameter takes.
#' @return An object of class 'HazardComparisonDashboard'.
#' @export
#' @import httr
#' @import crayon
#' @import clisymbols
#' @import methods
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
