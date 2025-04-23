#'@title get results
#'
#'@description Helper function to get ClassyFire query results
#'
#'@details This is a helper function used by [query_classyfire()]. It generally
#'  should not be called directly by the user.
#'
#'  # Return value details
#'
#'  - `id`: Character: the string giving the ClassyFire query id. This is parsed from function argument `url`. If `url` is not of the expected structure, `id` will be `NA_character_`.
#'  - `url`: Character: The URL of the request, as provided in function argument `url`.
#'  - `query_status`: Character: The HTTP status message of the GET request (if one was available), or the error message if no HTTP status was available.
#'  - `label`: Character: The text label assigned to the query, if any.
#'  - `classification_status`: Character: the status of the ClassyFire classification. One of `'Done'`, `'In Queue'`, `'Processing'`, or `'Failed'`.
#'  - `number_of_elements`: Integer: the number of elements returned. If something went wrong, the placeholder value is `NA_integer_`.
#'  - `number_of_pages`: Integer: the number of pages of results. If something went wrong, the placeholder value is 1.
#'  - `invalid_entities`: A `data.frame` of entities in the query whose structures were invalid and no classification could be made. May be an empty `list` if there were no invalid entities, or if classification failed. If not empty, will be a `data.frame` with variables `structure` (character, containing the input structure) and `report` (character, explaining what went wrong).
#'  - `entities`: A nested `data.frame` with information about entities in the query whose structures were valid and which could be classified.
#'
#'  `entities` should be a `data.frame` with the following variables:
#'
#'  - `identifier` Character: the entity identifier.
#'  - `smiles` Character: the entity SMILES string.
#'  - `inchikey` Character: the entity InChIKey.
#'  - `kingdom` A `data.frame` with the entity's kingdom-level classification.
#'  - `superclass` A `data.frame` with the entity's superclass-level classification.
#'  - `class` A `data.frame` with the entity's class-level classification.
#'  - `subclass` A `data.frame` with the entity's subclass-level classification.
#'  - `intermediate_nodes` A `list` of `data.frame`s, one list element for each entity. Each `data.frame` will contain information about any intermediate-level classifications between subclass-level and direct parent.
#'  - `direct_parent` A `data.frame` with the entity's direct parent classification (i.e. the most specific level of classification available).
#'  - `classification_version` Character: describes the version of ClassyFire.
#'
#'  If something went wrong, the following placeholders are used for `entities`:
#'  all the non-`data.frame` elements will be present but assigned
#'  `character(0)`. `entities$intermediate_nodes` will be present and assigned
#'  an empty list.
#'
#'  The `data.frame`s in `entities$kingdom`, `entities$superclass`,
#'  `entities$class`, `entities$subclass`, and `entities$direct_parent` -- and
#'  each of the `data.frame`s in the list `entities$intermediate_nodes`, if any
#'  -- should have the following variables:
#'
#'  - `name` Character: the name of the classification at this level
#'  - `description` Character: Description of the classification at this level.
#'  - `chemont_id` Character: the ChemOnt node ID of the classification at this level.
#'  - `url` Character: The URL for the ChemOnt node of the classification at this level.
#'
#'  If something went wrong, the placeholder for `entities$kingdom`,
#'  `entities$superclass`, `entities$class`, `entities$subclass`, and
#'  `entities$direct_parent` is a `data.frame` with the above-listed variables,
#'  but zero rows (i.e., each variable was assigned `character(0)`). The
#'  placeholder for `entities$intermediate_nodes` is just an empty list.
#'
#'@param url As for [query_classyfire()]. Character: A URL for a GET query.
#'  Usually, this should be of format
#'  `'http://classyfire.wishartlab.com/queries/NNNNN.json'`, where NNNNN are one
#'  or more digits 0-9 denoting a query ID number.
#'@param label As for [query_classyfire()]. Character: A label for the query.
#'  Default `"query"`.
#'@param retry_get_times As for [query_classyfire()]. Integer: number of times
#'  to retry GET query if HTTP error was returned (except for errors in
#'  `terminate_on`, which will terminate immediately without retrying).
#'@param wait_sec As for [query_classyfire()]. Integer: number of seconds to
#'  wait between retries. Minimum 5.
#'@param terminate_on Integer vector: List of HTTP status codes which will
#'  immediately terminate with no more retries. Default: ` c(400:407, 409:418,
#'  421:428, 431, 451, 500:511)`
#'@return A list containing the parsed JSON results (if query was successful),
#'  or placeholder values for the expected elements of the parsed JSON results.
#'  Will contain the following named elements: `id`, `query_url`, `query_status`, `label`,
#'  `classification_status`, `number_of_elements`, `number_of_pages`,
#'  `invalid_entities`, and `entities`. See Details.
#'@author Caroline Ring, Paul Kruse
#' @examples
#' # No example
#'

get_results <- function(url,
                        label = "query",
                        retry_get_times = 3,
                        wait_sec = 5,
                        terminate_on = c(400:407, #but not 408
                                         409:418,
                                         421:428, #but not 429
                                         431, 451,
                                         500:511),
                        type = "structure" #or "inchikey"
                        ){

  json_parse <- tryCatch(
    {
      #try getting results for query
      resp <- tryCatch(
        httr::RETRY(verb = "GET",
                    url = url,
                    encode = "json",
                    times = retry_get_times,
                    pause_min = wait_sec,
                    terminate_on = terminate_on
        ),
        error = function(err){
          err
        }
      )

      #if httr::RETRY() itself threw an error, then resp will be of class "simpleError"
      #otherwise, resp will be of class "response"
      if("simpleError" %in% class(resp)){
        message(paste("ClassyFire query failed. Error message:",
                      paste0(resp$message, "."),
        ))
        json_res <- NULL
      }else{
        #check for http error
        had_error <- tryCatch(
          #will throw error if resp isn't a reasonable http response
          httr::http_error(resp),
          error = function(err){
            #if httr::http_error itself throws an error,
            #then set had_error to TRUE
            TRUE
          }
        )

        #get http status
        request_status <- tryCatch(
          #will throw error if resp has unknown status code
          #or is generally not something reasonable
          httr::http_status(resp),
          error = function(err){
            err
          }
        )

        if(had_error %in% TRUE){
          #the following will work whether httr::http_status() returned a status or
          #an error:
          message(
            paste("ClassyFire query failed. HTTP error status:\n",
                  request_status$message)
          )
          json_res <- NULL
        }else{ #if request successful
          json_res <- httr::content(x = resp, as = "text", encoding = "UTF-8")

        }
      } #end else (if resp was not "simpleError")


      #json_res should be NULL if anything failed
      #if it's non-NULL, try to parse it
      if(!is.null(json_res)){
        #if json_res is not parseable,
        #then just return the request status and the JSON parse error
        json_parse <- tryCatch({
          tmp <- jsonlite::fromJSON(json_res)
          tmp$query_status <- request_status$message
          tmp
        }
          ,
          error = function(err){
            list("query_status" = paste(request_status$message,
                                  ". JSON parsing failed with error:",
                                  err$message),
                 "classification_status" = "Failed")
          }
        )

      }else{
        #if json_res is NULL, something failed
        #just record the status (error message)
        json_parse <- list("query_status" = request_status$message,
                           "classification_status" = "Failed")
      }
      json_parse
    }, #end main try/catch expression to try
    error = function(err){
      #if all else fails, just return whatever error was thrown
      list("query_status" = paste0("Error:",
                             err$message),
           "classification_status" = "Failed")
    }
  ) #end tryCatch

#initialize a placeholder json_parse in case everything else fails

  #For type = "structure",
#the following items are expected:
# [1] "id"                    "label"                 "classification_status"
# [4] "number_of_elements"    "number_of_pages"       "invalid_entities"
# [7] "entities"
  if(type %in% "structure"){
m <- regexec(text=url, pattern = "(\\d+)\\.json")
query_id <- regmatches(x = url, m = m)[[1]][2]
json_parse_default <- list(
  "id" = query_id,
  "query_url" = url,
  "query_status" = NA_character_,
  "label" = label,
  "classification_status" = "Unknown",
  "number_of_elements" = NA_real_,
  "number_of_pages" = NA_real_,
  "invalid_entities" = list(),
  "entities" = data.frame("identifier" = character(0),
                          "smiles" = character(0),
                          "inchikey" = character(0),
                          "kingdom" = data.frame("name" = character(0),
                                                 "description" = character(0),
                                                 "chemont_id" = character(0),
                                                 "url" = character(0)),
                          "superclass" = data.frame("name" = character(0),
                                                    "description" = character(0),
                                                    "chemont_id" = character(0),
                                                    "url" = character(0)),
                          "class" = data.frame("name" = character(0),
                                               "description" = character(0),
                                               "chemont_id" = character(0),
                                               "url" = character(0)),
                          "subclass" = data.frame("name" = character(0),
                                                  "description" = character(0),
                                                  "chemont_id" = character(0),
                                                  "url" = character(0)),
                          "intermediate_nodes" = list(),
                          "direct_parent" = data.frame("name" = character(0),
                                                       "description" = character(0),
                                                       "chemont_id" = character(0),
                                                       "url" = character(0)),
                          "alternative_parents" = list(),
                          "molecular_framework" = character(0),
                          "substituents" = list(),
                          "description" = character(0),
                          "external_descriptors" = list(),
                          "ancestors" = list(),
                          "predicted_chebi_terms" = list(),
                          "predicted_lipidmaps_terms" = list(),
                          "classification_version" = character(0))
)
  }else if(type %in% "inchikey"){
    #basically the same, except everything in "entities" is moved up one level,
    #and there's no number of elements, number of pages, or invalid entities.
    json_parse_default <- list(
      "query_url" = url,
      "query_status" = NA_character_,
      "classification_status" = "Unknown",
     "identifier" = character(0),
                              "smiles" = character(0),
                              "inchikey" = character(0),
                              "kingdom" = data.frame("name" = character(0),
                                                     "description" = character(0),
                                                     "chemont_id" = character(0),
                                                     "url" = character(0)),
                              "superclass" = data.frame("name" = character(0),
                                                        "description" = character(0),
                                                        "chemont_id" = character(0),
                                                        "url" = character(0)),
                              "class" = data.frame("name" = character(0),
                                                   "description" = character(0),
                                                   "chemont_id" = character(0),
                                                   "url" = character(0)),
                              "subclass" = data.frame("name" = character(0),
                                                      "description" = character(0),
                                                      "chemont_id" = character(0),
                                                      "url" = character(0)),
                              "intermediate_nodes" = list(),
                              "direct_parent" = data.frame("name" = character(0),
                                                           "description" = character(0),
                                                           "chemont_id" = character(0),
                                                           "url" = character(0)),
                              "alternative_parents" = list(),
                              "molecular_framework" = character(0),
                              "substituents" = list(),
                              "description" = character(0),
                              "external_descriptors" = list(),
                              "ancestors" = list(),
                              "predicted_chebi_terms" = list(),
                              "predicted_lipidmaps_terms" = list(),
                              "classification_version" = character(0)
    )
  }else{
  stop("'type' must be either 'structure' or 'inchikey'")
}

#fill in any elements in json_parse_default not in json_parse
#otherwise, keep elements in json_parse
json_parse <- c(
  json_parse,
  json_parse_default[
    setdiff(
      names(json_parse_default),
      names(json_parse)
    )
  ]
)

return(json_parse)
}
