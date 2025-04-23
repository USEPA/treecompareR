
#'@title InChIKey classification
#'
#'@description This function uses the ClassyFire API to classify chemicals from
#'  an input data.table using the InChIKey chemical identifier.
#'
#'@details This function queries ClassyFire's lookup table of pre-classified
#'  InChiKeys to get classifications, using the ClassyFire API.
#'
#' \itemize{
#'  \item{identifier: The input InCHiKey that was queried. For example,
#'  "XGQJGMGAMHFMAO-UHFFFAOYSA-N"}
#'  \item{smiles: The corresponding SMILES returned by ClassyFire, if any}
#'  \item{inchikey: The InCHiKey returned by ClassyFire, if any. Will be of
#'  format "InCHiKey=XGQJGMGAMHFMAO-UHFFFAOYSA-N"}
#'  \item{classification_version: The version number returned by ClassyFire.
#'  For example, "2.1"}
#'  \item{level: The names of all levels in this InCHiKey's classification.
#'  Will be elements of \code{tax_level_labels}.}
#'  \item{name: The name or label of each level in this InCHiKey's
#'  classification.}
#'  \item{report: A text string reporting the status of the classification,
#'  "ClassyFire returned a classification" if successful; otherwise the report
#'  returned by ClassyFire, or a report about an internal error.}
#' }
#'
#'  Note that the data.frame is in "long" format, with multiple rows for each
#'  InCHiKey. There is one row for each level of classification for each
#'  InCHiKey. For example, the classification for InCHiKey
#'  "XGQJGMGAMHFMAO-UHFFFAOYSA-N" terminates at level 5, so there will be five
#'  rows for that InCHiKey. The classification for InCHiKey
#'  "PZNXLZZWWBSQQK-UHFFFAOYSA-N" terminates at level 4 (subclass), so there
#'  will be four rows for that InCHiKey.
#'
#'  For use with \code{treecompareR} functions that expect a `data.frame` of
#'  classified entities, this `data.frame` will need to be reshaped into wider
#'  format, with one row for each InCHiKey and one column for each level of
#'  classification. This can be done, e.g., using \code{tidyr::pivot_wider(dat,
#'  names_from = "level", values_from = "name")} (where \code{dat} is the
#'  returned `data.frame`.) However, be on the lookout for pathological cases
#'  where an InCHiKey is listed with two different labels at the same level.
#'  These occur rarely, but they do occur. \code{tidyr::pivot_wider()} will
#'  throw a warning if this happens -- pay attention to it!
#'
#'@param inchikeys Character: A vector of InCHiKeys to be classified
#'@param tax_level_labels By default, the list of taxonomy levels for
#'  ClassyFire: \code{kingdom, superclass, class, subclass, level5, ...
#'  level11}.
#'@param wait_sec Numeric: A parameter controlling how many seconds between
#'  qqueries sent to the ClassyFire API server. Default 5, to respect the limit
#'  of 12 requests per second. Setting it any lower than 5 will result in failed
#'  queries.
#'@return A `data.frame`; see Details.
#'
#'
#' @examples
#' #an example with good, bad, NA, and blank inputs
#' classify_inchikeys(
#' c("SEMRCUIXRUXGJX-UHFFFAOYSA-N",
#' "PLDWAJLZAAHOGG-UHFFFAOYSA-N",
#' "BAD-INCHIKEY-X",
#' NA_character_,
#' "    ")
#' )
#'
#'@export
#'@references \insertRef{djoumbou2016classyfire}{treecompareR}
#'
#'
#'@seealso \code{\link{classify_structures}}
#'

classify_inchikeys <- function(inchikeys,
                               tax_level_labels = chemont_tax_levels,
                               wait_sec = 5){

  #placeholder output for everything
  #do this for all original inputs
  output <- data.frame("identifier" = inchikeys)
  #add placeholder columns for ClassyFire classifications at all levels
  output[tax_level_labels] <- NA_character_
  #add placeholder columns for ClassyFire-derived smiles, inchikey
  output$smiles <- NA_character_
  output$inchikey <- NA_character_


  INCHIKEYS <- unique(inchikeys) #save time by removing duplicates
  INCHIKEYS <- INCHIKEYS[!is.na(INCHIKEYS) & #remove NAs
                           nzchar(trimws(INCHIKEYS)) & #remove blanks
                           is_inchikey(INCHIKEYS) #remove any non-validly-formatted inchikeys
  ]

  #get classification results as a list, one for each INCHIKEY
  entities_list <- lapply(
    INCHIKEYS,
    function(this_inchikey){
      #if inchikey is valid, query ClassyFire API
      Sys.sleep(wait_sec) #pause minimum number of seconds before querying again
      this_output <- query_inchikey(inchikey = this_inchikey,
                                               wait_sec = wait_sec)
      this_output$identifier <- this_inchikey
      return(this_output)
    } #end function to apply to each INCHIKEY
  )

  class_list <- lapply(entities_list,
                       parse_classification,
                       query_type = "inchikey",
                       tax_level_labels = tax_level_labels)

  #now rowbind the list of data.frames to get one big data.frame
  inchi_class <- dplyr::bind_rows(class_list)

  more_output <- sapply(
    c("alternative_parents",
      "molecular_framework",
      "substituents",
      "description",
      "external_descriptors",
      "ancestors",
      "predicted_chebi_terms",
      "predicted_lipidmaps_terms"),
    function(this_item){
      lapply(entities_list,
             function(ii) parse_item(entities = ii,
                                          item = this_item,
                                          query_type = "inchikey")) |>
        dplyr::bind_rows() |>
        dplyr::mutate()
    },
    simplify = FALSE,
    USE.NAMES = TRUE)

  inchi_out <- c(list("classified" = inchi_class),
                 more_output)


  inchi_out <- lapply(inchi_out,
                      function(this_out){
                        #get any input inchikeys that were excluded
                        #or are otherwise missing from the classified output
                        #Get any input structures not handled by ClassyFire
                        #e.g. any duplicates, or batches that failed, etc.
                        #keep their placeholder values.
                        if("identifier" %in% names(this_out)){
                          handled_id <- this_out$identifier
                        }else{
                          handled_id <- character(0)
                        }
                        missing_entities <- output |>
                          dplyr::filter(!(identifier %in% handled_id)) |>
                          dplyr::mutate(entity_type = "not handled by ClassyFire")

                        #Output: bind things handled by ClassyFire and things not handled
                        this_out <- dplyr::bind_rows(missing_entities,
                                                        this_out) |>
                          as.data.frame()  #convert from tibble to data.frame

                      })

  return(inchi_out)
}

#' @title Query ClassyFire by structure
#'
#' @description This function takes a vector of structural identifiers (SMILES
#'   strings or InChi strings) and queries the ClassyFire API to get
#'   classifications for each one.
#'
#' @details
#'
#' This function returns a list of `data.frame`s with the following names.
#'
#' - `classification`: The main ClassyFire classification for the input structures.
#' - `alternative_parents`: Alternative parents for the input structures.
#' - `molecular_framework`: Molecular framework descriptors for the input structures.
#' - `substituents`: Substituents for the input structures.
#' - `description`: Descriptions for the input structures.
#' - `external_descriptors`: External descriptors for the input structures, if any.
#' - `ancestors`: All ancestors for the input structures.
#' - `predicted_chebi_terms`: Predicted ChEBI terms for the input structures.
#' - `predicted_lipidmaps_terms`: Predicted LipidMaps terms for the input structures.
#'
#' @param input A character vector of structural identifiers: SMILES strings or
#'   InChi strings. May optionally be named. If so, the names will be returned
#'   as a column named  \code{identifier} in the output data.frame. If not
#'   named, the vector itself will be returned as a column named
#'   \code{identifier} in the output data.frame
#' @param tax_level_labels By default, the list of taxonomy levels for
#'   ClassyFire: \code{kingdom, superclass, class, subclass, level5, ...
#'   level11}.
#' @param ... Other arguments as for \code{\link{query_structure}}.
#' @return A list of `data.frame`s with ClassyFire classification results for each input
#'   structural identifier.
#'
#'   Will contain variables `identifier`, `structure`,
#'   one variable named for each taxonomy level in `tax_level_labels` containing
#'   the classification at each level), `smiles`,  `inchikey`, and then several
#'   variables with information about the ClassyFire query: `report`,
#'   `classification_version`, `id`, `label`, `classification_status`, `query_url` and
#'   `query_status`. If classification could not be completed for one or more
#'   structures, they will be listed in the `data.frame` with `NA_character_`
#'   for all taxonomy levels. Check variables `report` and `status` for more
#'   information about what may have gone wrong.
#'
#' @examples
#'  classify_structures(
#'input = c("COC1=CC(Br)=CC=C1",
#'          "SCCSCCS",
#'          NA_character_,
#'          "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O",
#'          "  ",
#'          "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
#'          "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
#'          "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
#')
#'
#'
#' @export
#' @seealso \code{\link{classify_inchikeys}}
#'
classify_structures <- function (input = NULL,
                                 tax_level_labels = chemont_tax_levels,
                                 sizelim = 1000,
                                 wait_sec = 5,
                                 ...){

  if(!is.character(input)){
    stop("input must be a character vector of structures (SMILES and/or InChI strings)")
  }

  if(any(is_inchikey(input))){
    warning(paste("One or more input structures appear to be InChIKeys.",
                  "Please use `classify_inchikeys()` to query by InChIKey.")
    )
  }

  #Names of input structures are their identifiers.
  #If no names, use the structures themselves as identifiers.
  if(is.null(names(input))){
    names(input) <- ifelse(!is.na(input),
                           input,
                           "NA")
  }

  #save original input
  input_orig <- input

  #placeholder output: identifier and structure for all original inputs
  output <- data.frame("identifier" = names(input_orig))
  output$structure <- input_orig[output$identifier]

  #filter inputs
  input_df <- data.frame("input1" = input,
                         "identifier1" = names(input)) |>
    #remove duplicates
    dplyr::distinct() |>
    #remove NAs and blanks
    dplyr::filter(!is.na(input1) &
                    nzchar(trimws(input1)))

  #use the filtered inputs
  input <- input_df$input1
  names(input) <- input_df$identifier1


  #Can't send a query with more than sizelim structures.
  #If there are more than sizelim, batch them.

  message(paste("Input contains", length(input),
                "unique, non-NA, non-blank structures."))

  if(length(input) > sizelim){
    message(paste("ClassyFire queries must be less than",
                  sizelim,
                  "structures.",
                  "Dividing input into",
                  ceiling(length(input)/sizelim),
                  "batches before querying."))
    batch_id <- cut(1:length(input),
                    breaks = ceiling(length(input)/sizelim),
                    labels = FALSE)
  }else{
    batch_id <- rep(1, length(input))
  }


  #loop over batches
  batch_output <- sapply(
    unique(batch_id),
    function(this_batch_id){
      #enforce waiting time between batches
      Sys.sleep(wait_sec)
      #select inputs in this batch
      this_batch_input <- input[batch_id %in% this_batch_id]
      #Now, query ClassyFire with these structures.
      json_parse <- do.call(query_structure,
                            args = c(list(input = this_batch_input,
                                          url = NULL),
                                     ...)
      )
      #json_parse will be a list, one element for each page of the results.
      #or at least one element even if there were no results.

      #If classification_status is Failed, In Queue, or Processing, then both
      #entities and invalid_entities will be empty lists.

      #Otherwise, if classification_status is Done, at least one of entities or
      #invalid entities should be non-empty.

      #Now try to parse the ClassyFire output in json_parse into a data.frame.

      #pull list of valid entities (i.e., those with classifications), if any.
      #there will be one for each page of output.
      entities_list <- lapply(json_parse,
                              function(x) x$entities)

      #Parse classified entities. If there were none, parse_classification()
      #will just return a data.frame with zero rows.
      classified <- lapply(entities_list,
                           parse_classification,
                           query_type = "structure") |>
        dplyr::bind_rows() |>
        dplyr::mutate(structure = this_batch_input[identifier])

      #Get the additional pieces of the classification output.
      more_output <- sapply(
        c("alternative_parents",
          "molecular_framework",
          "substituents",
          "description",
          "external_descriptors",
          "ancestors",
          "predicted_chebi_terms",
          "predicted_lipidmaps_terms"),
        function(this_item){
          #parse this item for each page
          this_item_parsed_list <- lapply(entities_list,
                        function(this_page){
                          tryCatch({
                            this_parse <- parse_item(
                              entities = this_page,
                              item = this_item,
                              query_type = "structure",
                              tax_level_labels = tax_level_labels
                            )
                            if(length(this_parse)==0){
                              this_parse <-  data.frame(
                                identifier = this_page$identifier
                              )
                            }

                            this_parse
                          },
                          error = function(err){
                            #if this item couldn't be parsed --
                            #then just return the identifiers
                            data.frame(
                              identifier = this_page$identifier
                            )
                          }
                          ) #end tryCatch
                        }
          ) #end lapply over entities_list

          #rowbind all the pages for this item
          dplyr::bind_rows(this_item_parsed_list) |>
            #add in structures
            dplyr::mutate(structure = this_batch_input[identifier])
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      ) #end sapply over additional items


      #Parse invalid entities, if any.
      #this will be the same for all pages, so just pull from the first.
      #This will be a data.frame or list.
      invalid <- tryCatch(
        json_parse[[1]]$invalid_entities,
        error = function(err){
          return(data.frame())
        }
      )
      if(length(invalid)==0){ #if no invalid entities,
        #return zero-row data.frame with the expected variables
        invalid <- data.frame(identifier = character(0),
                              structure = character(0),
                              report = character(0)
        )
      }else{
        invalid <- tidyr::unnest(invalid, cols = report)

        #Add identifier column:
        #(the names of the inputs that were invalid structures)
        invalid$identifier <- names(this_batch_input)[this_batch_input %in% invalid$structure]
      } # if(length(invalid)>0)

      #Pack all classification outputs into a list

      class_out <- c(list("classification" =  classified),
        more_output)

      #For each piece of classification output,
      #row-bind the invalid entities.
      class_out <- sapply(class_out,
             function(this_item){
               this_new <- dplyr::bind_rows(
                 list(
                   "valid" = this_item,
                   "invalid" = invalid
                 ),
                 .id = "entity_type"
               )

               #add informational columns
               if(nrow(this_new)>0){
                 this_new[c("id",
                               "label",
                               "classification_status",
                               "query_url",
                               "query_status")] <- json_parse[[1]][c("id",
                                                               "label",
                                                               "classification_status",
                                                               "query_url",
                                                               "query_status")]
               }

               return(this_new)
             },
             simplify = FALSE,
             USE.NAMES = TRUE
             )

      return(class_out)
    }, #end function to apply for each batch ID
    simplify = FALSE,
    USE.NAMES = TRUE
  ) #end sapply over batches


  #combine the output from batches
  all_output <- sapply(
    c("classification",
      "alternative_parents",
      "molecular_framework",
      "substituents",
      "description",
      "external_descriptors",
      "ancestors",
      "predicted_chebi_terms",
      "predicted_lipidmaps_terms"),
    function(this_item){
      #extract this element from all items in batch output
      #the result will be a list of data.frames,
      #one for each batch
      this_item_list <- sapply(batch_output,
                               function(x) x[[this_item]],
                               simplify = FALSE,
                               USE.NAMES = TRUE)
      #bind rows to combine all batches
      this_item_df <- dplyr::bind_rows(this_item_list, .id = "batch")

      #Get any input structures not handled by ClassyFire
      #e.g. any duplicates, or batches that failed, etc.
      #keep their placeholder values.
      if("identifier" %in% names(this_item_df)){
        handled_id <- this_item_df$identifier
      }else{
        handled_id <- character(0)
      }
      missing_entities <- output |>
        dplyr::filter(!(identifier %in% handled_id)) |>
      dplyr::mutate(entity_type = "not handled by ClassyFire")

      #Output: bind things handled by ClassyFire and things not handled
      tmp_output <- dplyr::bind_rows(missing_entities,
                                     this_item_df) |>
        as.data.frame()  #convert from tibble to data.frame

    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  return(all_output)
}

#'@title Query ClassyFire by InChIKey
#'
#'@description Query ClassyFire API by InChIKey.
#'
#'@details This is a helper function that used by [classify_inchikeys()]. It
#'  should generally not be called directly by the user.
#'
#'@param inchikey Character: one InChIKey to be queried.
#'@param retry_get_times Integer: The number of times to retry the query before
#'  giving up if it returns an HTTP error. Default 3.
#'@param wait_sec Numeric: The number of seconds to pause between retry
#'  attempts. Default 5, because ClassyFire rate-limits to 12 requests per
#'  minute. If set to less than 5, it will be reset to 5, with a message.
#'@param base_url The base URL for ClassyFire API queries. Default
#'  `"http://classyfire.wishartlab.com/entities"`. If you change this, the query
#'  is highly unlikely to work!
#'
#'@return A list consisting of the parsed JSON for each page of results. A completed
#'  query will have named elements `id`, `query_url`, `query_status`, `label`,
#'  `classification_status`, `number_of_elements`, `number_of_pages`,
#'  `invalid_entities`, and `entities`. A failed query will have named elements
#'  `id`, `query_url`, `query_status`, `label`, `classification_status`, and
#'  `number_of_pages`. `id` gives the numerical query ID (assigned by
#'  ClassyFire). `query_url` gives the queried URL. `label` gives the
#'  user-supplied query label. `query_status` gives the HTTP status of the
#'  request (200 means successful). `classification_status` reports the
#'  ClassyFire status: "Done" means classification was successfully completed;
#'  "In Queue" means the query timed out while it was still queued; "In
#'  Progress" means the query timed out while it was still processing; "Failed"
#'  means the query failed before ClassyFire could report a status (e.g. HTTP
#'  status code other than 200, or ClassyFire returned some content that did not
#'  include a classification status). `number_of_elements` gives the number of
#'  classified entities. `number_of_pages` gives the total number of pages for
#'  this query. `invalid_entities`, if present, is a `data.frame` listing queried
#'  entity identifiers that ClassyFire found invalid. `entities`, if present, is
#'  a nested `data.frame` giving classifications.
#'
#' @examples
#' query_inchikey(inchikey = "PLDWAJLZAAHOGG-UHFFFAOYSA-N")
#'
query_inchikey <- function(inchikey,
                                      retry_get_times = 3,
                                      wait_sec = 5,
                                      terminate_on = c(400:407, #but not 408
                                                       409:418,
                                                       421:428, #but not 429
                                                       431, 451,
                                                       500:511),
                                      base_url = "http://classyfire.wishartlab.com/entities"){

  url <- paste0(base_url,
                "/",
                inchikey,
                ".json")

  json_parse <- get_results(url = url,
                            retry_get_times = retry_get_times,
                            wait_sec = wait_sec,
                            terminate_on = terminate_on,
                            type = "inchikey")

  #assign identifier as inchikey
  json_parse$identifier <- inchikey

  return(json_parse)
}


#' Query ClassyFire by structure
#'
#' This function takes a vector of structural identifiers (SMILES strings or
#' InChi strings) and queries the ClassyFire API to get classifications for each
#' one.
#'
#' See ClassyFire API details at http://classyfire.wishartlab.com/access.
#'
#' @param input A character vector of structural identifiers: SMILES strings or
#'   InChi strings. May optionally be named. Default NULL to use `url` instead.
#'   Either `input` or `url` must be provided.
#' @param url Optional: a string giving an existing ClassyFire query URL to
#'   retry. Default NULL to use `input` instead. Either `input` or `url` must be
#'   provided.
#' @param label An optional text label for the ClassyFire query. Default
#'   \code{"query"}.
#' @param type String giving ClassyFire query type. Default \code{"structure"}.
#' @param retry_get_times Max number of times to retry GET command to ClassyFire
#'   API to get status of query (with wait time in between tries). Default 10.
#' @param wait_sec Minimum number of seconds to wait before retrying any GET
#'   command or query. Default 5 seconds (because ClassyFire rate-limits
#'   requests to 12 per minute). If less than 5 seconds, will be reset to 5
#'   seconds, with a warning. Wait times for multiple retries
#'   use exponential backoff with full jitter (wait time is \code{runif(1,
#'   wait_sec, wait_sec * 2^(attempt_number)}).
#' @param retry_query_times Max number of attempts to retry an "In Queue" or
#'   "Processing" query (with wait time in between tries). Default 10.
#' @param processing_wait_per_input Minimum number of seconds to wait *per input
#'   string* before retrying to retrieve results for a query whose status is
#'   "Processing." Default 0.1. Uses exponential backoff: effective wait time
#'   per input = \code{processing_wait_per_input * 2^(attempt number)}. Total
#'   wait time is either the effective wait time per input times the number of
#'   input structures, or \code{wait_sec} seconds, whichever is greater. (This
#'   means `processing_wait_per_input` effectively can never be less than
#'   `wait_sec/length(input)`, even if you set it to a smaller number.)
#' @param terminate_on Integer vector: List of HTTP status codes which will
#'   immediately terminate with no more retries. Default: ` c(400:407, 409:418,
#'   421:428, 431, 451, 500:511)`
#' @param base_url The base URL for ClassyFire API queries. Default
#'   `"http://classyfire.wishartlab.com/entities"`. If you change this, the
#'   query is highly unlikely to work!
#' @return A list of lists. The outer list has one element for each page of the
#'   JSON output from the ClassyFire query (there is one page for every ten
#'   entities). If the ClassyFire query failed or timed out, the list will have
#'   one element. Each list element is itself a list, consisting of the parsed
#'   JSON result for each page. A completed query will have named elements `id`,
#'   `query_url`, `query_status`, `label`, `classification_status`,
#'   `number_of_elements`, `number_of_pages`, `invalid_entities`, and
#'   `entities`. A failed query will have named elements `id`, `query_url`,
#'   `query_status`, `label`, `classification_status`, and `number_of_pages`.
#'   `id` gives the numerical query ID (assigned by ClassyFire). `query_url`
#'   gives the queried URL. `label` gives the user-supplied query label.
#'   `query_status` gives the HTTP status of the request (200 means successful).
#'   `classification_status` reports the ClassyFire status: "Done" means
#'   classification was successfully completed; "In Queue" means the query timed
#'   out while it was still queued; "In Progress" means the query timed out
#'   while it was still processing; "Failed" means the query failed before
#'   ClassyFire could report a status (e.g. HTTP status code other than 200, or
#'   ClassyFire returned some content that did not include a classification
#'   status). `number_of_elements` gives the number of classified entities.
#'   `number_of_pages` gives the total number of pages for this query.
#'   `invalid_entities`, if present, is a `data.frame` listing queried entity
#'   identifiers that ClassyFire found invalid. `entities`, if present, is a
#'   nested `data.frame` giving classifications.
#'
#'   @examples
#'      query_structure(input = c("COC1=CC(Br)=CC=C1", "SCCSCCS"))
#'
#'
query_structure <- function(input = NULL,
                             url = NULL,
                             label = "query",
                             type = "STRUCTURE",
                             retry_get_times = 10,
                             wait_sec = 5,
                             retry_query_times = 10,
                             processing_wait_per_input = 0.1,
                             terminate_on = c(400:407, #but not 408
                                              409:418,
                                              421:428, #but not 429
                                              431, 451,
                                              500:511),
                             base_url = "http://classyfire.wishartlab.com/queries"){

  #if both input and url are NULL, throw error
  if(is.null(input) &
     is.null(url)){
    stop("Either `input` or `url` must be provided, but both were NULL.")
  }

  if(!is.null(input) &
     !is.null(url)){
    stop("Either `input` or `url` must be provided, but not both. You provided both.")
  }

  if(wait_sec < 5){
    warning(paste("ClassyFire rate-limits requests to 12 per minute.",
                  "You supplied wait_sec =", wait_sec,
                  "; it will be reset to 5 seconds."))
    wait_sec <- 5
  }

  #if input provided, construct new query
  if(!is.null(input)){
    if(is.null(processing_wait_per_input)){
      processing_wait_per_input <- wait_sec/length(input)
    }

    if(is.null(names(input))){
      names(input) <- input
    }

    #create tab-separated input
    query_input <- paste(names(input),
                         input,
                         sep = "\t",
                         collapse = "\n")
    #convert to JSON
    q <- rjson::toJSON(list(label = label,
                            query_input = query_input,
                            query_type = type)
    )
    #construct POST request
    resp <- httr::POST(url = base_url,
                       body = q,
                       httr::content_type_json(),
                       httr::accept_json(),
                       httr::timeout(getOption("timeout"))
    )
    post_cont <- httr::content(resp)
    #get URL for query
    url <- paste0(base_url, "/", post_cont$id, ".json")
  } #end if(!is.null(input))

  #check url format (whether newly constructed or user-provided)
  if(length(url) %in% 1){
    if(is.character(url)){
      url_good <- grepl(x = url,
                        pattern = paste0(base_url, "/", "\\d+", ".json"))
    }else{
      url_good <- FALSE
    }
  }else{
    url_good <- FALSE
  }

  if(!(url_good %in% TRUE)){
    stop(paste0("Error in treecompareR::query_structure():",
                "URL is not of valid format. ",
                "URL is ", url,
                " and expected format is ",
                base_url,
                "/", "NNNNN.json",
                " where NNNNN are one or more digits 0-9 denoting a query ID number."))
  }

  json_parse <- get_results(url = url,
                            label = label,
                            retry_get_times = retry_get_times,
                            wait_sec = wait_sec,
                            terminate_on = terminate_on,
                            type = "structure")

  #retry until Done, Failed, or until max number of retries is reached
  retry_count <- 0
  while(!(json_parse$classification_status %in% c("Done",
                                                  "Failed")) &
        retry_count < retry_query_times){
    if(json_parse$classification_status %in% "In Queue"){
      #wait for query to come out of queue
      #this does not depend on size of input
      #calculate exponential backoff wait time
      wait_time <- runif(1, wait_sec, wait_sec*2^(retry_count))
      message(paste0("Classification status is In Queue",
                     "; waiting ",
                     round(wait_time, digits = 2),
                     " seconds and retrying"))
    }else{ #query is processing
      #processing time *does* depend on size of input
      #calculate exponential backoff wait time
      proc_eff <- runif(1,
                        processing_wait_per_input,
                        processing_wait_per_input * 2^(retry_count)
      )
      #wait_time at least 5 seconds
      wait_time <- ceiling(max(wait_sec,  proc_eff*length(input)))
      message(paste0("Classification status is ",
                     json_parse$classification_status,
                     "; waiting ",
                     round(wait_time, digits = 2),
                     " seconds (greater of ",
                     wait_sec,
                     " seconds, or ",
                     proc_eff,
                     " seconds per input, with ",
                     length(input),
                     " inputs) ",
                     "and retrying"))
    } #end if/else to check classification status
    #wait and retry
    Sys.sleep(wait_time)
    json_parse <- get_results(url = url,
                              label = label,
                              retry_get_times = retry_get_times,
                              wait_sec = wait_sec,
                              terminate_on = terminate_on,
                              type = "structure")
    retry_count <- retry_count + 1
  } #end while loop

  #If classification was done, then results will be paginated
  if(json_parse$classification_status %in% "Done"){

    #Each page = 10 entities
    #Page URLs are constructed like http://classyfire.wishartlab.com/queries/XXX.json?page=1
    #get results for each individual page

    message("Classification done!")
    if(json_parse$number_of_pages > 0){
      #get vector of individual page URLs
      url_pages <- paste0(url,
                          "?page=",
                          seq(from = 1,
                              to = json_parse$number_of_pages,
                              by = 1))
      output <- sapply(url_pages,
                       function(url_pg){
                         #wait between querying pages
                         Sys.sleep(wait_sec)
                         message(paste("Getting page",
                                       url_pg,
                                       "of",
                                       json_parse$number_of_pages,
                                       "..."))
                         #grab the results for this page.
                         get_results(url = url_pg,
                                     label = label,
                                     retry_get_times = retry_get_times,
                                     wait_sec = wait_sec,
                                     terminate_on = terminate_on,
                                     type = "structure")
                       },
                       simplify = FALSE,
                       USE.NAMES = TRUE
      )
    }else{
      #if classification status was Done but there were zero pages, then it
      #means all entities were invalid, in which case it will return not a
      #list of lists of parsed json, but only one list. Pack it into a nested
      #list as expected by classify_structures().
      json_parse <- list(json_parse)
      output <- json_parse
    }

  }else{
    #this block should only occur if it either failed, or status was still "In
    #Queue" or "Processing" when the max number of retries was reached. In which
    #case it will return not a list of lists of parsed json, but only one list.
    #Pack it into a nested list as expected by classify_structures.
    json_parse <- list(json_parse)
    output <- json_parse
  }
  return(output)
}

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



#' @title Parse classified entities
#'
#' @description Parse ClassyFire classification results to a data.frame.
#'
#' @details
#'
#' This is a helper function used by [classify_structures()] and
#' [classify_inchikeys()]. It should generally not be called directly by the
#' user.
#'
#' @param entities JSON element of classified entities as returned by ClassyFire
#'   API query (e.g, from the results of [query_structure()])
#' @param tax_level_labels ChemOnt taxonomy level labels (default
#'   \code{\link{chemont_tax_levels}})
#' @return A data.frame consisting of rows corresponding to each classified
#'   entry from the `entities` input. A data.frame with one row for each entity
#'   identifier, and variables `identifier`, `smiles`, `inchikey`, `kingdom`,
#'   `superclass`, `class`, `subclass`, `level5`, `level6`, ... `level11`,
#'   `classification_version`, and `report`.
#' @examples
#' #get a set of results to parse
#' my_results <- query_structure(input = c("COC1=CC(Br)=CC=C1", "SCCSCCS"))
#' #parse classifications
#' parse_classification(entities = my_results[[1]]$entities)
#'
parse_classification <- function(entities,
                                      query_type,
                                      tax_level_labels = chemont_tax_levels){

  #check to see whether classifications actually exist for these
  if(length(entities)==0){
    #if no classifications, return empty data.frame,
    #but with the expected variable names
    classifications <- data.frame(identifier = character(0),
                                      smiles = character(0),
                                      inchikey = character(0),
                                      classification_version = character(0),
                                      level = character(0),
                                      name = character(0),
                                      description = character(0),
                                      chemont_id = character(0),
                                      url = character(0)
                                      )
  }else{ #if length(entities)>0
    classifications <- lapply(c("kingdom",
             "superclass",
             "class",
             "subclass",
             "intermediate_nodes",
             "direct_parent"),
           function(this_item){
             parse_item(entities = entities,
                             item = this_item,
                             query_type = query_type,
                             tax_level_labels = tax_level_labels)
           }) |>
      dplyr::bind_rows() |>
      #remove any with NA names
      dplyr::filter(!is.na(name))


    #To remove duplicates,
    #get level numbers for all labels in the ChemOnt tree,
    #merge them into classifications,
    #and take unique rows only.
    chemont_tree_df <- get_tree_df(chemont_tree)

    classifications <- dplyr::left_join(classifications,
                                 chemont_tree_df[,
                                                 c("Name",
                                                   "level")],
                                 by = c("name" = "Name")) |>
      #keep only unique rows
      dplyr::distinct() |>
      #arrange by entity and then by level
      dplyr::arrange(identifier, level) |>
      #Change the levels from numbers to taxonomy level names
      #e.g. 1, 2, 3, 4 -> kingdom, superclass, class, subclass
      dplyr::mutate(level = tax_level_labels[level])

  return(classifications)
  }
}

#' @title Parse list item
#'
#' @description Parse an item from a list of classified entities.
#'
#' @details This is a helper function used by [classify_structures()] and
#'   [classify_inchikeys()]. It generally should not be called directly by the
#'   user.
#'
#'   The main classification is parsed by [parse_classification()], but
#'   ClassyFire returns several additional pieces of classification output.
#'   These include:
#'
#' - Levels of classification: `kingdom`, `superclass`, `class`, `subclass`, `intermediate_nodes`, `direct_parent`
#' - Alternative parents (`alternative_parents`)
#' - Molecular framework (`molecular_framework`)
#' - Substituents (`substituents`)
#' - Description (`description`)
#' - External descriptors (`external_descriptors`)
#' - Ancestors (`ancestors`)
#' - Predicted ChEbI terms (`predicted_chebi_terms`)
#' - Predicted LipidMaps terms (`predicted_lipidmaps_terms`)
#'
#'   This function parses any of these additional items into a `data.frame`.
#'
#'   The resulting `data.frame` will have variables `identifier`, `smiles`,
#'   `inchikey`, `classification_version`, and the following additional
#'   variables depending on what item was requested:
#'
#' - `kingdom`, `superclass`, `class`, `subclass`, `intermediate_nodes`, `direct_parent`, and `alternative_parents`: `name`, `description`, `chemont_id`, `url`
#' - `external_descriptors`: If any entities had these, `source` (the external source, e.g. "CHEBI", "KEGG", or "LIPID MAPS"); `source_id` (the ID of the descriptor in the external source, e.g. "CHEBI:15756" from ChEBI, or "LMFA01010001" from LipidMaps); and `annotations`, the annotation of the external descriptor (e.g. "long-chain fatty acid" or "Straight chain fatty acids")
#' - everything else: a single additional variable named after the requested item (i.e., `molecular_framework`, `substituents`, etc.)
#'
#'
#' @param entities A `data.frame` of classified entities from one page of
#'   ClassyFire results, e.g, the `entities` element of one element of the list
#'   output by [query_structure()].
#' @param item Character: The name of an item in `entities`. One of `'alternative_parents'`,
#'   `'molecular_framework'`, `'substituents'`, `'description'`,
#'   `'external_descriptors'`, `'ancestors'`, `'predicted_chebi_terms'`,
#'   `'predicted_lipidmaps_terms'`. See Details.
#'   These are all the named items currently in `entities` as returned by
#'   ClassyFire, so providing any other string here will result in an empty
#'   `data.frame` being returned.
#' @return A `data.frame` containing the information from the specified item.
#'   Contains variables `identifier`, `smiles`, `inchikey`,
#'   `classification_version`, and other variables depending on what item was
#'   specified. If `item = 'alternative_parents'`, additional variables are
#'   `name`, `description`, `chemont_id`, `url`. For any of the other possible
#'   values of `item` listed above, there will be
#'   one additional variable named for `item` (e.g., if `item =
#'   'molecular_framework'`, the additional variable will be named
#'   `molecular_framework`). If
#'   `item` is some string other than those listed above, i.e., not matching any
#'   of the named items in `entities`, then an empty `data.frame` will be
#'   returned (zero rows and no variables).
#' @examples
#' #get a set of results to parse
#' my_results <- query_structure(input = c("COC1=CC(Br)=CC=C1", "SCCSCCS"))
#' #get alternative parents
#' parse_item(entities = my_results[[1]]$entities, item = "alternative_parents")
#' #get molecular framework
#' parse_item(entities = my_results[[1]]$entities, item = "molecular_framework")
#' @author Caroline Ring
parse_item <- function(entities,
                            item,
                            query_type,
                            tax_level_labels = chemont_tax_levels){

  if(length(item) > 1){
    stop(paste("Error in treecompareR::parse_item():",
         "`item` must be a single string"))
  }

  if(length(entities) == 0){
    #if entities is empty, return an empty data.frame.
    return(data.frame())
  }else{

    if(item %in% c("identifier",
                   "smiles",
                   "inchikey",
                   "classification_version",
                   "query_url",
                   "query_status")){
      return(entities[[item]])
    }

      entities <- as.list(entities)
      for(ii in c("identifier",
                    "smiles",
                    "inchikey",
                    "classification_version")){
        if(length(entities[[ii]])==0){
          entities[[ii]] <- NA_character_
        }
      }

      base_dat <- data.frame(entities[c("identifier",
                                  "smiles",
                                  "inchikey",
                                  "classification_version")])

      #each item may be one of several formats.

      # For structure queries:

      #it may be a list with one element for each identifier.
      #if so, each list element may be a data.frame, or a vector.
      #the data.frames may be empty, or not.

      #lists of data frames are in: intermediate_nodes, alternative_parents
      #for intermediate nodes, if some identifiers had them and other did not,
      #the data.frames for the identifiers with no intermediate nodes will have 0 columns and 0 rows,
      #and the data.frames for the identifiers with intermediate nodes will have 4 columns and 1 row.

      #lists of vectors are in: substituents, ancestors,
      #predicted_chebi_terms, predicted_lipidmaps_terms (if present)

      #sometimes the list elements are themselves empty lists, if no entity had a value.
      #This is the case for external_descriptors and predicted_lipidmaps_terms, for everything I've tried so far.
      #It is also the case for intermediate nodes, if no entity had intermediate nodes.

      #it may be a single data.frame with one row for each identifier.
      #this is the case for kingdom, superclass, class, and subclass, and direct_parent.
      #if one identifier has a value for these and another does not,
      #then there will still be one row for each identifier, but it'll be filled with NA if necessary.

      #for kingdom, superclass, class, and subclass, if no identifiers had a value,
      #it will be a vector of NAs with one element for each identifier.

      #it may be a vector with one element for each identifier. In this case,
      #any identifier without a value should have NA inserted.

      #For single-inchikey queries:

      #kingdom, superclass, class, subclass, and direct_parent will be lists,
      #unless that inchikey does not have a classification at that level,
      #in which case the corresponding item will be NULL.

      #intermediate_nodes will be a data.frame, unless that inchikey does not
      #have intermediate nodes, in which case it wil be an empty list.

      #alternative_parents is a data.frame with a varying number of rows.

      #molecular_framework and description are 1-element vectors
      #substituents is a vector of varying length
      #ancestors is a vector of varying length
      #predicted_chebi_terms is a vector of varying length
      this_item <- entities[[item]]

      if(query_type %in% "structure"){
      if(is.list(this_item) &
         !is.data.frame(this_item)){
        #if a list and not a single data.frame:
        #convert each list element to a data.frame
        this_item <- lapply(this_item,
                            as.data.frame)
        #
        #name the list elements after the identifiers
        names(this_item) <- entities$identifier
        #rowbind with identifiers as a new variable
        this_out <- dplyr::bind_rows(this_item,
                                     .id = "identifier")

        #if any list-columns (this happens for external_descriptors),
        #then unnest.
        list_cols <- names(this_out)[sapply(this_out, is.list)]
        if(length(list_cols)>0){
          this_out <- tidyr::unnest(this_out,
                                    cols = tidyr::all_of(list_cols))
        }
        #handle renaming if necessary
        this_out <- df_check(this_out = this_out,
                             item = item,
                             entities_identifier = entities$identifier)
      } else if(is.data.frame(this_item)){
        if(nrow(this_item)>0){
          #this_item should have as many rows as identifiers.
          #add the identifiers as a column.
          this_out <- cbind(data.frame(identifier = entities$identifier,
                                       this_item))
        }else{
          this_out <- data.frame()
        }

        this_out <- df_check(this_out)

      }else{
        #if it's not a list and not a data.frame
        #convert it to a data.frame
        #and add identifier as a new column, if possible
        if(length(this_item)>0){
        this_out <- cbind(data.frame(identifier = entities$identifier,
                                      as.data.frame(this_item)))
        }else{
          this_out <- data.frame()
        }

        #if any list-columns (this happens for external_descriptors),
        #then unnest.
        list_cols <- names(this_out)[sapply(this_out, is.list)]
        if(length(list_cols)>0){
          this_out <- tidyr::unnest(this_out,
                                    cols = tidyr::all_of(list_cols))
        }

        this_out <- df_check(this_out = this_out,
                             item = item,
                             entities_identifier = entities$identifier)
      }
   }else if(query_type %in% "inchikey"){
     #it turns out the steps are the same no matter what the type of the item.

     this_out <- as.data.frame(this_item)

     #if any list-columns (this happens for external_descriptors),
     #then unnest.
     if(length(this_out) > 0){
     list_cols <- names(this_out)[sapply(this_out, is.list)]
     if(length(list_cols)>0){
       this_out <- tidyr::unnest(this_out,
                                 cols = tidyr::all_of(list_cols))
     }
     }

     if(nrow(this_out)>0){
       #this_item should have as many rows as identifiers.
       #add the identifiers as a column.
       this_out <- cbind(data.frame(identifier = entities$identifier,
                                    this_out))
     }else{
       this_out <- data.frame()
     }
     this_out <- df_check(this_out = this_out,
                          item = item,
                          entities_identifier = entities$identifier)
   }else{
     stop("query_type should be either 'structure' or 'inchikey'")
   }

      this_out <- dplyr::left_join(base_dat,
                                   this_out,
                                   by = "identifier")

      return(this_out)
   } #end if/else length(entities)==0
}

#' @title Parsed data frame check
#'
#' @description Helper function used by [parse_item()]
#'
#' @details This is a helper function used by [parse_item()] and generally should not be called directly by the user.
#'
#' @param this_out A `data.frame` consisting of parsed output in one element of the ClassyFire JSON `entities` element
#' @param item Character: the name of the element of the ClassyFire JSON `entities` element
#' @param entities_identifier Character vector: all identifiers in `entities$identifier`
#' @return A `data.frame` either the same as the original `this_out`, or renamed, or with NAs filled in if `this_out` was empty
#' @author Caroline Ring
#' @examples
#' # There is no example.
#'
df_check <- function(this_out, item, entities_identifier){
  non_id_names <- setdiff(names(this_out),
                          "identifier")
  #if the result has "identifier" and columns "name", "description",
  #"chemont_id", and "url", then keep as-is.

  #if the result has "identifier" and one other column, then name that other column after item,
  #except if it should have "name", "description", "chemont_id", and "url"
  if(length(non_id_names) == 1){
    if(item %in% c("kingdom",
                   "superclass",
                   "class",
                   "subclass",
                   "intermediate_nodes",
                   "direct_parent",
                   "alternative_parents")){
      this_out <- data.frame(identifier = entities_identifier,
                             name = NA_character_,
                             description = NA_character_,
                             chemont_id = NA_character_,
                             url = NA_character_)
    }else if(item %in% "external_descriptors"){
      this_out <- data.frame(identifier = entities_identifier,
                             source = NA_character_,
                             source_id = NA_character_,
                             annotations = NA_character_)
    }else{
      this_out <- setNames(this_out, c("identifier", item))
    }

  }
  #if the result has only the "identifier" column and zero rows,
  #then it means all list elements were empty.
  if(length(non_id_names) == 0 &
     nrow(this_out) == 0){
    #if this is for kingdom, superclass, class, subclass,
    #intermediate_nodes, direct_parent, or alternative_parents,
    #then fill in with "name", "description", "chemont_id", and "url" set to NA for all identifiers.
    if(item %in% c("kingdom",
                   "superclass",
                   "class",
                   "subclass",
                   "intermediate_nodes",
                   "direct_parent",
                   "alternative_parents")){
      this_out <- data.frame(identifier = entities_identifier,
                             name = NA_character_,
                             description = NA_character_,
                             chemont_id = NA_character_,
                             url = NA_character_)
    }else if(item %in% "external_descriptors"){
      this_out <- data.frame(identifier = entities_identifier,
                             source = NA_character_,
                             source_id = NA_character_,
                             annotations = NA_character_)
    }else{
      #Otherwise, fill in with a column named after item, set to NA for all identifiers.
      this_out <- data.frame(identifier = entities_identifier,
                             NA_character_)
      this_out <- setNames(this_out,
                           c("identifier", item))
    }
  } #end if(length(non_id_names) == 0 & nrow(this_out) == 0)

  return(this_out)
}

#' Is InChIKey
#'
#' Check if a string is a validly-formatted InChIKey
#'
#' @details The criteria are the following, as described in Heller et al. (2015)
#'   (see References).
#'
#' - must have 27 characters
#' - the first 14 characters must be uppercase letters
#' - the 15th character must be a hyphen
#' - the 16th-23rd characters must be uppercase letters
#' - the 24th character must be either "S" or "N"
#' - the 25th character must be "A"
#' - the 26th character must be a hyphen
#' - the 27th character must be an uppercase letter
#'
#'
#' @param x Character: A vector of one or more strings to check for InChIKey
#'   format.
#' @return A logical vector the same length as `x`, containing `TRUE` if valid
#'   InChIKey format, and `FALSE` otherwise.
#' @examples
#' is_inchikey("PLDWAJLZAAHOGG-UHFFFAOYSA-N") #returns TRUE
#' is_inchikey("PLDWAJLZAAHOGG-UHFFFAOYSB-N") #returns FALSE
#' is_inchikey("BAD-INCHIKEY-X") #returns FALSE
#' is_inchikey(5.234) #returns FALSE
#' is_inchikey(NA_character_) #returns FALSE
#' is_inchikey(NULL) #returns FALSE
#'
#' @export
#' @author Caroline Ring
#' @references Heller, S.R., McNaught, A., Pletnev, I. et al. InChI, the IUPAC
#'   International Chemical Identifier. J Cheminform 7, 23 (2015).
#'   https://doi.org/10.1186/s13321-015-0068-4

is_inchikey <- function(x){
  if(is.character(x)){
    return(
      grepl(x=x,
            pattern = "^[[:upper:]]{14}\\-[[:upper:]]{8}[SN]A\\-[[:upper:]]$")
    )
  }else{
    return(rep(FALSE, length(x)))
  }
}

