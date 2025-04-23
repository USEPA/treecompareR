
#'@title InChIKey classification
#'
#'@description This function uses the ClassyFire API to classify chemicals from
#'  an input data.table using the InChIKey chemical identifier.
#'
#'@details This function queries ClassyFire's lookup table of pre-classified
#'  InChiKeys to get classifications, using the ClassyFire API.
#'
#'@param inchikeys Character: A vector of InCHiKeys to be classified
#'@param tax_level_labels By default, the list of taxonomy levels for
#'  ClassyFire: \code{kingdom, superclass, class, subclass, level5, ...
#'  level11}.
#'@param wait_sec Numeric: A parameter controlling how many seconds between qqueries sent
#'  to the ClassyFire API server. Default 5, to respect the limit of 12 requests
#'  per second. Setting it any lower than 5 will result in failed queries.
#'@return A `data.frame with the following variables`: \itemize{
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
#'  }
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
#'  Note also that the returned data.frame includes only unique, valid
#'  InCHiKeys. Any duplicates, blanks, NAs, or anything that is not a valid
#'  InCHiKey is not queried, and is not
#'  included in the output.
#'
#'  If you have a source `data.frame` with duplicate, missing, or invalid
#'  InCHiKeys, you can merge the returned `data.frame` with it (or pivot wider,
#'  then merge). For example, if your source `data.frame` is called
#'  \code{source_dat} with variable \code{"INCHIKEY"} containing the InCHiKeys,
#'  and the returned `data.frame` is in variable \code{dat}, the following code
#'  will do the merge: \code{dplyr::left_join(source_dat, dat, by = "INCHIKEY" =
#'  "identifier")}. In that case, classification columns will be filled with NA
#'  for any missing or invalid InCHiKeys.
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
      this_output <- query_classyfire_inchikey(inchikey = this_inchikey,
                                               wait_sec = wait_sec)
      this_output$identifier <- this_inchikey
      return(this_output)
    } #end function to apply to each INCHIKEY
  )

  class_list <- lapply(entities_list,
                       parse_classified_entities,
                       tax_level_labels = tax_level_labels)

  #now rowbind the list of data.frames to get one big data.frame
  inchi_class <- dplyr::bind_rows(class_list)

  #get any input inchikeys that were excluded
  #or are otherwise missing from the classified output
  #Get any input structures not handled by ClassyFire
  #e.g. any duplicates, or batches that failed, etc.
  #keep their placeholder values.
  if("identifier" %in% names(this_item_df)){
    handled_id <- inchi_class$identifier
  }else{
    handled_id <- character(0)
  }
  missing_entities <- output |>
    dplyr::filter(!(identifier %in% handled_id)) |>
    dplyr::mutate(entity_type = "not handled by ClassyFire")

  #Output: bind things handled by ClassyFire and things not handled
  inchi_class <- dplyr::bind_rows(missing_entities,
                                 inchi_class) |>
    as.data.frame()  #convert from tibble to data.frame

  inchi_class <- dplyr::bind_rows(inchi_class,
                                  missing)

  return(inchi_out)
}

#' Classify data.table
#'
#' @param data A data.table with rows corresponding to chemicals. A column named
#' `INCHIKEY` containing the InChIKey for each chemical is required.
#'
#' @return A data.table containing the original data as well as the ClassyFire
#' classification data for each chemical, when available.
#' @export
#'
#' @examplesIf FALSE
#' # Chemical information for Bisphenol A
#' bpa <- data.table(PREFERRED_NAME = 'Bisphenol A',
#'                   CASRN = '80-05-7',
#'                   INCHIKEY = 'IISBACLAFKSPIT-UHFFFAOYSA-N')
#' bpa_classified <- classify_datatable(bpa)
#' bpa_classified
classify_datatable <- function(data) {
  INCHIKEY <- NULL
  if (!is.data.frame(data)){
    stop("Input data must be a data.table!")
  }

  if (!data.table::is.data.table(data)){
    warning('Casting input data.frame to data.table')
    datatable <- data.table::data.table(data)
  } else {
    datatable <- data.table::copy(data)
  }

  if (!('INCHIKEY' %in% names(datatable))){
    stop('Input data must include a column of InChIKeys name `INCHIKEY`!')
  }

  # Grab unique inchikeys
  inchikeys <- datatable[, unique(INCHIKEY)]

  # Classify inchikeys
  inchikeys_classified <- classify_inchikeys(inchikeys = inchikeys)

  # Cast as a data.table
  inchikeys_classified <- data.table::data.table(inchikeys_classified)

  # Default format of ClassyFire results
  default_table <- data.table::data.table(identifier = character(),
                                          smiles = character(),
                                          inchikey = character(),
                                          report = character(),
                                          level = character(),
                                          name = character(),
                                          description = character(),
                                          chemont_id = character(),
                                          url = character())
  default_names <- names(default_table)

  if (dim(inchikeys_classified)[[1]] == 0){
    return(data)
  }

  # Convert from long to wide format
  wide_classifications <- dcast(
    inchikeys_classified,
    formula = identifier + smiles + inchikey + report ~ level,
    value.var = 'name'
  )

  # Bind with default results
  temp_table <- rbindlist(list(default_table, wide_classifications),
                          fill = TRUE)
  setcolorder(temp_table, neworder = default_names)

  final_table <- data.table::merge.data.table(x = datatable,
                                              y = temp_table,
                                              by.x = 'INCHIKEY',
                                              by.y = 'identifier',
                                              all = TRUE)
  return(final_table)
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
#' @param ... Other arguments as for \code{\link{query_classyfire}}.
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
      json_parse <- do.call(query_classyfire,
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

      #Parse classified entities. If there were none, parse_classified_entities()
      #will just return a data.frame with zero rows.
      classified <- lapply(entities_list,
                           parse_classified_entities) |>
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
                            this_parse <- parse_list_item(
                              entities = this_page,
                              item = this_item,
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

#' Send ClassyFire query
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
#'   command or query. Default 5 seconds (because ClassyFire requests that POST
#'   requests be limited to 12 per second). If less than 5 seconds, will be
#'   reset to 5 seconds, with a warning. Wait times for multiple retries
#'   use exponential backoff with full jitter (wait time is \code{runif(1,
#'   wait_sec, wait_sec * 2^(attempt_number)}).
#' @param retry_query_times Max number of attempts to retry an "In Queue" or
#'   "Processing" query (with wait time in between tries). Default 10.
#' @param processing_wait_per_input Minimum number of seconds to wait *per input
#'   string* before retrying to retrieve results for a query whose status is
#'   "Processing." Default NULL results in \code{wait_sec/length(input)}. Uses
#'   exponential backoff: effective wait time per input =
#'   \code{processing_wait_per_input * 2^(attempt number)}. Total wait time is
#'   either the effective wait time per input times the number of input
#'   structures, or \code{wait_sec} seconds, whichever is greater.
#' #'@param terminate_on Integer vector: List of HTTP status codes which will
#'  immediately terminate with no more retries. Default: ` c(400:407, 409:418,
#'  421:428, 431, 451, 500:511)`
#'@param base_url The base URL for ClassyFire API queries. Default
#'  `"http://classyfire.wishartlab.com/entities"`. If you change this, the query
#'  is highly unlikely to work!
#' @return A list of lists. The outer list has one element for each page of the
#'   JSON output from the ClassyFire query (there is one page for every ten
#'   entities). If the ClassyFire query failed or timed out, the list will have
#'   one element. Each list element is itself a list, consisting of the parsed
#'   JSON result for each page. A completed query will have named elements `id`,
#'   `query_url`, `query_status`, `label`, `classification_status`,
#'   `number_of_elements`, `number_of_pages`, `invalid_entities`, and
#'   `entities`. A failed query will have named elements `id`, `query_url`,
#'   `query_status`, `label`, `classification_status`, and `number_of_pages`.
#'   `id` gives the numerical query ID (assigned by ClassyFire). `query_url` gives the
#'   queried URL. `label` gives the user-supplied query label. `query_status`
#'   gives the HTTP status of the request (200 means successful).
#'   `classification_status` reports the ClassyFire status: "Done" means
#'   classification was successfully completed; "In Queue" means the query timed
#'   out while it was still queued; "In Progress" means the query timed out
#'   while it was still processing; "Failed" means the query failed
#'   before ClassyFire could report a status (e.g. HTTP status code other than
#'   200, or ClassyFire returned some content that did not include a
#'   classification status). `number_of_elements` gives the number of classified
#'   entities. `number_of_pages` gives the total number of pages for this query.
#'   `invalid_entities`, if present, is a `data.frame` listing queried entity
#'   identifiers that ClassyFire found invalid. `entities`, if present, is a
#'   nested `data.frame` giving classifications.
#'
#'   @examples
#'      query_classyfire(input = c("COC1=CC(Br)=CC=C1", "SCCSCCS"))
#'
#'
query_classyfire <- function(input = NULL,
                             url = NULL,
                             label = "query",
                             type = "STRUCTURE",
                             retry_get_times = 10,
                             wait_sec = 5,
                             retry_query_times = 10,
                             processing_wait_per_input = NULL,
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
    stop(paste0("Error in treecompareR::query_classyfire():",
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
#'   API query (e.g, from the results of [query_classyfire()])
#' @param tax_level_labels ChemOnt taxonomy level labels (default
#'   \code{\link{chemont_tax_levels}})
#' @return A data.frame consisting of rows corresponding to each classified
#'   entry from the `entities` input. A data.frame with one row for each entity
#'   identifier, and variables `identifier`, `smiles`, `inchikey`, `kingdom`,
#'   `superclass`, `class`, `subclass`, `level5`, `level6`, ... `level11`,
#'   `classification_version`, and `report`.
#' @examples
#' #get a set of results to parse
#' my_results <- query_classyfire(input = c("COC1=CC(Br)=CC=C1", "SCCSCCS"))
#' #parse classifications
#' parse_classified_entities(entities = my_results[[1]]$entities)
#'
parse_classified_entities <- function(entities,
                                      tax_level_labels = chemont_tax_levels){

  #check to see whether classifications actually exist for these
  if(length(entities)==0){
    #if no classifications, return empty data.frame,
    #but with the expected variable names
    classified_entities <- data.frame(identifier = character(0),
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
    #expected named items in entities and their expected classes:
    # identifier                    smiles                  inchikey
    # "character"               "character"               "character"
    # kingdom                superclass                     class
    # "data.frame"              "data.frame"              "data.frame"
    # subclass        intermediate_nodes             direct_parent
    # "data.frame"                    "list"              "data.frame"
    # alternative_parents       molecular_framework              substituents
    # "list"               "character"                    "list"
    # description      external_descriptors                 ancestors
    # "character"                    "list"                    "list"
    # predicted_chebi_terms predicted_lipidmaps_terms    classification_version
    # "list"                    "list"               "character"


    #get identifiers, smiles, inchikey
    #these will be character vectors identifying entities

    #convert from data.frame to list, if not already
    entities <- as.list(entities)

    #note that some of these may be character(0) -- handle
    #by filling with NAs
    for(item in c("identifier",
                  "smiles",
                  "inchikey",
                  "classification_version")){
      if(length(entities[[item]])==0){
        entities[[item]] <- NA_character_
      }
    }

    cf <- data.frame(entities[c("identifier",
                                "smiles",
                                "inchikey",
                                "classification_version")])

    #kingdom, superclass, class, subclass are all data.frames
    #with one row for each identifier,
    #and variables name, description, chemont_id, url
    #giving the relevant ChemOnt taxonomy label for each identifier
    #(i.e., "name" in the "kingdom" element gives the kingdom for each identifier,
    #"name" in the "superclass" element gives the superclass for each identifier,
    #etc.)
    #rowbind these

    #note: if no classification for one of these four levels, the item may be NA or NULL.
    #handle accordingly.
    for(item in c("kingdom",
                  "superclass",
                  "class",
                  "subclass")){
      #check if NA, NULL, or length-0 -- indicates no classification
      if(all(is.na(entities[[item]])) |
         is.null(entities[[item]]) |
         length(entities[[item]]) %in% 0){
        #replace with empty data.frame with the required names
        entities[[item]] <- data.frame(name = character(0),
                                       description = character(0),
                                       chemont_id = character(0),
                                       url = character(0))
      }else if(!is.data.frame(entities[[item]])){
        #if not a data.frame --
        #e.g. if the query was for an inchikey rather than a structure,
        #kingdom/superclass/class/subclass are provided as lists --
        #then coerce it to a data.frame.
        entities[[item]] <- as.data.frame(entities[[item]])
      }

      #add variable with identifiers
      if(nrow(entities[[item]])>0){
        entities[[item]][["identifier"]] <- entities$identifier
      }else{
        entities[[item]]$identifier <- character(0)
      }
    }

    cf_class1 <- dplyr::bind_rows(entities[c("kingdom",
                                             "superclass",
                                             "class",
                                             "subclass")])

    #"intermediate_nodes"
    #is a list with one element for each item in "input",
    #a 1-row data.frame with variables name, description, chemont_id, url
    #giving more-specific levels of classification, if any.
    #if no intermediate nodes for an entity,
    #its corresponding element in "intermediate_nodes" is an empty data.frame.
    #if it is an empty list (as happens for single-InChiKey queries without intermediate nodes)
    #then substitute an empty data-frame
    cf_class2 <- entities$intermediate_nodes
    if(is.data.frame(cf_class2)){
      #for an inchikey query, it will be a single data.frame, not a list of them.
      #just add the identifier column.
      cf_class2_df <- cf_class2
      cf_class2_df$identifier <- entities$identifier
    }else{
      #for a structure query, it will be a list of data.frames
      #(even if there was only one entity -- it'll be a one-element list)
      if(length(cf_class2)>0){
        names(cf_class2) <- entities$identifier
        cf_class2_df <- dplyr::bind_rows(cf_class2, .id = "identifier")
      }else{
        cf_class2_df <- data.frame(name = character(0),
                                   description = character(0),
                                   chemont_id = character(0),
                                   url = character(0),
                                   identifier = character(0))
      }
    }
    rm(cf_class2)

    #rowbind the kingdom-thru-superclass labels,
    #and the "intermediate nodes" labels.
    cf_class <- dplyr::bind_rows(cf_class1,
                                 cf_class2_df)

    # element 'direct_parent' is a data.frame
    #with one row for each identifier,
    #and variables name, description, chemont_id, url
    direct <- entities$direct_parent
    if(!is.data.frame(direct)){
      #for an inchikey query it'll be a list rather than a data.frame,
      #so coerce it
      direct <- as.data.frame(direct)
    }
    if(nrow(direct) > 0){
      direct$identifier <- entities$identifier
    }else{
      direct$identifier <- character(0)
    }
    cf_class <- dplyr::bind_rows(cf_class, direct)

    #"direct parent" is the terminal label,
    #and may be a duplicate of kingdom, superclass, class, subclass,
    #or something in "intermediate nodes".

    #To remove duplicates,
    #get level numbers for all labels in the ChemOnt tree,
    #merge them into cf_class,
    #and take unique rows only.
    chemont_tree_df <- get_tree_df(chemont_tree)

    cf_class <- dplyr::left_join(cf_class,
                                 chemont_tree_df[,
                                                 c("Name",
                                                   "level")],
                                 by = c("name" = "Name")) |>
      dplyr::select(identifier, level, name) |>
      dplyr::distinct() #keep only unique rows

    #Change the levels from numbers to taxonomy level names
    #e.g. 1, 2, 3, 4 -> kingdom, superclass, class, subclass
    cf_class <- cf_class |>
      dplyr::mutate(level = tax_level_labels[level])

    #We have found at least one instance where the ChemOnt tree listed something at the wrong level
    #(see data-raw/make_chemont_taxonomy.R).
    #Though we manually fixed that instance, there may be others.
    #If so, cf_class will have more than one "name" at a given "level."

    # #find cases of multiple name at given level, if any
    #
    # dup_rows <- cf_class |>
    #   dplyr::summarise(n = dplyr::n(), .by = c(identifier, level)) |>
    #   dplyr::filter(n > 1L) |>
    #   as.data.frame()
    #
    # if(nrow(dup_rows)>0){
    #   warning(paste("There are multiple labels listed at the same level",
    #                 "for one or more identifiers. Duplicates:\n",
    #                 paste(capture.output(print(dup_rows)), collapse = "\n")
    #   ))
    #
    #   #reshape wider
    #   suppressWarnings(cf_class_wide <- cf_class |>
    #                      tidyr::pivot_wider(id_cols = "identifier",
    #                                         names_from = "level",
    #                                         values_from = "name") |>
    #                      tidyr::unnest(cols = dplyr::pick(
    #                        tidyr::all_of(tax_level_labels)
    #                      )
    #                      )
    #   )
    # }else{
    #   cf_class_wide <- cf_class |>
    #     tidyr::pivot_wider(id_cols = "identifier",
    #                        names_from = "level",
    #                        values_from = "name")
    # }
    #
    # classified_entities <- cf_class_wide

    classified_entities <- cf_class

    #merge in smiles and inchikey
    classified_entities <- dplyr::left_join(cf,
                                            classified_entities,
                                            by = "identifier") |>
      as.data.frame()
  } #end if length(entities)>0

  return(classified_entities)
}

#' @title Parse list item
#'
#' @description Parse an item from a list of classified entities.
#'
#' @details This is a helper function used by [classify_structures()] and
#'   [classify_inchikeys()]. It generally should not be called directly by the
#'   user.
#'
#'   The main classification is parsed by [parse_classified_entities()], but
#'   ClassyFire returns several additional pieces of classification output.
#'   These include:
#'
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
#' - `alternative_parents`: `name`, `description`, `chemont_id`, `url`
#' - everything else: a single additional variable named after the requested item (i.e., `molecular_framework`, `substituents`, etc.)
#'
#'
#' @param entities A `data.frame` of classified entities from one page of
#'   ClassyFire results, e.g, the `entities` element of one element of the list
#'   output by [query_classyfire()].
#' @param item Character: The name of an item in `entities` that is a list the
#'   same length as `entities$identifier`. One of `'alternative_parents'`,
#'   `'molecular_framework'`, `'substituents'`, `'description'`,
#'   `'external_descriptors'`, `'ancestors'`, `'predicted_chebi_terms'`,
#'   `'predicted_lipidmaps_terms'`. See Details. If you provide
#'   `'classification'`, this function will call [parse_classified_entities()].
#'   These are all the named items currently in `entities` as returned by
#'   ClassyFire, so providing any other string here will result in an empty
#'   `data.frame` being returned.
#' @return A `data.frame` containing the information from the specified item.
#'   Contains variables `identifier`, `smiles`, `inchikey`,
#'   `classification_version`, and other variables depending on what item was
#'   specified. If `item = 'alternative_parents'`, additional variables are
#'   `name`, `description`, `chemont_id`, `url`. For any of the other possible
#'   values of `item` listed above (except `'classification'`), there will be
#'   one additional variable named for `item` (e.g., if `item =
#'   'molecular_framework'`, the additional variable will be named
#'   `molecular_framework`). If `item = 'classification'`, see
#'   [parse_classified_entities()] for details on the returned `data.frame`. If
#'   `item` is some string other than those listed above, i.e., not matching any
#'   of the named items in `entities`, then an empty `data.frame` will be
#'   returned (zero rows and no variables).
#' @examples
#' #get a set of results to parse
#' my_results <- query_classyfire(input = c("COC1=CC(Br)=CC=C1", "SCCSCCS"))
#' #get alternative parents
#' parse_list_item(entities = my_results[[1]]$entities, item = "alternative_parents")
#' #get molecular framework
#' parse_list_item(entities = my_results[[1]]$entities, item = "molecular_framework")
#' @author Caroline Ring
parse_list_item <- function(entities,
                            item,
                            tax_level_labels = chemont_tax_levels){

  if(length(item) > 1){
    stop(paste("Error in treecompareR::parse_list_item():",
         "`item` must be a single string"))
  }

  if(item %in% "classification"){
    return(
      parse_classified_entities(
        entities = entities,
        tax_level_labels = tax_level_labels
      )
    )
  }

  if(length(entities) == 0){
    #if entities is empty, return an empty data.frame.
    return(data.frame())
  }else{

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

    if(!(item %in% names(entities))){
      item_out <- NULL
    }else{
      #check to see if the specified item is a list (not a single data.frame)
      if(is.list(entities[[item]]) &
         !is.data.frame((entities[[item]]))){
        #if a list, is it the same length as entities$identifier?
        if(length(entities[[item]] %in% length(entities$identifier))){
          entities[[item]] <- lapply(
            entities[[item]],
            function(this_element){
              if(!is.data.frame(this_element)){
                #if the elements of entities[[item]] are not themselves data frames,
                #but are vectors, matrixes, or arrays,
                #coerce them to data.frames.
                tmp <- as.data.frame(this_element)
                if(length(tmp)==1){
                  #i.e. if data.frame with one variable
                  tmp <- setNames(tmp,
                           item)
                }else if(length(tmp)==0){
                  #if data.frame with no variables,
                  #i.e. if elements were empty lists,
                  #substitute placeholder NA
                  tmp <- data.frame(NA_character_)
                  tmp <- setNames(tmp,
                           item)
                }else{
                  #if not a data.frame,
                  #but neither a vector nor empty,
                  #e.g. if element is a list, matrix, or array,
                  #then coerce to a data.frame.
                  as.data.frame(tmp)
                }
              }else{
                #if the elements of entities[[item]] are data frames,
                #just keep them
                this_element
              }
            }
          )
           #entities[[item]] should now be a list of data.frames
          names(entities[[item]]) <- entities$identifier

          item_out <- dplyr::bind_rows(entities[[item]],
                                       .id = "identifier") |>
            dplyr::left_join(base_dat,
                             by = "identifier") |>
            dplyr::relocate(identifier,
                            smiles,
                            inchikey,
                            classification_version) |>
            as.data.frame()
        }else if(length(entities[[item]]) == 0){
          #if the item is an empty list,
          #return an empty data.frame
          item_out <- data.frame()
          }else{
          #if the item is a non-empty list a different length from identifiers,
          #then something has probably gone wrong
          #just return it as-is
          warning(paste("Item", paste0("'", item, "'"), "in `entities`",
                        "is a list, but doesn't seem to be the same length as 'identifier'.",
                        "This probably means something has gone wrong.",
                        "Row-binding it into a single data.frame."))
          item_out <- dplyr::bind_rows(entities[[item]])
        }
      }else if(is.data.frame(entities[[item]])){
        #if it is a single data.frame, rather than a list of data.frames,
        #then if it has the same number of rows as identifiers,
        #or if there is only one identifier,
        #just add the identifiers.
        if(nrow(entities[[item]]) %in% length(entities$identifier) |
           length(entities$identifier) == 1){
          item_out <- entities[[item]] |>
            dplyr::mutate(identifier = entities$identifier) |>
            dplyr::left_join(base_dat,
                             by = "identifier") |>
            dplyr::relocate(identifier,
                            smiles,
                            inchikey,
                            classification_version) |>
            as.data.frame()
        }else if(nrow(entities[[item]]) == 0){
          item_out <- entities[[item]]
          }else{ #if it's non-empty but has a different number of rows from identifiers
          #just return it as-is
          message(paste("Item", paste0("'", item, "'"), "in `entities`",
                        "is a data.frame but doesn't seem to have",
                        "the same number of rows as 'identifiers'.",
                        "This probably means something has gone wrong."))
          item_out <- entities[[item]]
        }
      }else{ #if neither a list nor a data.frame
        #try to coerce it to a data.frame
        item_out <- as.data.frame(entities[[item]])
        if(length(item_out)==1){
          #i.e. if data.frame with one variable
          item_out <- setNames(item_out,
                          item)
        }
        if(nrow(item_out) > 0){
        if(nrow(item_out) %in% length(entities$identifier) |
           length(entities$identifier) == 1){
          item_out$identifier <- entities$identifier
          #move identifier column to the front
          item_out <- item_out |>
            dplyr::left_join(base_dat,
                             by = "identifier") |>
            dplyr::relocate(identifier,
                            smiles,
                            inchikey,
                            classification_version) |>
            as.data.frame()
        }else{
          message(paste("Item", paste0("'", item, "'"), "in `entities`",
                        "could be coerced to a data.frame, but it doesn't seem to have",
                        "the same number of rows as 'identifiers'.",
                        "This probably means something has gone wrong."))
        }
        }
      } #end if/else for type of item: list, data.frame, or other
    } #end if/else item in names(entities)
    return(item_out)
  } #end if/else length(entities)==0
}

#'@title Query ClassyFire InChIKey
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
#' query_classyfire_inchikey(inchikey = "PLDWAJLZAAHOGG-UHFFFAOYSA-N")
#'
query_classyfire_inchikey <- function(inchikey,
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
