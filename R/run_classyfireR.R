
#'InChIKey classification
#'
#'This function uses the ClassyFire API to classify chemicals from an input
#'data.table using the InChIKey chemical identifier.
#'
#'This function queries ClassyFire's lookup table of pre-classified InChiKeys
#'
#'@param inchikeys A vector of InCHiKeys to be classified
#'@param tax_level_labels By default, the list of taxonomy levels for
#'  ClassyFire: \code{kingdom, superclass, class, subclass, level5, ...
#'  level11}.
#'@param wait_min A parameter controlling how many seconds between qqueries sent
#'  to the ClassyFire API server.
#'@return A data.frame with the following variables: \itemize{
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
#'  InCHiKey according to `webchem::is.inchikey()`is not queried, and is not
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
                               wait_min = 0.5){


  INCHIKEYS <- unique(inchikeys) #save time by removing duplicates

  #get classifications as a list of data.tables, one for each INCHIKEY
  entities_list <- lapply(
    INCHIKEYS,
    function(this_inchikey){
      #if inchikey is NA, blank, or not in valid InChiKey format,
      #do not query API; return NA classification instead
      if(is.na(this_inchikey) || #if InChiKey is NA
         !nzchar(trimws(this_inchikey)) || #if InChiKey is blank
         !webchem::is.inchikey(this_inchikey, #if not valid InChiKey format
                               type = "format")
         ){
        #placeholder blank output list with all expected elements
        output <- list()
      }else{ #if inchikey is valid, query ClassyFire API
        Sys.sleep(wait_min) #pause minimum number of seconds before querying again
        output <- query_classyfire_inchikey(inchikey = this_inchikey,
                                            wait_min = wait_min)
        output$identifier <- this_inchikey
      }
      return(output)
    } #end function to apply to each INCHIKEY
  )

  class_list <- lapply(entities_list,
                       parse_classified_entities,
                       tax_level_labels = tax_level_labels)

  #now rowbind the list of data.frames to get one big data.frame
  inchi_class <- dplyr::bind_rows(class_list)

  #now order as for input -- re-inserting any duplicates
  # new_class <- inchi_class[match(inchikeys, inchi_class$identifier),
  #                          ]

  #result will be long-format

  return(inchi_class)
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
                                          kingdom = character(),
                                          superclass = character(),
                                          class = character(),
                                          subclass = character(),
                                          level5 = character(),
                                          level6 = character(),
                                          level7 = character(),
                                          level8 = character(),
                                          level9 = character(),
                                          level10 = character(),
                                          level11 = character())
  default_names <- names(default_table)

  if (dim(inchikeys_classified)[[1]] == 0){
    return(data)
  }

  # Convert from long to wide format
  wide_classifications <- dcast(inchikeys_classified, formula = identifier + smiles + inchikey + report ~ level, value.var = 'name')

  # Bind with default results
  temp_table <- rbindlist(list(default_table, wide_classifications), fill = TRUE)
  setcolorder(temp_table, neworder = default_names)

  final_table <- data.table::merge.data.table(x = datatable,
                                              y = temp_table,
                                              by.x = 'INCHIKEY',
                                              by.y = 'identifier',
                                              all = TRUE)
  return(final_table)
}


#' Query ClassyFire by structure
#'
#' This function takes a vector of structural identifiers (SMILES strings or
#' InChi strings) and queries the ClassyFire API to get classifications for each
#' one.
#'
#' @param input A character vector of structural identifiers: SMILES strings or
#'   InChi strings. May optionally be named. If so, the names will be returned as a
#'   column named  \code{identifier} in the output data.frame. If not named, the
#'   vector itself will be returned as a column named  \code{identifier} in the
#'   output data.frame
#' @param tax_level_labels By default, the list of taxonomy levels for
#'   ClassyFire: \code{kingdom, superclass, class, subclass, level5, ...
#'   level11}.
#' @param ... Other arguments as for \code{\link{query_classyfire}}.
#' @return A data frame with ClassyFire classifications for each input
#'   structural identifier. Will contain columns \code{identifier},
#'   \code{smiles}, \code{inchikey}, one column for each taxonomy level (defined
#'   in argument \code{tax_level_labels}), and \code{report}. The final
#'   \code{report} column explains what happened if a classification could not
#'   be obtained.
#'
#' @export
#' @importFrom magrittr %>%
#' @seealso \code{\link{classify_inchikeys}}
#'
  classify_structures <- function (input = NULL,
                                   tax_level_labels = chemont_tax_levels,
                                   ...){
    identifier <- NULL

    #Names of input structures are their identifiers.
    #If no names, use the structures themselves as identifiers.
    if(is.null(names(input))){
      names(input) <- input
    }


    #Now, query ClassyFire with these structures.
    json_parse <- do.call(query_classyfire,
                          args = c(list(input = input,
                                        url = NULL),
                                   ...)
    )

    #json_parse will be a list, one element for each page of the results.
    #or at least one element even if there were no results.

    #default output to return if classification failed:
    #prepare a named list of NA columns, one for each taxonomy level
    cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
    names(cols_list) <- tax_level_labels
    #put these as columns of a data.frame, along with input structures.
    output <- do.call(data.frame,
                      c(list("identifier" = names(input),
                             "structure" = input),
                        cols_list
                        )
                      )
    #add columns for ClassyFire-derived smiles, inchikey, and report.
    output$smiles <- NA_character_
    output$inchikey <- NA_character_
    output$report <- NA_character_

    #Now try to parse the ClassyFire output in json_parse into a data.table.

        #first pull list of invalid entities, if any.
        #this will be the same for all pages, so just pull from the first.
    #This will be a data.frame or list.
        invalid <- json_parse[[1]]$invalid_entities
        if(length(invalid)==0){ #if no invalid entities
          #return empty data.frame but with the expected variables
          invalid <- data.frame(identifier = character(0),
                                structure = character(0),
                                smiles = character(0),
                                inchikey = character(0),
                                kingdom = character(0),
                                superclass = character(0),
                                class = character(0),
                                subclass = character(0),
                                level5 = character(0),
                                level6 = character(0),
                                level7 = character(0),
                                level8 = character(0),
                                level9 = character(0),
                                level11 = character(0),
                                report = character(0))
        }else{
          #Add identifier column (the names of the inputs that were invalid structures)
          invalid$identifier <- names(input)[input %in% invalid$structure]
          #Add NA columns for taxonomy levels, to denote no classification
          invalid[tax_level_labels] <- NA_character_
          invalid$smiles <- NA_character_
          invalid$inchikey <- NA_character_
          invalid$report2 <- sapply(invalid$report, paste, collapse = " | ")
          invalid$report <- NULL
          invalid$report <- invalid$report2
          invalid$report2 <- NULL
          invalid <- as.data.frame(invalid) #in case it was a list
        } # if(length(invalid)>0)

        #pull list of valid entities (i.e., those with classifications)
        entities_list <- lapply(json_parse,
                                function(x) x$entities)

        #parse classified entities JSON
        classified <- lapply(entities_list,
               parse_classified_entities) %>%
        dplyr::bind_rows() %>%
          dplyr::mutate(structure = input[identifier])

        #Start constructing output:
        #entities handled by ClassyFire (either classified or reported invalid)
        cf_entities <- dplyr::bind_rows(classified,
                                   invalid)

        #Get any entities not handled by ClassyFire
        missing_entities <- dplyr::anti_join(output,
                                             cf_entities,
                                             by = c("identifier",
                                                    "structure"))
        # debugging the code
        print("missing entities")
        print(missing_entities)
        print("\n")
        print("classified entities")
        print(cf_entities)
        print("\n")
        # debugging the code

        #Output: bind classified, invalid, and missing entities
        output <- dplyr::bind_rows(missing_entities,
                                   cf_entities) %>%
          as.data.frame() #convert from tibble to data.frame

        # debugging the code
        print('output')
        print(output)
        # debugging the code

        #set order the same as input
        rownames(output) <- NULL
        output <- output[match(names(input), output$identifier), ]

        # debugging the code
        print('output, again')
        print(output)
        # debugging the code

        #add informational columns
        output[c("id",
                 "label",
                 "classification_status",
                 "url",
                 "status_code")] <- json_parse[[1]][c("id",
                          "label",
                          "classification_status",
                          "url",
                          "status_code")]

    return(output)
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
#' @param type String giving ClassyFire query type. Default \code{"STRUCTURE"}.
#' @param retry_get_times Max number of times to retry GET command to ClassyFire
#'   API to get status of query. Default 10.
#' @param wait_min Minimum time to wait before retrying any GET command or
#'   query. Default 5 seconds (because ClassyFire requests that POST requests be
#'   limited to 12 per second). Highly recommend not changing this to be less
#'   than 5 seconds, to be respectful of the ClassyFire server. All wait times
#'   use exponential backoff with full jitter (\code{wait time = runif(1,
#'   wait_min, wait_min * 2^(attempt number)}).
#' @param retry_query_times Max number of attempts to retry an "In Queue" or
#'   "Processing" query (with wait time in between tries). Default 3.
#' @param processing_wait_per_input Minimum number of seconds to wait per input
#'   string before retrying to retrieve results for a query whose status is
#'   "Processing." Default NULL results in \code{wait_min/length(input)}. Uses
#'   exponential backoff: effective wait time per input =
#'   \code{processing_wait_per_input * 2^(attempt number)}. Total wait time is
#'   either the effective wait time per input times the number of input
#'   structures, or \code{wait_min} seconds, whichever is greater.
#' @return A list of lists. The outer list has one element for each page of the
#'   JSON output from the ClassyFire query (there is one page for every ten
#'   entities). If the ClassyFire query failed or timed out, the list will have
#'   one element. Each list element is itself a list, consisting of the parsed
#'   JSON result for each page. A completed query will have named elements `id`,
#'   `url`, `status_code`, `label`, `classification_status`,
#'   `number_of_elements`, `number_of_pages`, `invalid_entities`, and
#'   `entities`. A failed query will have named elements `id`, `url`,
#'   `status_code`, `label`, `classification_status`, and `number_of_pages`.
#'   `id` gives the numerical query ID (assigned by ClassyFire). `url` gives the
#'   queried URL. `label` gives the user-supplied query label. `status_code`
#'   gives the HTTP status of the request (200 means successful).
#'   `classification_status` reports the ClassyFire status: "Done" means
#'   classification was successfully completed; "In Queue" means the query timed
#'   out while it was still queued; "In Progress" means the query timed out
#'   while it was still processing; "No classification" means the query failed
#'   before ClassyFire could report a status (e.g. HTTP status code other than
#'   200, or ClassyFire returned some content that did not include a
#'   classification status). `number_of_elements` gives the number of classified
#'   entities. `number_of_pages` gives the total number of pages for this query.
#'   `invalid_entities`, if present, is a data.frame listing queried entity
#'   identifiers that ClassyFire found invalid. `entities`, if present, is a
#'   nested data.frame giving classifications. Use function
#'   \code{\link{parse_classified_entities}} to parse these results into a data.frame
#'   of classified entities, suitable for use with the tree visualization,
#'   similarity analysis, or similarity visualization functions. Or call
#'   \code{\link{classify_structures}} which is a wrapper for this function and
#'   \code{\link{parse_classified_entities}}.
#' @export
#'
 query_classyfire <- function(input = NULL,
                              url = NULL,
                              label = "query",
                              type = "STRUCTURE",
                              retry_get_times = 10,
                              wait_min = 5,
                              retry_query_times = 3,
                              processing_wait_per_input = NULL){

   base_url <- "http://classyfire.wishartlab.com/queries"
   #if both input and url are NULL, throw error
   if(is.null(input) &
      is.null(url)){
     stop("Either `input` or `url` must be provided, but both were NULL.")
   }

   if(!is.null(input) &
      !is.null(url)){
     stop("Either `input` or `url` must be provided, but not both. You provided both.")
   }

 #if input provided, construct new query
   if(!is.null(input)){
     if(is.null(processing_wait_per_input)){
       processing_wait_per_input <- wait_min/length(input)
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
   url_good <- grepl(x = url,
                     pattern = paste0(base_url, "/", "\\d+", ".json"))
   if(!url_good){
     stop(paste0("URL is not of valid format. ",
                 "URL is ", url,
                 " and expected format is ",
                 base_url,
                 "/", "NNNNN.json",
                 " where NNNNN are one or more digits 0-9 denoting a query ID number."))
   }

   #try getting results for query
   resp <- httr::RETRY(verb = "GET",
                       url = url,
                       encode = "json",
                       times = retry_get_times,
                       pause_min = wait_min,
                       terminate_on = c(404))

   json_res <- httr::content(resp, "text")
   #json_res will be NULL if no content was returned
   if(!is.null(json_res)){
     json_parse <- jsonlite::fromJSON(json_res)
   }else{
     json_parse <- NULL
   }

   #initialize a placeholder json_parse in case everything else fails
   m <- regexec(text=url, pattern = "(\\d+)\\.json")
   query_id <- regmatches(x = url, m = m)[[1]][2]
   json_parse_default <- list("id" = query_id,
        "url" = url,
        "status_code" = resp$status_code,
        "label" = label,
        "classification_status" = "No classification",
        "number_of_elements" = NA_real_,
        "number_of_pages" = 1,
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
                                "classification_version" = character(0))
   )
   if(!is.null(input)){
     json_parse_default$number_of_elements <- length(input)
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

   #retry until Done or until max number of retries is reached
     retry_count <- 0
     while(!(json_parse$classification_status %in% "Done") &
           retry_count < retry_query_times){
       if(json_parse$classification_status %in% "In Queue"){
         #wait for query to come out of queue
         #this does not depend on size of input
         #calculate exponential backoff wait time
         wait_time <- runif(1, wait_min, wait_min*2^(retry_count))
         message(paste0("Query status is In Queue",
                        "; waiting ",
                        round(wait_time, digits = 2),
                        " seconds and retrying"))
       }else{ #query is processing
         #processing time *does* depend on size of input
         #calculate exponential backoff wait time
         proc_eff <- ceiling(runif(1,
                                   wait_min/length(input),
                                   processing_wait_per_input * 2^(retry_count)
                                   )
                             )
         wait_time <- ceiling(max(wait_min,  proc_eff*length(input)))
         message(paste0("Query status is ",
                        json_parse$classification_status,
                        "; waiting ",
                        round(wait_time, digits = 2),
                        " seconds (greater of ",
                        wait_min,
                        " seconds, or ",
                        proc_eff,
                        " seconds per input, with ",
                        length(input),
                        " inputs) ",
                        "and retrying"))
       } #end if/else to check classification status

       #wait and retry
       Sys.sleep(wait_time)
       resp <- httr::RETRY(verb = "GET",
                           url = url,
                           encode = "json",
                           times = retry_get_times,
                           pause_min = wait_min,
                           terminate_on = c(404))
       json_res <- httr::content(resp, "text")
       #json_res will be NULL if no content was returned
       if(!is.null(json_res)){
         json_parse <- jsonlite::fromJSON(json_res)
       }else{
         json_parse <- NULL
       }
       json_parse_default$status_code <- resp$status_code
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
       retry_count <- retry_count + 1
     } #end while loop
     #If classification was done, then results will be paginated
     #Each page = 10 entities
     #Page URLs are constructed like http://classyfire.wishartlab.com/queries/XXX.json?page=1

     #check if this URL ends in a page specification or not
     #if so, return the parsed JSON for this page
     url_onepage <- grepl(pattern = "?page=\\d+$",
                          x = url)
     if(url_onepage){
     output <- json_parse
     }else{ #if this is not a specific page, but an overall query,
       #recursively call this function to get results for each individual page
       #get vector of individual page URLs
       url_pages <- paste0(url,
                           "?page=",
                           seq(from = 1,
                               to = json_parse$number_of_pages,
                               by = 1))
       output <- lapply(url_pages,
                        function(url_pg){
                          query_classyfire(input = NULL,
                                           url = url_pg,
                                           label = label,
                                           type = type,
                                           retry_get_times = retry_get_times,
                                           wait_min = wait_min,
                                           retry_query_times = retry_query_times,
                                           processing_wait_per_input = processing_wait_per_input)
                        }
       )
     }
       return(output)
     }


#' Parse classified entities
#'
#' @param entities JSON element of classified entities as returned by ClassyFire
#'   API query
#' @param tax_level_labels ChemOnt taxonomy level labels (default
#'   \code{\link{chemont_tax_levels}})
# @param
#' @return A data.frame consisting of rows corresponding to each classified
#' entry from the `entities` input. A data.frame with one row for each entity identifier, and variables
#'  `identifier`, `smiles`, `inchikey`, `kingdom`, `superclass`, `class`,
#'   `subclass`, `level5`, `level6`, ... `level11`, `classification_version`,
#'    and `report`.
parse_classified_entities <- function(entities,
                                      tax_level_labels = chemont_tax_levels){
  identifier <- NULL
  level <- NULL
  name <- NULL
  #check to see whether classifications actually exist for these
  if(length(entities)==0){
    #if no classifications, return empty data.frame,
    #but with the expected variable names
    classified_entities <- data.frame(identifier = character(0),
                                      smiles = character(0),
                                      inchikey = character(0),
                                      #kingdom = character(0),
                                      #superclass = character(0),
                                      #class = character(0),
                                      #subclass = character(0),
                                      #level5 = character(0),
                                      #level6 = character(0),
                                      #level7 = character(0),
                                      #level8 = character(0),
                                      #level9 = character(0),
                                      #level11 = character(0),
                                      classification_version = character(0),
                                      level = character(0),
                                      name = character(0),
                                      report = character(0))
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

    cf <- as.data.frame(entities[c("identifier",
                         "smiles",
                         "inchikey",
                         "classification_version")])

    #kingdom, superclass, class, subclass are all data.frames
    #with one row for each identifier,
    #and variables name, description, chemont_id, url
    #giving the relevant ChemOnt taxonomy label for each identifier
    #(i.e. "name" in the "kingdom" element gives the kingdom for each identifier,
    #"name" in the "superclass" element gives the superclass for each identifier,
    #etc.)
    #rowbind these

    #note: if no classification for one of these four levels, the item may be NULL.
    #handle accordingly.
    for(item in c("kingdom",
                  "superclass",
                  "class",
                  "subclass")){
      #check if NULL
      if(is.null(entities[[item]])){
        #replace with list with NAs
        entities[[item]] <- list(name = NA_character_,
                                 description = NA_character_,
                                 chemont_id = NA_character_,
                                 url = NA_character_)
      }

      #check if character(0) and replace with NA
      for (subitem in names(entities[[item]])){
        if(length(entities[[item]][[subitem]])==0){
          entities[[item]][[subitem]] <- NA_character_
        }
      }
    }

    cf_class1 <- dplyr::bind_rows(as.list(entities[c("kingdom",
                                                       "superclass",
                                                       "class",
                                                       "subclass")]))
    #add a column for the identifiers -- repeat for each taxonomy level
    cf_class1$identifier <- rep(entities$identifier,
                                4)

    #"intermediate_nodes"
    #is a list with one element for each item in "input",
    #a 1-row data.frame with variables name, description, chemont_id, url
    #giving more-specific levels of classification, if any.
    #if no intermediate nodes for an entity,
    #its corresponding element in "intermediate_nodes" is an empty data.frame.
    #if it is an empty list (as happens for single-InChiKey queries without intermediate nodes)
    #then handle accordingly
    cf_class2 <- entities$intermediate_nodes
    if(length(cf_class2)>0){
    names(cf_class2) <- entities$identifier
    cf_class2_df <- dplyr::bind_rows(cf_class2, .id = "identifier")
    }else{
    cf_class2_df <- data.frame()
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
    direct$identifier <- entities$identifier
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
                                 by = c("name" = "Name")) %>%
      dplyr::select(identifier, level, name) %>%
      dplyr::distinct() #keep only unique rows

    #Change the levels from numbers to taxonomy level names
    #e.g. 1, 2, 3, 4 -> kingdom, superclass, class, subclass
    cf_class <- cf_class %>%
      dplyr::mutate(level = tax_level_labels[level])

    #We have found at least one instance where the ChemOnt tree listed something at the wrong level
    #(see data-raw/make_chemont_taxonomy.R).
    #Though we manually fixed that instance, there may be others.
    #If so, cf_class will have more than one "name" at a given "level."
    #This will cause `tidyr::pivot_wider()` to return a list-type column,
    #which then breaks downstream use of `dplyr::bind_rows()`.
    #I don't know how to handle this situation automatically.
    #I think we just have to continue returning the long format data.frame,
    #and make the user deal with it as they see fit.

    # #reshape to wide format, one column for each level
    # cf_class_wide <- cf_class %>%
    #   tidyr::pivot_wider(id_cols = "identifier",
    #                      names_from = "level",
    #                      values_from = "name")
    #
    # #add NA columns for any taxonomy levels not assigned for these entities
    # cols_add <- setdiff(tax_level_labels,
    #                     names(cf_class_wide))
    # classified_entities <-  cf_class_wide %>% as.data.frame()
    # if(!is.na(cols_add) &
    #    length(cols_add > 0)){
    # classified_entities[cols_add] <- NA_character_
    # }


    classified_entities <- cf_class
    #add a "report" column
    if(!"report" %in% names(entities)){
    classified_entities$report <- "ClassyFire returned a classification"
    }else{
      classified_entities$report <- entities$report
    }

    #merge in smiles and inchikey
    classified_entities <- dplyr::left_join(cf,
                                            classified_entities,
                                            by = "identifier") %>%
      as.data.frame()
  } #end if length(entities)>0

  return(classified_entities)
}


#' Query ClassyFire InChIKey
#'
#' This is a helper function that is used to communicate with the ClassyFire API
#' in order to retrieve classifications.
#'
#' @param inchikey A character string encoding the InChIKey of the chemical
#'   under query.
#' @param retry_get_times The number of times to retry the query, a positive
#'   integer with default value 3.
#' @param wait_min The number of seconds to pause in between attempts, used in
#'   the \code{\link{httr}{RETRY}} function.
#'
#' @return A JSON file with the classification data.
query_classyfire_inchikey <- function(inchikey,
                                      retry_get_times = 3,
                                      wait_min = 0.5){

  base_url <- "http://classyfire.wishartlab.com/entities"
  url <- paste0(base_url,
                "/",
                inchikey,
                ".json")
  #try getting results for query
  resp <- httr::RETRY(verb = "GET",
                      url = url,
                      encode = "json",
                      times = retry_get_times,
                      pause_base = wait_min,
                      pause_min = wait_min,
                      terminate_on = c(404))

  json_res <- httr::content(resp, "text")
  #json_res will be NULL if no content was returned
  if(!is.null(json_res)){
    json_parse <- jsonlite::fromJSON(json_res, simplifyDataFrame = FALSE)
  }else{
    json_parse <- NULL
  }

  #initialize a placeholder json_parse in case everything else fails
  #the items that are named lists are usually actually data.frames,
  #but lists are coerceable to data.frames and making it a list
  #helps with making it more general.
  json_parse_default <- list("identifier" = character(0),
                                   "smiles" = character(0),
                                   "inchikey" = character(0),
                                   "kingdom" = list("name" = character(0),
                                                          "description" = character(0),
                                                          "chemont_id" = character(0),
                                                          "url" = character(0)),
                                   "superclass" = list("name" = character(0),
                                                             "description" = character(0),
                                                             "chemont_id" = character(0),
                                                             "url" = character(0)),
                                   "class" = list("name" = character(0),
                                                        "description" = character(0),
                                                        "chemont_id" = character(0),
                                                        "url" = character(0)),
                                   "subclass" = list("name" = character(0),
                                                           "description" = character(0),
                                                           "chemont_id" = character(0),
                                                           "url" = character(0)),
                                   "intermediate_nodes" = list(),
                                   "direct_parent" = list("name" = character(0),
                                                                "description" = character(0),
                                                                "chemont_id" = character(0),
                                                                "url" = character(0)),
                                   "classification_version" = character(0))

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

  #assign identifier as inchikey
  json_parse$identifier <- inchikey


  return(json_parse)
}
