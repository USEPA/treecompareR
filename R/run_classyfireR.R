
#' InChIKey classification
#'
#' This function uses the ClassyFire API to classify chemicals from an input
#' data.table using the InChIKey chemical identifier.
#'
#' This function queries ClassyFire's lookup table of pre-classified InChiKeys
#'
#' @param inchikeys A vector of InCHiKeys to be classified
#' @param tax_level_labels By default, the list of taxonomy levels for ClassyFire:
#'   \code{kingdom, superclass, class, subclass, level5, ... level11}.
#' @return A data.frame with the same number of rows as the length of input
#'   argument \code{inchikeys}, and columns consisting of "INCHIKEY" (containing input
#'   argument \code{inchikeys}) and one column for each of the ClassyFire taxonomy levels
#'   given in the input argument \code{tax_level_labels}
#' @export
#' @references
#' \insertRef{djoumbou2016classyfire}{treecompareR}
#'
#'
#' @seealso \code{\link{classify_structure}}
#'

classify_inchikeys <- function(inchikeys){
  if (!requireNamespace("classyfireR", quietly = TRUE)) {
    warning(paste("Package \"classyfireR\" must be installed to use this function.\n", "Returning input data.table!"))
    return(datatable)
  }

  tax_level_labels = c("kingdom",
                       "superclass",
                       "class",
                       "subclass",
                       paste0("level", 5:11))

  INCHIKEYS <- unique(inchikeys) #save time by removing duplicates

  #get classifications as a list of data.tables, one for each INCHIKEY
  class_list <- sapply(
    INCHIKEYS,
    function(this_inchikey){
      #if inchikey is NA, blank, or not in valid InChiKey format,
      #do not query API; just return NA
      if(is.na(this_inchikey) || #if InChiKey is NA
         !nzchar(trimws(this_inchikey)) || #if InChiKey is blank
         !webchem::is.inchikey(this_inchikey, #if not valid InChiKey format
                               type = "format")
         ){
        outlist <- rep(list(NA_character_),
                       times = length(tax_level_labels))
        names(outlist) <- tax_level_labels
        output <- as.data.frame(outlist)
      }else{ #if inchikey is valid, query ClassyFire API
      cf <- classyfireR::get_classification(this_inchikey)
      #if no classification available, return NA for classification
      if (is.null(cf) ||
          !("Classification" %in% names(classyfireR::classification(cf))
          )
      ){
        #return data frame with columns for all taxonomy levels, but NAs
        outlist <- rep(list(NA_character_),
                       times = length(tax_level_labels))
        names(outlist) <- tax_level_labels
        output <- as.data.frame(outlist)
      }else{ #if classification is available

        #get classification
        cf_class <- classyfireR::classification(cf)
        #select relevant columns
        cf_class_select <- dplyr::select(cf_class, Level, Classification)
        #strip whitespace from "level 5", "level 6", etc.
        cf_class_select <- cf_class_select %>%
          dplyr::mutate(Level = gsub(pattern = "\\s",
                                     replacement = "",
                                     x = Level))
        #reshape to wider format -- one column for each level
        #and format it as a data table
        output <- as.data.frame(
          tidyr::pivot_wider(cf_class_select,
                             names_from = "Level",
                             values_from = "Classification")
        )
        #if missing any taxonomy levels, fill them with NAs
        missing_levels <- setdiff(tax_level_labels,
                                  names(output))
        if(length(missing_levels)>0){
          output[missing_levels] <- NA_character_
        }
        #reorder columns to match tax_level_labels
        output <- output[, tax_level_labels]
      }
      }
      return(output)
    }, #end function to apply to each INCHIKEY
    simplify = FALSE, #sapply should return a list, not a vector/matrix
    USE.NAMES = TRUE #sapply should name the list elements after the input INCHIKEY
  )

  #now rowbind the list of data.frames to get one big data.frame
  inchi_class <- plyr::rbind.fill(class_list)
  #add InChiKey identifier column
  inchi_class$INCHIKEY <- INCHIKEYS #same as names(class_list)

  #now order as for input -- re-inserting any duplicates
  new_class <- inchi_class[match(inchikeys, inchi_class$INCHIKEY),
                           c("INCHIKEY",
                             tax_level_labels)]

  return(new_class)
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
  classify_structures <- function (input,
                                   tax_level_labels = chemont_tax_levels,
                                   ...){
    if (!requireNamespace("classyfireR", quietly = TRUE)) {
      warning(paste("Package \"classyfireR\" must be installed to use this function.\n", "Returning input data.table!"))
      return(input)
    }

    if(is.null(names(input))){
      names(input) <- input
    }


    #Initialize output data.frame.
    #Structure of output data.frame: One row for each input structure
    #Variables: c("identifier", "structure", "smiles", "inchikey",
    #tax_level_labels, "report")

    #prepare a named list of NA columns, one for each taxonomy level
    cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
    names(cols_list) <- tax_level_labels
    #put these as columns of a data.frame, along with input structure.
    output <- do.call(data.frame,
                      c(list("identifier" = names(input),
                             "structure" = input),
                        cols_list))
    output$smiles <- NA_character_
    output$inchikey <- NA_character_
    output$report <- NA_character_

    #If ClassyFire query fails, we'll just return the above.

    #Now, query ClassyFire.
    json_parse <- do.call(query_classyfire,
                          args = c(list(input = input),
                                   ...)
    )
    #We expect the JSON output to have the following names
    expected_json <- c("id",
                       "label",
                       "classification_status",
                       "number_of_elements",
                       "number_of_pages",
                       "invalid_entities",
                       "entities")
    #If these names are not present, throw a warning and return NAs
    if(!(all(expected_json %in% names(json_parse)))){
      warning(paste("The JSON output was not in the expected format and cannot be further parsed."))
      output$report <- "The JSON output was not in the expected format and cannot be further parsed."
      return(output)
    }else{ #if there was a valid response from ClassyFire
      if(json_parse$classification_status %in% "Done"){
        #pull list of valid entities (i.e., those with classifications)
        classified <- json_parse$entities

        #check to see whether classifications actually exist for these
        if(length(classified)==0){
          #if no classifications, return empty data.frame,
          #but with the expected variable names
          classified_entities <- data.frame(identifier = character(0),
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
        }else{ #if length(classified)>0
          #expected named items in classified and their expected classes:
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
          cf <- classified[, c("identifier",
                               "smiles",
                               "inchikey")]

          #kingdom, superclass, class, subclass are all data.frames
          #with one row for each identifier,
          #and variables name, description, chemont_id, url
          #giving the relevant ChemOnt taxonomy label for each identifier
          #(i.e. "name" in the "kingdom" element gives the kingdom for each identifier,
          #"name" in the "superclass" element gives the superclass for each identifier,
          #etc.)
          #rowbind these
          cf_class1 <- dplyr::bind_rows(as.list(classified[c("kingdom",
                                                     "superclass",
                                                     "class",
                                                     "subclass")]))
          #add a column for the identifiers -- repeat for each taxonomy level
          cf_class1$identifier <- rep(classified$identifier,
                                      4)

         #"intermediate_nodes"
          #is a list with one element for each item in "input",
          #a 1-row data.frame with variables name, description, chemont_id, url
          #giving more-specific levels of classification, if any.
          #if no intermediate nodes for an entity,
          #its corresponding element in "intermediate_nodes" is an empty data.frame.
          cf_class2 <- classified$intermediate_nodes
          names(cf_class2) <- classified$identifier
          cf_class2_df <- dplyr::bind_rows(cf_class2, .id = "identifier")
          rm(cf_class2)

          #rowbind the kingdom-thru-superclass labels,
          #and the "intermediate nodes" labels.
          cf_class <- dplyr::bind_rows(cf_class1,
                                       cf_class2_df)

          #"classified" element 'direct_parent' is a data.frame
          #with one row for each identifier,
          #and variables name, description, chemont_id, url
          direct <- classified$direct_parent
          direct$identifier <- classified$identifier
          cf_class <- dplyr::bind_rows(cf_class, direct_final)

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

          #reshape to wide format, one column for each level
          cf_class_wide <- cf_class %>%
            tidyr::pivot_wider(id_cols = "identifier",
                               names_from = "level",
                               values_from = "name")

          #add NA columns for any taxonomy levels not assigned for these entities
          cols_add <- setdiff(tax_level_labels,
                              names(cf_class_wide))
          classified_entities <-  cf_class_wide %>% as.data.frame()
          classified_entities[cols_add] <- NA_character_

          #add a "report" column
          classified_entities$report <- "ClassyFire returned a classification"

          #merge in smiles and inchikey
          classified_entities <- dplyr::left_join(cf,
                                                  classified_entities,
                                                  on = "identifier")
        } #end if length(classified)>0

        #Add the rows for any invalid entities (without classifications)
        invalid <- json_parse$invalid_entities
        if(length(invalid)==0){
          #report empty data.frame but with the expected variables
          invalid <- data.frame(identifier = character(0),
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
          #Add identifier column (the names of the inputs)
          invalid$identifier <- names(input)[input %in% invalid$structure]

          #Add NA columns for taxonomy levels, to denote no classification
          invalid[tax_level_labels] <- NA_character_
          invalid$smiles <- NA_character_
          invalid$inchikey <- NA_character_
          invalid$report <- "ClassyFire reports invalid entity"
        } # if(length(invalid)>0)

        #Start constructing output:
        #entities handled by ClassyFire (either classified or reported invalid)
        cf_entities <- dplyr::bind_rows(classified_entities,
                                   invalid)
        #add structure column
        cf_entities$structure <- input[cf_entities$identifier]

        #Get any entities not handled by ClassyFire
        missing_entities <- dplyr::anti_join(output,
                                             cf_entities,
                                             by = "identifier")

        #Output: bind classified, invalid, and missing entities
        output <- dplyr::bind_rows(cf_entities,
                                   missing_entities) %>%
          as.data.frame() #convert from tibble to data.frame

        #set order the same as input
        rownames(output) <- NULL
        output <- output[match(names(input), output$identifier), ]
      } #end if(json_parse$classification_status %in% "Done")
    } #end if(length(json_parse)>0)

    rownames(output) <- NULL
    output <- output[, c("identifier",
                         "structure",
                         "smiles",
                         "inchikey",
                         tax_level_labels,
                         "report")]
    return(output)
  }

#' Send ClassyFire query
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
#' @param label An optional text label for the ClassyFire query. Default
#'   \code{"query"}.
#' @param type String giving ClassyFire query type. Default \code{"STRUCTURE"}.
#' @param retry_get_times Max number of times to retry GET command to ClassyFire
#'   API to get status of query. Default 10.
#' @param wait_min Minimum time to wait before retrying any GET command or
#'   query. Default 5 seconds (because ClassyFire requests that POST requests be
#'   limited to 12 per second). Highly recommend not changing this to be less
#'   than 5 seconds, to be respectful of the ClassyFire server. All wait
#'   times use exponential backoff with full jitter (\code{wait time = runif(1, wait_min, wait_min *
#'   2^(attempt number)}).
#' @param retry_query_times Max number of attempts to retry an "In Queue" or
#'   "Processing" query (with wait time in between tries). Default 3.
#' @param processing_wait_per_input Minimum number of seconds to wait per input
#'   string before retrying to retrieve results for a query whose status is
#'   "Processing." Default NULL results in \code{wait_min/length(input)}. Uses exponential
#'   backoff: effective wait time per input = \code{processing_wait_per_input *
#'   2^(attempt number)}. Total wait time is either the effective wait time per
#'   input times the number of input structures, or \code{wait_min} seconds,
#'   whichever is greater.
#' @return A named list object: the JSON output from the ClassyFire query, parsed into list form.
#' @export
#'
 query_classyfire <- function(input,
                              label = "query",
                              type = "STRUCTURE",
                              retry_get_times = 10,
                              wait_min = 5,
                              retry_query_times = 3,
                              processing_wait_per_input = NULL){

   if(is.null(processing_wait_per_input)){
     processing_wait_per_input <- wait_min/length(input)
   }

   base_url <- "http://classyfire.wishartlab.com/queries"
   #create tab-separated input
   query_input <- paste(names(input),
                        input,
                        sep = "\t",
                        collapse = "\n")
   #convert to JSON
   q <- rjson::toJSON(list(label = label,
                           query_input = query_input,
                           query_type = "STRUCTURE")
                      )
   #construct POST request
   resp <- httr::POST(url = base_url,
                      body = q,
                      httr::content_type_json(),
                      httr::accept_json(),
                      httr::timeout(getOption("timeout"))
                      )
   post_cont <- httr::content(resp)
   url <- paste0(base_url, "/", post_cont$id, ".json")
   resp <- httr::RETRY(verb = "GET",
                       url = url,
                       encode = "json",
                       times = retry_get_times,
                       pause_min = wait_min,
                       terminate_on = c(404))

   if (!(resp$status_code %in% 200)){
     message(paste0("ClassyFire API response status code: ",
                   resp$status_code),
             ". Returning empty list.")
     #return an empty list
     json_parse <- list()
   }else{
     #get and parse the JSON response
     json_res <- httr::content(resp, "text")
     json_parse <- jsonlite::fromJSON(json_res)

     #check classification status

     #if in queue, wait and retry, up to retry_query_times
     retry_count <- 0
     while(!(json_parse$classification_status %in% c("Done")) &
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
                        proc_eff, digits = 2,
                        " seconds per input, with ",
                        length(input),
                        " inputs) ",
                        "and retrying"))
       } #end if/else to check classification status

       #wait and retry
       Sys.sleep(wait_time)
       resp <- httr::GET(url = url, httr::accept_json())
       if (resp$status_code == 200){
       json_res <- httr::content(resp, "text")
       json_parse <- jsonlite::fromJSON(json_res)
       }
       retry_count <- retry_count + 1
     } #end while loop


     message(paste("Classification status", json_parse$classification_status))
   }

     return(json_parse)
 }
