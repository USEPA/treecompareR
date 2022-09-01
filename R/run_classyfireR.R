
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
#' @param tax_level_labels By default, the list of taxonomy levels for
#'   ClassyFire: \code{kingdom, superclass, class, subclass, level5, ...
#'   level11}.
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
                                   label = "query",
                                   type = "STRUCTURE",
                                   retry_get_times = 10,
                                   wait_min = 5,
                                   retry_query_times = 3,
                                   processing_wait_per_input = NULL,
                                   tax_level_labels = chemont_tax_levels){
    if (!requireNamespace("classyfireR", quietly = TRUE)) {
      warning(paste("Package \"classyfireR\" must be installed to use this function.\n", "Returning input data.table!"))
      return(input)
    }

    if(is.null(names(input))){
      names(input) <- input
    }

    if(is.null(processing_wait_per_input)){
      processing_wait_per_input <- wait_min/length(input)
    }

    base_url <- "http://classyfire.wishartlab.com/queries"
    query_input <- paste(names(input), input, sep = "\t", collapse = "\n")
    q <- rjson::toJSON(list(label = label, query_input = query_input,
                            query_type = "STRUCTURE"))
    resp <- httr::POST(url = base_url, body = q, httr::content_type_json(),
                       httr::accept_json(), httr::timeout(getOption("timeout")))
    post_cont <- httr::content(resp)
    url <- paste0(base_url, "/", post_cont$id, ".json")
    resp <- httr::RETRY("GET",
                        url = url,
                        encode = "json",
                        times = retry_get_times,
                        pause_min = wait_min)
    print(resp)
    if (resp$status_code == 200) {
      json_res <- classyfireR::get_query(query_id = post_cont$id, format = "json")
      json_parse <- jsonlite::fromJSON(json_res)
      #check classification status
      #if in queue, wait and retry, up to retry_query_times
      retry_count <- 0
      while(!(json_parse$classification_status %in% c("Done",
                                                      "In progress")) &
            retry_count < retry_query_times){
        if(json_parse$classification_status %in% "In Queue"){
          #exponential backoff
          wait_time <- runif(1, wait_min, wait_min*2^(retry_count))
          message(paste0("Query status is In Queue",
                         "; waiting ",
                         wait_time,
                         " seconds and retrying"))
        }else if(json_parse$classification_status %in% "Processing"){
          proc_eff <- ceiling(runif(1, wait_min/length(input), processing_wait_per_input * 2^(retry_count)))
          wait_time <- ceiling(max(wait_min,  proc_eff*length(input)))
          message(paste0("Query status is Processing",
                         "; waiting ",
                         wait_time,
                         " seconds (greater of ",
                         wait_min,
                         " seconds, or ",
                         proc_eff,
                         " seconds per input, with ",
                         length(input),
                         " inputs) ",
                         "and retrying"))
        }else{
          proc_eff <- ceiling(runif(1, wait_min/length(input), processing_wait_per_input * 2^(retry_count)))
          wait_time <- ceiling(max(wait_min,  proc_eff*length(input)))
          message(paste0("Query status is ",
                         json_parse$classification_status,
                         "; waiting ",
                         wait_time,
                         " seconds (greater of ",
                         wait_min,
                         " seconds, or ",
                         proc_eff,
                         " seconds per input, with ",
                         length(input),
                         " inputs) ",
                         "and retrying"))
        }


        #wait and retry
        Sys.sleep(wait_time)
        json_res <- classyfireR::get_query(query_id = post_cont$id, format = "json")
        json_parse <- jsonlite::fromJSON(json_res)
        retry_count <- retry_count + 1
      }


      message(paste("Classification status", json_parse$classification_status))

      if(json_parse$classification_status %in% c("Done",
                                                 "In progress")){

        #pull list of valid entities (i.e., those with classifications)
        classified <- json_parse$entities

        #check to see whether classifications actually exist for these
        if(all(c("identifier",
                 "smiles",
                 "inchikey") %in% names(classified)) &
           length(classified)>0){

          #if any valid entities, take relevant portions of classification data.frame
          #get identifiers, smiles, inchikey
          cf <- classified[, c("identifier",
                               "smiles",
                               "inchikey")]
          #loop over levels kingdom, superclass, class, subclass
          #and extract the corresponding labels from the corresponding list elements
          cf_class1 <- sapply(c("kingdom",
                                "superclass",
                                "class",
                                "subclass"),
                              function(x) {
                                tmp <- classified[, x][, c("name",
                                                           "chemont_id")]
                                tmp$identifier <- classified$identifier
                                tmp
                              },
                              simplify = FALSE,
                              USE.NAMES = TRUE) %>%
            dplyr::bind_rows(.id = "level")

          #add level number to cf_class1
          cf_class1 <- cf_class1 %>%
            dplyr::group_by(identifier) %>%
            dplyr::mutate(level_num = dplyr::row_number())

          #anything more specific than "subclass" will be in "intermediate nodes"
          #(or in "direct parent" but we'll handle that later)
          #which is a list with one element for each item in "input"
          #if no intermediate nodes, the element will be an empty list
          cf_class2 <- classified$intermediate_nodes
          names(cf_class2) <- classified$identifier
          # #for any empty data frames, give a name of NA
          # cf_class2 <- sapply(cf_class2,
          #                     function(x) if(length(x)==0){
          #                       data.frame(name = NA_character_,
          #                                  chemont_id = NA_character_)
          #                       }else{
          #                         x
          #                       },
          #                     simplify = FALSE,
          #                   USE.NAMES = TRUE)
          cf_class2_df <- dplyr::bind_rows(cf_class2, .id = "identifier")

          if(nrow(cf_class2_df) > 0){
            #now add level -- starting with "level 5" and going up as many rows as exist for each identifier
            cf_class2_df <- cf_class2_df %>%
              dplyr::group_by(identifier) %>%
              dplyr::mutate(level_num = dplyr::row_number() + 4,
                            level = paste0("level", level_num)) %>%
              dplyr::select(identifier, name, chemont_id, level, level_num)
          }


          #rowbind
          cf_class <- dplyr::bind_rows(cf_class1,
                                       cf_class2_df)

          #now take the direct ancestor of each one
          direct <- classified$direct_parent
          direct$identifier <- classified$identifier
          direct <- direct %>% dplyr::select(identifier, name, chemont_id)
          #if the name and chemont_id do not already appear in cf_class for each identifier,
          #i.e. if this ancestor is not already listed,
          #then add them as a final level
          direct_final <- dplyr::anti_join(direct, cf_class,
                                           by = c("identifier", "name", "chemont_id"))

          cf_class <- dplyr::bind_rows(cf_class, direct_final)

          #if final level was added, add the appropriate level number
          #which will be 1+ the max level number for this identifier
          cf_class <- cf_class %>%
            dplyr::mutate(level_num2 = dplyr::if_else(is.na(level_num),
                                                      max(level_num, na.rm = TRUE)+1,
                                                      as.numeric(level_num))) %>%
            dplyr::mutate(level_num = NULL) %>%
            dplyr::rename(level_num = level_num2) %>%
            dplyr::mutate(level = dplyr::if_else(is.na(level), #paste level number to create level label
                                                 paste0("level", level_num),
                                                 level))

          #reshape to wide format, one column for each level
          cf_class_wide <- cf_class %>%
            tidyr::pivot_wider(id_cols = c("identifier"),
                               names_from = "level",
                               values_from = "name")

          #add NA columns for any unused levels in tax_level_labels
          cols_add <- setdiff(tax_level_labels,
                              names(cf_class_wide))
          output <-  cf_class_wide %>% as.data.frame()
          output[cols_add] <- NA_character_

          #merge in the smiles and inchikeys
          output <- merge(cf, output, by = "identifier")

          #add a "report" column
          output$report <- "Classification returned"

          #add a "structure" column
          output$structure <- input[output$identifier]

          #Add the rows for any invalid entities (without classifications)
          invalid <- json_parse$invalid_entities
          if(length(invalid)>0){
            #Add identifier column (the names of the inputs)
            invalid$identifier <- names(input)[input %in% invalid$structure]

            #Add NA columns for taxonomy levels, to denote no classification
            invalid[tax_level_labels] <- NA_character_
            invalid$smiles <- NA_character_
            invalid$inchikey <- NA_character_
          }

          #rowbind
          output <- dplyr::bind_rows(output, invalid)

          #for any structures not otherwise handled, return NAs and return classification status as report
          cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
          names(cols_list) <- tax_level_labels
          missing <- do.call(data.frame,
                             c(list("identifier" = names(input[!(input %in% output$structure)]),
                                    "structure" = input[!(input %in% output$structure)]),
                               cols_list,
                               list("check.names" = FALSE)))
          missing$report <- paste0("No classification.")
          missing$smiles <- NA_character_
          missing$inchikey <- NA_character_

          output <- dplyr::bind_rows(output, missing)


          #set order the same as input
          rownames(output) <- NULL
          output <- output[match(names(input), output$identifier), ]
        }else{ #if length(classified)==0
          #this could mean either ClassyFire timed out, or no entities could be classified
          #return invalid entities if any
          invalid <- json_parse$invalid_entities
          if(length(invalid)>0){
            #Add identifier column (the names of the inputs)
            invalid$identifier <- names(input)[match(invalid$structure, input)]

            #Add NA columns for taxonomy levels, to denote no classification
            invalid[tax_level_labels] <- NA_character_
          }

          #for any structures not in invalid entities, return NAs and return classification status as report
          cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
          names(cols_list) <- tax_level_labels
          missing <- do.call(data.frame,
                             c(list("identifier" = names(input[!(input %in% invalid$structure)]),
                                    "structure" = input[!(input %in% invalid$structure)]),
                               cols_list,
                               list("check.names" = FALSE)))
          missing$report <- "No classification."
          output <- dplyr::bind_rows(invalid, missing)

          output$smiles <- NA_character_
          output$inchikey <- NA_character_

          #set order the same as input
          rownames(output) <- NULL
          output <- output[match(names(input), output$identifier), ]

        }
      }else{ #if not (json_parse$classification_status %in% "Done")
        #return a data frame with all the same column names, but NAs for labels
        cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
        names(cols_list) <- tax_level_labels
        output <- do.call(data.frame,
                          c(list("identifier" = names(input),
                                 "structure" = input),
                            cols_list))
        output$smiles <- NA_character_
        output$inchikey <- NA_character_

      }

      rownames(output) <- NULL
      output <- output[, c("identifier",
                           "structure",
                           "smiles",
                           "inchikey",
                           tax_level_labels,
                           "report")]
    }else{ #if resp$status_code != 200
      #return a data frame with all the same column names, but NAs for labels
      cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
      names(cols_list) <- tax_level_labels
      output <- do.call(data.frame,
                        c(list("identifier" = names(input),
                               "structure" = input),
                          cols_list))
      output$smiles <- NA_character_
      output$inchikey <- NA_character_

    }

    rownames(output) <- NULL
    output <- output[, c("identifier",
                         "structure",
                         "smiles",
                         "inchikey",
                         tax_level_labels,
                         "report")]
    return(output)
  }
