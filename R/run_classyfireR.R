
#' InChIKey classification
#'
#' This function uses the ClassyFire API to classify chemicals from an input
#' data.table using the InChIKey chemical identifier.
#'
#' @param inchikeys A vector of InCHiKeys to be classified
#' @return A data.frame with the same number of rows as the length of input
#'   argument \code{inchikeys}, and columns consisting of "INCHIKEY" (containing input
#'   argument \code{inchikeys}) and one column for each of the ClassyFire taxonomy levels
#'   given in the input argument \code{tax_level_labels}
#' @export
#'
#' @seealso \code{\link{classify_structure}}
#'
classify_inchikeys <- function(inchikeys,
                               tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                                    'level 5', 'level 6', 'level 7', 'level 8',
                                                    'level 9', 'level 10', 'level 11')){

  INCHIKEYS <- unique(inchikeys) #save time by removing duplicates

  #get classifications as a list of data.tables, one for each INCHIKEY
  class_list <- sapply(
    INCHIKEYS,
    function(this_inchikey){
      #if inchikey is NA or not in valid InChiKey format,
      #do not query API; just return NA
      if(is.na(this_inchikey) ||
         !webchem::is.inchikey(this_inchikey, type = "format")
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
#' InChiKeys) and queries the ClassyFire API to get classifications for each
#' one.
#'
#' @param input A character vector of SMILES strings or InChiKeys. May
#'   optionally be named. If so, names will be returned as "identifier" column
#'   in output. If not named, the structural identifiers themselves will be
#'   returned as "identifier" column in output.
#' @param label An optional text label for the query. Default "query".
#' @param type String giving ClassyFire query type. Default "STRUCTURE".
#' @param retry_times Number of times to retry GET command. Default 100.
#' @param queued_wait Number of seconds to wait before retrying to retrieve
#'   results for a queued query. Default 2.
#' @param processing_wait_per_input Number of seconds to wait per input string
#'   before retrying to retrieve results for a query whose status is
#'   "Processing." Default 1.
#' @param retry_query_times Max number of attempts to retry an "In Queue" or
#'   "Processing" query (with wait time in between tries). Default 100.
#' @param tax_level_labels The default list of taxonomy levels for ClassyFire.
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
                           retry_times = 100,
                           retry_query_times = 100,
                           queued_wait = 2,
                           processing_wait_per_input = 1,
                           tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                                'level 5', 'level 6', 'level 7', 'level 8',
                                                'level 9', 'level 10', 'level 11'))
{
  if(is.null(names(input))){
    names(input) <- input
  }
  base_url <- "http://classyfire.wishartlab.com/queries"
  query_input <- paste(names(input), input, sep = "\t", collapse = "\n")
  q <- rjson::toJSON(list(label = label, query_input = query_input,
                          query_type = "STRUCTURE"))
  resp <- httr::POST(url = base_url, body = q, httr::content_type_json(),
                     httr::accept_json(), httr::timeout(getOption("timeout")))
  post_cont <- httr::content(resp)
  url <- paste0(base_url, "/", post_cont$id, ".json")
  resp <- httr::RETRY("GET", url = url, encode = "json", times = retry_times)
  print(resp)
  if (resp$status_code == 200) {
    json_res <- classyfireR::get_query(query_id = post_cont$id, format = "json")
    json_parse <- jsonlite::fromJSON(json_res)
    #check classification status
    #if in queue, wait 1 second and retry, up to retry_query_times
    retry_count <- 0
    while(!(json_parse$classification_status %in% "Done") &
          retry_count < retry_query_times){
      if(json_parse$classification_status %in% "In Queue"){
        wait_time <- queued_wait
        message(paste0("Query status is In Queue",
                      "; waiting ",
                      wait_time,
                      " seconds and retrying"))
      }else if(json_parse$classification_status %in% "Processing"){
        wait_time <- processing_wait_per_input * length(input)
        message(paste0("Query status is Processing",
                      "; waiting ",
                      wait_time,
                      "seconds (",
                      processing_wait_per_input,
                      " seconds per input, with ",
                      length(input),
                      " inputs) ",
                      "and retrying"))
      }else{
        wait_time <- 1
        message(paste0("Query status is ",
        json_parse$classification_status,
                       "; waiting ",
                       wait_time,
                       " seconds and retrying"))
      }


      #wait one second and retry
      Sys.sleep(wait_time)
      json_res <- classyfireR::get_query(query_id = post_cont$id, format = "json")
      json_parse <- jsonlite::fromJSON(json_res)
      retry_count <- retry_count + 1
    }


message(paste("Classification status", json_parse$classification_status))

   #pull list of valid entities (i.e., those with classifications)
  classified <- json_parse$entities
  #if any valid entities, take relevant portions of classification data.frame
  if(length(classified)>0){
    cf <- classified[, c("identifier",
                                  "smiles",
                                  "inchikey")]
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


    #anything more specific than "subclass" will be in "intermediate nodes"
    #which is a list with one element for each item in "input"
    cf_class2 <- classified$intermediate_nodes
    names(cf_class2) <- classified$identifier
    #for any empty data frames, give a name of NA
    cf_class2 <- sapply(cf_class2,
                        function(x) if(nrow(x)==0){
                          data.frame(name = NA_character_)
                          }else{
                            x
                          },
                        simplify = FALSE,
                      USE.NAMES = TRUE)
    cf_class2_df <- dplyr::bind_rows(cf_class2, .id = "identifier")

    #now add level -- starting with "level 5" and going up as many rows as exist for each identifier
cf_class2_df <- cf_class2_df %>%
  dplyr::group_by(identifier) %>%
  dplyr::mutate(level_num = dplyr::row_number() + 4,
                level = paste("level", level_num)) %>%
  dplyr::select(identifier, name, chemont_id, level, level_num)

#add level number to cf_class
cf_class1 <- cf_class1 %>%
  dplyr::group_by(identifier) %>%
  dplyr::mutate(level_num = dplyr::row_number())

#rowbind
cf_class <- dplyr::bind_rows(cf_class1,
                             cf_class2_df)

#now take the direct ancestor of each one
direct <- classified$direct_parent
direct$identifier <- classified$identifier
direct <- direct %>% dplyr::select(identifier, name, chemont_id)
#if the name and chemont_id do not already appear in cf_class for each identifier,
#then add them as a final level
direct_final <- dplyr::anti_join(direct, cf_class,
                                 by = c("identifier", "name", "chemont_id"))

cf_class <- dplyr::bind_rows(cf_class, direct_final)

#if final level was added, add the appropriate level number
#which will be 1+ the max level number for this identifier
cf_class <- cf_class %>%
  dplyr::mutate(level_num2 = if_else(is.na(level_num),
                              max(level_num, na.rm = TRUE)+1,
                              level_num)) %>%
  dplyr::mutate(level_num = NULL) %>%
  dplyr::rename(level_num = level_num2) %>%
  dplyr::mutate(level = if_else(is.na(level), #paste level number to create level label
                         paste("level", level_num),
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
#Add identifier column (the names of the inputs)
invalid$identifier <- names(input)[input %in% invalid$structure]

#Add NA columns for taxonomy levels, to denote no classification
invalid[tax_level_labels] <- NA_character_
invalid$smiles <- NA_character_
invalid$inchikey <- NA_character_

#rowbind
output <- dplyr::bind_rows(output, invalid)

#for any structures not otherwise handled, return NAs and return classification status as report
cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
names(cols_list) <- tax_level_labels
missing <- do.call(data.frame,
                  c(list("identifier" = names(input[!input %in% output$structure]),
                         "structure" = input[!input %in% output$structure]),
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

    #Add identifier column (the names of the inputs)
    invalid$identifier <- names(input)[match(invalid$structure, input)]

    #Add NA columns for taxonomy levels, to denote no classification
    invalid[tax_level_labels] <- NA_character_

    #for any structures not in invalid entities, return NAs and return classification status as report
cols_list <- as.list(rep(NA_character_, length(tax_level_labels)))
names(cols_list) <- tax_level_labels
missing <- do.call(data.frame,
                  c(list("identifier" = names(input[!input %in% invalid$structure]),
                         "structure" = input[!input %in% invalid$structure]),
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
