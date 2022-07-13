
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
#' @seealso \code{\link{classify_by_smiles}}
#'
classify_inchikeys <- function(inchikeys,
                               tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                                    'level5', 'level6', 'level7', 'level8',
                                                    'level9', 'level10', 'level11')){

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
      }else{ #if inchikey is not NA or blank, query ClassyFire API
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

# This function completes a second pass after the classify_datatable function by
# looking at the rows missing a classification and attempting to classify
# via the associated SMILES string if present

#' SMILES string classification
#'
#' This function takes a data.table of chemicals and classifications from
#' \code{\link{classify_datatable}} and attempts to classify chemicals from the
#' input data.table using SMILES strings that are missing classification data.
#'
#' @param datatable A data.table that is the output from classify_datatable.
#' @return A data.table object with classification information for each row.
#' @export
#' @import data.table
#' @importFrom purrr map2
#' @import classyfireR
#'
#' @seealso \code{\link{classify_datatable}}
#'
classify_by_smiles <- function(datatable){
  INCHIKEY <- NULL
  SMILES <- NULL
  kingdom <- NULL
  if (!data.table::is.data.table(datatable)){
    stop('Input must be a data.table object!')
  }
  col_names <- c("PREFERRED_NAME", "CASRN", "INCHIKEY", "kingdom", "superclass",
                 "class", "subclass", "level5", "level6", "level7", "level8",
                 "level9", "level10", "level11")
  present_names <- which(col_names %in% names(datatable))
  if(length(present_names) < 14){
    missing_names <- col_names[setdiff(1:14, present_names)]
    stop('Input data.table missing the col(s) \n', paste(missing_names, collapse = '\n'), '!')
  }

  # Copy data.table
  new_table <- copy(datatable)

  # Get unique SMILES strings and remove NA values, '' values
  SMILES_str <- new_table[is.na(INCHIKEY) | kingdom == '', unique(SMILES)]
  if (length(SMILES_str)==0) return(new_table)
  SMILES_str <- SMILES_str[!is.na(SMILES_str)]
  if (length(SMILES_str)==0) return(new_table)
  SMILES_str <- SMILES_str[sapply(SMILES_str, function(t) {t != ''})]

  if (length(SMILES_str)==0) return(new_table)

  # Set names for list of SMILES
  names(SMILES_str) <- paste0('MOL', 1:length(SMILES_str))

  classifications <- purrr::map2(SMILES_str, names(SMILES_str), function(s, n){
    data <- tryCatch({classyfireR::submit_query(label = n,
                                   input = s,
                                   structure)},
                     error = function(e) {
                       print(s)
                       print(e$message)
                       return(NULL)}
    )
  })

  # Populate data.table with Classification data
  for (i in seq_along(SMILES_str)){
    # Handle case where INCHIKEY submitted returned no classification or
    # the classification returned was empty but not null
    if (is.null(classifications[[i]]) || !("Classification" %in% names(classifications[[i]]@classification))){
      classifiers <- rep('', 11)
    } else {
      temp <- classyfireR::classification(classifications[[i]])$Classification


      classifiers <- c(temp, rep('', 11 - length(temp)))

    }
    new_table[SMILES == SMILES_str[[i]], c("kingdom", "superclass", "class", "subclass", "level5", "level6", "level7", "level8", "level9", "level10", "level11") := as.list(classifiers)]
  }
  return(new_table)
}
