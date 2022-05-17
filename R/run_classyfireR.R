
# This function takes in a data.table with columns "PREFERRED_NAME", "CASRN", and "INCHIKEY"
# and runs the the classyFireR API call for each chemical in the data.table.

#' Classify chemicals from input data.table using InChIKey
#'
#' @param datatable A data.table with columns 'PREFERRED_NAME', 'CASRN', 'INCHIKEY
#' @return A data.table object with classifications attached to each row of input data.table
#' @export
#' @import data.table
#' @importFrom purrr map
#' @import classyfireR
classify_datatable <- function(datatable){
  INCHIKEY <- NULL
  if (!data.table::is.data.table(datatable)){
    stop('Input must be a data.table object!')
  }
  col_names <- c("PREFERRED_NAME", "CASRN", "INCHIKEY")
  present_names <- which(col_names %in% names(datatable))
  if(length(present_names) < 3){
    missing_names <- col_names[setdiff(1:3, present_names)]
    stop('Input data.table missing the col(s) \n', paste(missing_names, collapse = '\n'), '!')
  }

  # Copy table
  new_table <- copy(datatable)

  # Get unique INCHIKEY values and remove NA values
  INCHIKEYS <- new_table[, unique(INCHIKEY)]
  INCHIKEYS <- INCHIKEYS[!is.na(INCHIKEYS)]


  # Run classyFireR on unique INCHIKEY values
  classifications <- purrr::map(INCHIKEYS, classyfireR::get_classification)

  # Populate data.table with Classification data
  for (i in seq_along(INCHIKEYS)){
    # Handle case where INCHIKEY submitted returned no classification or
    # the classification returned was empty but not null
    if (is.null(classifications[[i]]) || !("Classification" %in% names(classifications[[i]]@classification))){
      classifiers <- rep('', 11)
    } else {
      temp <- classyfireR::classification(classifications[[i]])$Classification
      classifiers <- c(temp, rep('', 11 - length(temp)))
    }
    new_table[INCHIKEY == INCHIKEYS[[i]], c("kingdom", "superclass", "class", "subclass", "level5", "level6", "level7", "level8", "level9", "level10", "level11") := as.list(classifiers)]
  }
  # Set rows with NA in INCHIKEY to empty classification
  new_table[is.na(INCHIKEY), c("kingdom", "superclass", "class", "subclass", "level5", "level6", "level7", "level8", "level9", "level10", "level11") := as.list(rep('',11))]

  return(new_table)
}

# This function completes a second pass after the classify_datatable function by
# looking at the rows missing a classification and attempting to classify
# via the associated SMILES string if present

#' Classify chemicals from input data.table using SMILES
#'
#' @param datatable A data.table that is the output from classify_datatable.
#' @return A data.table object with classification information for each row.
#' @export
#' @import data.table
#' @importFrom purrr map2
#' @import classyfireR
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
  SMILES_str <- SMILES_str[!is.na(SMILES_str)]
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
                       message(e)
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
