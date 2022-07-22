#' Add terminal label
#'
#' This function takes in a data.table of chemicals and their classifications
#' and determines the terminal classification label for each chemical.
#'
#' @param data A data.table of classified chemicals.
#' @param tip An alternate parameter controlling whether to include only tip
#'   labels.
#' @param labels An alternate parameter giving a list of the classification
#'   levels when not using ClassyFire classified data.
#' @param tree An alternate parameter giving a different tree structure if not
#'   using ChemOnt taxonomy.
#' @return A new data.table, augmenting the input data.table with a column
#'   consisting of the terminal labels.
#' @export
#' @import data.table
add_terminal_label <- function(data,
                               tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                                    'level5', 'level6', 'level7', 'level8',
                                                    'level9', 'level10', 'level11')){
  labels <- tidyr::pivot_longer(data, #reshape to longer format
                                cols = tidyselect::all_of(tax_level_labels),
                                names_to = "tax_level",
                                values_to = "label") %>%
    dplyr::filter(!is.na(label)) %>% #remove any unused levels
    dplyr::group_by( #group by item (e.g. chemical)
      dplyr::across(
        dplyr::all_of(setdiff(names(data),
                              tax_level_labels)
        )
      )
    ) %>%
    dplyr::slice_tail() %>%  #take most-specific label for each item
    dplyr::rename(terminal_label = label,
                  terminal_tax_level = tax_level)

  data <- merge(data, labels, by = setdiff(names(data),
                                           tax_level_labels))

  return(data)
}

calc_number_overlap <- function(base_data,
                                compare_data){
  #get terminal labels if not already there

}

#' Terminal label function
#'
#' This is a helper function to return the terminal label for
#' \code{\link{add_terminal_label}}. It returns the terminal label for an
#' input data.table consisting of a single row, a chemical with its
#' classification.
#'
#' @param t A data.table with classification data.
#' @param tip Alternate parameter for determining whether to only return tip
#'   labels.
#' @param tax_level_labels Parameter giving classification levels.
#' @param tree Alternate parameter for giving a taxonomy tree different from
#'   ChemOnt.
#' @return A string giving the terminal label.
#' @import data.table
terminal_function <- function(t, tip = FALSE, tax_level_labels, tree = NULL){
  labels <- t[1, .SD, .SDcols = which(names(t) %in% tax_level_labels)]

  if (is.null(tree)){
    tree <- chemont_tree
  }
  if (tip){
    if (any(labels %in% tree$tip.label)){
      return(labels[[which(labels %in% tree$tip.label)]])
    }
    return(NA_character_)
  }


  index <- which(unname(sapply(labels, function(t) {return(is.na(t) | t == '')})))

  if(length(index) == 0){
    return(labels[length(tax_level_labels)])
  }
  if(length(index) == length(tax_level_labels)){
    return(NA_character_)
  }

  return(labels[[(index[[1]] - 1)]])
}

#' Label level
#'
#' This function takes in a data.table of chemicals with classification data and
#' a taxonomy level, and returns all the unique labels for the given taxonomy
#' level.
#'
#' @param data A data.table with data that has been classified by some taxonomy.
#' @param level_label A string indicating a taxonomy level of the classified
#'   data.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @return The unique labels corresponding to the given level and classified
#'   data.
#' @import data.table
get_label_level <- function(data, level_label, tax_level_labels = NULL){
  if (is.null(tax_level_labels)){
    tax_level_labels <- c('kingdom', 'superclass', 'class', 'subclass',
                          'level5', 'level6', 'level7', 'level8',
                          'level9', 'level10', 'level11')
  }
  if (!(level_label %in% names(data) | !(level_label %in% tax_level_labels)))
    stop(paste('Please input a valid label!', level_label))

  # Collect the labels
  labels <- data[, sapply(.SD, unique), .SDcols = c(level_label)]

  # Remove NA from list of labels
  labels <- labels[!is.na(labels)]

  # Remove '' from list of labels
  labels <- labels[which(labels != '')]

  return(labels)
}


#' Get labels
#'
#' This is a helper function for retrieving labels in a data.frame of classified
#' chemicals.
#'
#' @param data A data.frame consisting of classified items. Rows are items; columns must include all names in \code{tax_level_labels}.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @return A list of terminal classification labels for each item in the input data.
get_labels <- function(data,
                       tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                            'level5', 'level6', 'level7', 'level8',
                                            'level9', 'level10', 'level11')){

  labels <- tidyr::pivot_longer(data, #reshape to longer format
                                cols = tidyselect::all_of(tax_level_labels),
                                names_to = "tax_level",
                                values_to = "label") %>%
    dplyr::filter(!is.na(label)) %>% #remove any unused levels
    dplyr::group_by( #group by item (e.g. chemical)
      dplyr::across(
        dplyr::all_of(setdiff(names(data),
                              tax_level_labels)
        )
      )
    ) %>%
    dplyr::slice_tail() %>% #take most-specific label for each item
    dplyr::pull(label) #return as a vector

  return(labels)
}

#' Number of labels
#'
#' This is a helper function to determine number of occurrences for each label
#' in a data.table containing chemicals and their classification data.
#'
#' @param data A data.table of chemicals with classifications.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not  using ClassyFire taxonomy
#' @return A named list of chemical labels and their number of occurrences.
get_number_of_labels <- function(data,
                                 tax_level_labels = c('kingdom', 'superclass',
                                                      'class', 'subclass',
                                                      'level5', 'level6',
                                                      'level7', 'level8',
                                                      'level9', 'level10',
                                                      'level11')){


  labels <- get_labels(data = data, tax_level_labels = tax_level_labels)

  N <- sum(sapply(labels, length))

  number_of_labels <- integer(N)
  names(number_of_labels) <- unlist(unname(labels))

  #print(names(number_of_labels))

  for (i in seq_along(tax_level_labels)){
    #print(paste('There are ', length(labels[[i]]), 'levels'))
    #print(labels[[i]])

    for (j in seq_along(labels[[i]])){
      index <- which(names(number_of_labels) == labels[[i]][[j]])
      #print(paste(tax_level_labels[[i]], labels[[i]][[j]]))
      #print(names(labels)[[i]])
      #print(labels[[i]][[j]])
      #print(index)
      #print(data[, .(names(labels)[[i]])])
      #print(data[, .SD, .SDcols = c(names(labels)[[i]])])
      number_of_labels[[index]] <- length(which(unname(unlist(data[, .SD, .SDcols = c(names(labels)[[i]])])) == labels[[i]][[j]]))
      #print(length(which(unname(unlist(data[, .SD, .SDcols = c(names(labels)[[i]])])) == labels[[i]][[j]])))
    }
  }

  return(number_of_labels)
}

#' Label length
#'
#' This is a helper function that returns the number of labels per taxonomy
#' level for a given data.table of chemicals and their classification data.
#'
#' @param label_list A named list of labels corresponding to taxonomy levels.
#' @return The number of labels per taxonomy level.
get_label_length <- function(label_list){
  lengths <- sapply(label_list, length)
  lengths
}

