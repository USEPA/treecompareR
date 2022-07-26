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
  data_orig <- copy(data) #save original input data

  #check that the input data.frame has been classified properly
  if(!any(tax_level_labels %in% names(data))){
    stop(paste("The input data.frame does not appear to be classified",
               "according to the taxonomy with levels defined in",
               "the input 'tax_level_labels' as",
               paste(tax_level_labels, collapse = ", "),
               "because the input data.frame does not contain any columns",
               "named for these taxonomy levels.",
               "Please check that the input data.frame is classified,",
               "and/or check that the input 'tax_level_labels'",
               "matches the taxonomy levels of the classified data.frame."))
  }

  #if data.frame already has terminal label data, throw a warning,
  #but proceed
  if("terminal_label" %in% names(data)){
    warning(paste("Column 'terminal_label' already exists",
    "in the input data.frame;",
    "it will be overwritten"))
    data[["terminal_label"]] <- NULL
    data_orig[["terminal_label"]] <- NULL
  }

  if("terminal_level" %in% names(data)){
    warning(paste("Column 'terminal_level' already exists",
                  "in the input data.frame;",
                  "it will be overwritten"))
    data[["terminal_level"]] <- NULL
    data_orig[["terminal_level"]] <- NULL
  }

  #if data.frame is missing one or more levels (but not all of them),
  #throw a warning and treat those levels as unused (i.e. all NA)
  if(!all(tax_level_labels %in% names(data))){
    missing_tax_levels <- setdiff(tax_level_labels,
                                  names(data))
    warning(paste("Input data.frame is missing columns for taxonomy levels",
            paste(missing_tax_levels, collapse = "; "),
            "These levels will be treated as though they were unused",
            "(i.e., as though those columns were present,",
            "but filled with NAs)."))
    #add the missing columns with NAs
    data[missing_tax_levels] <- rep(NA_character_, nrow(data))
  }

  #sort taxonomy level columns in order as given in tax_level_labels
  #this will ensure that most-specific (terminal) label comes last
  data <- data[c(setdiff(names(data), #all other columns come first
                         tax_level_labels),
                 tax_level_labels)]

  labels <- tidyr::pivot_longer(data, #reshape to longer format
                                cols = tidyselect::all_of(tax_level_labels),
                                names_to = "tax_level",
                                values_to = "label") %>%
    dplyr::filter(!is.na(label)) %>% #remove any unused levels
    dplyr::group_by( #group by item (e.g. chemical)
      dplyr::across( #assumed identified by everything *except* tax_level_labels
        dplyr::all_of(setdiff(names(data),
                              tax_level_labels)
        )
      )
    ) %>%
    dplyr::slice_tail() %>%  #take most-specific label for each item (i.e. last row)
    dplyr::rename(terminal_label = label, #rename cols to refer to "terminal"
                  terminal_tax_level = tax_level) %>%
    dplyr::mutate(terminal_level = match(terminal_tax_level,
                                         tax_level_labels)) %>%
    dplyr::mutate(terminal_tax_level = NULL)

  #merge terminal label & terminal level info back into original
  data_out <- merge(data_orig,
                    labels,
                    by = setdiff(names(data_orig),
                                           tax_level_labels))

  return(data_out)
}


#' calculate overlap between two classified datasets at the individual entity level
#'
#'@param data_1 A data.frame of classified entities to be used as the base set
#'@param compare_dat A data.frame of classified entitites to be used as hte comp
calc_number_overlap <- function(data_1,
                                data_2,
                                entity_id_col,
                                at_level = "terminal",
                                tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                                     'level5', 'level6', 'level7', 'level8',
                                                     'level9', 'level10', 'level11')){
  if(at_level %in% "terminal"){
  #get terminal labels if not already there
  if(!("terminal_label" %in% names(data_1))){
    data_1 <- add_terminal_label(data_1,
                                    tax_level_labels = tax_level_labels)
  }

  if(!("terminal_label" %in% names(data_2))){
    data_2 <- add_terminal_label(data_2,
                                    tax_level_labels = tax_level_labels)
  }
    group_col <- "terminal_label"
  }else if(is.numeric(at_level)){
    if(at_level > length(tax_level_labels)){
      stop(
        paste("Cannot find overlap at level  'at_level' =",
              at_level,
              "because it is greater than the max level",
              "defined by the length of 'tax_level_labels' =",
              paste(tax_level_labels, collapse = ", ")
        )
      )
    }else{
      #pull the corresponding taxonomy level label
      group_col <- tax_level_labels[at_level]
    }
  }else if(is.character(at_level)){
    if(!(at_level %in% tax_level_labels)){
      stop(
        paste("Cannot find overlap at level 'at_level' =",
              at_level,
              "because it is not one of the levels defined in",
              "'tax_level_labels' =",
              paste(tax_level_labels, collapse = ", ")
        )
      )
    }else{
    #interpret as an explicit taxonomy level label
    group_col <- at_level
    }
  }

  count_1 <- data_1 %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) %>%
    dplyr::summarise(n_1 = dplyr::n_distinct(
      dplyr::across(
        dplyr::all_of(entity_id_col)
        )
      )
      )

  count_2 <- data_2 %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) %>%
    dplyr::summarise(n_2 = dplyr::n_distinct(
      dplyr::across(
        dplyr::all_of(entity_id_col)
      )
    )
    )

  df_count <- merge(count_1,
                    count_2,
                    by = group_col,
                    all = TRUE)


  #replace NAs with zeros for labels with no entities in one data set
  df_count$n_1[is.na(df_count$n_1)] <- 0
  df_count$n_2[is.na(df_count$n_2)] <- 0

  get_overlap <-  function(grouplab,
                           group_col,
                           entity_id_col){
    id_1 <- data_1[data_1[[group_col]] %in% grouplab, entity_id_col]
    id_2 <- data_2[data_2[[group_col]] %in% grouplab, entity_id_col]
    n_intersect <- length(intersect(id_1, id_2))
    n_union <-  length(union(id_1, id_2))
    simil <- n_intersect/n_union
    data.frame("n_intersect" = n_intersect,
               "n_union" = n_union,
               "simil" = simil)
  }

  df_list <- sapply(df_count[[group_col]],
                    get_overlap,
                    group_col = group_col,
                    entity_id_col = entity_id_col,
                    simplify = FALSE,
                    USE.NAMES = TRUE)

  overlap_df <- dplyr::bind_rows(df_list,
                                 .id = group_col)

  outdf <- merge(df_count,
                 overlap_df,
                 by = group_col,
                 all = TRUE)

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

