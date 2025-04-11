#' Add terminal label
#'
#' Add terminal classification labels
#'
#' Adds the terminal classification label for each entity to a classified data
#' set.
#'
#' The terminal classifications will be added in a new variable
#' `terminal_level`. If a variable with that name already exists in the data,
#' then it will be overwritten (with a warning).
#'
#' @param data A `data.frame` of classified entities (or something that can be
#'   coerced to one using [base::as.data.frame()]). Should have variables
#'   corresponding to one or more of the items in `tax_level_labels`.
#' @param entity_id_col Character: a name, or vector of names, of variables in
#'   `data` that uniquely identifies entities. Default `NULL` assumes that each
#'   row is a unique entity.
#' @param tax_level_labels A character vector of all levels in the taxonomy.
#'   Default is \code{\link{chemont_tax_levels}} to use the levels of the
#'   ChemOnt taxonomy.
#' @return The input `data.frame`, with new variables `terminal_label`
#'   containing the terminal label for each entity, and `terminal_level`
#'   containing the terminal taxonomy level for each entity (as an integer).
#' @export
#' @importFrom magrittr `%>%`
#' @examples
#' #add terminal labels to the first ten chemicals in the `biosolids_class` data set
#'  add_terminal_label(data = biosolids_class[1:10, ])
add_terminal_label <- function(data,
                               entity_id_col = NULL,
                               tax_level_labels = chemont_tax_levels){

  data <- as.data.frame(data)

  #if no entity ID column specified,
  #create one with row numbers
  id_null <- FALSE
  if(is.null(entity_id_col)){
    id_null <- TRUE
    #add a new variable to the data
    #ensure it does not conflict with any of the existing variable names
    entity_id_col <- rev(
      make.names(
        names = c(names(data),
                  "id"
        ),
        unique = TRUE
      )
    )[1]
    data[[entity_id_col]] <- 1:nrow(data)
    }



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
  }

  if("terminal_level" %in% names(data)){
    warning(paste("Column 'terminal_level' already exists",
                  "in the input data.frame;",
                  "it will be overwritten"))
    data[["terminal_level"]] <- NULL
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

  dat_orig <- data #save original input data


  labels <- data %>%
    dplyr::select(dplyr::all_of(c(entity_id_col,
                                  tax_level_labels))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(tax_level_labels),
                        names_to = "tax_level",
                        values_to = "label") %>%
    dplyr::filter(!is.na(label)) %>% #remove unused levels
    dplyr::distinct( #group by entity
      dplyr::pick(
        dplyr::all_of(entity_id_col),
        label,
        tax_level
      )
    ) %>%
    dplyr::group_by(
      dplyr::pick(
        dplyr::all_of(entity_id_col)
      )
    ) %>%
    dplyr::slice_tail() %>%  #take most-specific label for each item
    #(i.e. last row)
   dplyr::rename(terminal_label = label, #rename cols to refer to "terminal"
                  terminal_tax_level = tax_level) %>%
    dplyr::mutate(terminal_level = match(terminal_tax_level,
                                         tax_level_labels)) %>%
    dplyr::mutate(terminal_tax_level = NULL)


  #merge terminal label & terminal level info back into original
  dat_out <- merge(dat_orig,
                    labels,
                    by = entity_id_col)

  #if a new ID variable was added, remove it
  if(id_null %in% TRUE){
    dat_out[[entity_id_col]] <- NULL
  }

  return(dat_out)
}


#'Calculate overlap
#'
#'Calculate overlap between two datasets
#'
#'Calculate overlap between two classified datasets at the individual entity
#'level
#'
#'@param data_1 A `data.frame` of classified entities (or something that can be
#'   coerced to one using [base::as.data.frame()])
#'@param data_2 A `data.frame` of classified entities (or something that can be
#'   coerced to one using [base::as.data.frame()])
#'@param entity_id_col Character: Name of variable in `data_1` and `data_2` that
#'  identifies unique entities. Must be the same for both \code{data_1} and
#'  \code{data_2}. If `NULL` (default), each row is assumed to be a unique
#'  entity.
#'@param at_level Taxonomy level at which to calculate overlap. Default
#'  \code{"terminal"} calculates overlap for terminal labels. Otherwise, may be
#'  one of \code{tax_level_labels} to calculate overlap at a more-general level
#'  of the taxonomy.
#'@param tax_level_labels Taxonomy levels. Default [chemont_tax_levels].
#'@return A `data.frame` with a number of rows equal to the number of unique
#'  labels at the specified level that occur either in \code{data_1} or
#'  \code{data_2}. The first variable is named with the value of
#'  \code{at_level}, and contains the unique labels at that level that occur
#'  either in \code{data_1} or \code{data_2}. The other variables are:
#'  \item{n_1}{The number of entities for this label in \code{data_1}}
#'  \item{n_2}{The number of entities for this label in \code{data_2}}
#'  \item{n_intersect}{The number of entities for this label that are in both
#'  \code{data_1} and \code{data_2}}
#'   \item{n_union}{The number of entities for
#'  this label that are in either \code{data_1} or \code{data_2}}
#'  \item{simil}{The Jaccard similarity of the sets of entities in \code{data_1}
#'  and \code{data_2} for each label, \code{n_intersect / n_union}}
#' @export
#' @author Caroline Ring, Paul Kruse
#' @examples
#' #show the overlap between BIOSOLIDS 2021 and USGS_WATER data sets
#' #for organic and inorganic chemicals
#' calc_number_overlap(data_1 = biosolids_class,
#'  data_2 = usgs_class,
#'  at_level = "kingdom",
#'  tax_level_labels = chemont_tax_levels)
#'
calc_number_overlap <- function(data_1,
                                data_2,
                                entity_id_col = NULL,
                                at_level = "terminal",
                                tax_level_labels = chemont_tax_levels){

  #if no entity ID column specified, then set it to be all shared variables between the two data sets
  if(is.null(entity_id_col)){
  entity_id_col <- setdiff(intersect(names(data_1),
                             names(data_2)),
                           tax_level_labels)
  }

  if(at_level %in% "terminal"){
  #get terminal labels if not already there
  if(!("terminal_label" %in% names(data_1))){
    data_1 <- add_terminal_label(data_1,
                                 entity_id_col = entity_id_col,
                                    tax_level_labels = tax_level_labels)
  }

  if(!("terminal_label" %in% names(data_2))){
    data_2 <- add_terminal_label(data_2,
                                 entity_id_col = entity_id_col,
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

  dat1_grp <- data_1 %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(group_col)
      )
    )

  dat2_grp <- data_2 %>%
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(group_col)
      )
    )

  if(is.null(entity_id_col)){
    count_1 <- dat1_grp %>%
      dplyr::count() %>%
      dplyr::rename(n_1 = n)
    count_2 <- dat2_grp %>%
      dplyr::count() %>%
      dplyr::rename(n_2 = n)

  }else{
    count_1 <- dat1_grp %>%
      dplyr::summarise(
        n_1 = dplyr::n_distinct(
          dplyr::across(
            dplyr::all_of(entity_id_col)
          )
        )
      )

    count_2 <- dat2_grp %>%
      dplyr::summarise(
        n_2 = dplyr::n_distinct(
          dplyr::across(
            dplyr::all_of(entity_id_col)
          )
        )
      )
  }

  df_count <- merge(count_1,
                    count_2,
                    by = group_col,
                    all = TRUE)


  #replace NAs with zeros for labels with no entities in one data set
  df_count$n_1[is.na(df_count$n_1)] <- 0
  df_count$n_2[is.na(df_count$n_2)] <- 0

  df_list <- sapply(df_count[[group_col]],
                    get_overlap,
                    group_col = group_col,
                    entity_id_col = entity_id_col,
                    data_1 = data_1,
                    data_2 = data_2,
                    simplify = FALSE,
                    USE.NAMES = TRUE)

  overlap_df <- dplyr::bind_rows(df_list)
  overlap_df <- setNames(overlap_df,
                         c(group_col,
                           "n_intersect",
                           "n_union",
                           "simil"))

  outdf <- merge(df_count,
                 overlap_df,
                 by = group_col,
                 all = TRUE)

}

#'Get overlap
#'
#'Helper function for [calculate_overlap()]
#'
#'Calculates number of overlapping entities in two data sets for a particular
#'classification label.
#'
#'@param grouplab Character: the group label
#'@param group_col Character: the variable name for the group label (i.e., the
#'  taxonomy level)
#'@param entity_id_col Character vector: variable name(s) in `data_1` and
#'  `data_2` that define the unique entity ID
#'@param data_1 A data.frame of classified entities
#'@param data_2 A data.frame of classified entities
#'@return A `data.frame` with variables `group` (containing `grouplab`),
#'  `n_intersect` (containing the number of intersecting entities in that group
#'  for the two datasets), `n_union` (containing the number of entities in the
#'  union of that group for the two datasets), and `simil` (containing the
#'  fractional similarity or overlap of entities in that group between the two
#'  datasets)
#' @examples
#' get_overlap(grouplab = "Organic compounds",
#' group_col = "kingdom",
#' data_1 = biosolids_class,
#' data_2 = usgs_class)
#'
get_overlap <-  function(grouplab,
                         group_col,
                         entity_id_col = NULL,
                         data_1,
                         data_2){
  #if no entity ID column specified, then set it to be all shared variables between the two data sets
  if(is.null(entity_id_col)){
    entity_id_col <- intersect(names(data_1), names(data_2))
  }

  id_1 <- data_1[data_1[[group_col]] %in% grouplab, ]
  id_2 <- data_2[data_2[[group_col]] %in% grouplab,]

  n_intersect <- dim(dplyr::inner_join(id_1,
                                       id_2,
                                       by = entity_id_col))[[1]]
  n_union <- dim(dplyr::full_join(id_1,
                                  id_2,
                                  by = entity_id_col))[[1]]
  simil <- n_intersect/n_union
  data.frame("group" = grouplab,
             "n_intersect" = n_intersect,
             "n_union" = n_union,
             "simil" = simil
  )
}


#' Get labels at a specified taxonomy level
#'
#' Get all the unique labels for a given taxonomy level in a given set of
#' classified entities.
#'
#' @param data A `data.frame` of classified entities (or something that can be
#'   coerced to one using [base::as.data.frame()])
#' @param level_label Character: A taxonomy level of the classified data.
#' @param tax_level_labels A character vector of all of the taxonomy levels, in order.
#'   Default is \code{\link{chemont_tax_levels}}, the levels of the ChemOnt
#'   taxonomy.
#' @return Character vector: The unique labels in `data` at the given taxonomy
#'   level.
#' @examples
#' #get all the superclass labels in the BIOSOLIDS2021 data
#' get_label_level(data = biosolids_class,
#' level_label = "superclass")
#'
get_label_level <- function(data,
                            level_label,
                            tax_level_labels = chemont_tax_levels){
  data <- as.data.frame(data)

  if (!(level_label %in% names(data) | !(level_label %in% tax_level_labels)))
    stop(paste('Please input a valid label!', level_label))

  # Collect the labels
  labels <- unique(data[[level_label]])

  # Remove NA from list of labels
  labels <- labels[!is.na(labels)]

  # Remove '' from list of labels
  labels <- labels[which(labels != '')]

  return(labels)
}

#' Get labels
#'
#' Retrieve all labels in a `data.frame` of classified entities, grouped by
#' taxonomy level.
#'
#' @param data A `data.frame` of classified entities (or something that can be
#'   coerced to one using [base::as.data.frame()])
#' @param tax_level_labels A character vector of all of the taxonomy levels, in
#'   order. Default is \code{\link{chemont_tax_levels}}, the levels of the
#'   ChemOnt taxonomy.
#' @return A `list` the same length as `tax_level_labels`: character vectors
#'   containing the unique classification labels in the data for each level of
#'   taxonomy.
#' @examples
#' #get labels for the first ten chemicals in the classified biosolids set
#' get_labels(data = biosolids_class[1:10, ])
#'
get_labels <- function(data,
                       tax_level_labels = chemont_tax_levels){
  labels <- sapply(tax_level_labels,
                   function(t) {
                     get_label_level(data = data,
                                     level_label = t,
                                     tax_level_labels = tax_level_labels)
                     }
                   )
  labels
}

#' Get terminal labels
#'
#' Retrieve terminal labels in a `data.frame` of classified entities.
#'
#' @param data A `data.frame` of classified entities. Variables must include at
#'   least one of the taxonomy level names in \code{tax_level_labels}.
#' @param entity_id_col Character: the name of the variable(s) in `data` with
#'   unique entity identifiers. If `NULL` (default), then each row is assumed to
#'   be a unique entity.
#' @param tax_level_labels A vector of taxonomy levels. Default is
#'   \code{\link{chemont_tax_levels}}, the levels of the ClassyFire taxonomy.
#' @return A vector of terminal classification labels, one for each entity in
#'   the input data.
#' @examples
#' #get terminal labels for the first ten chemicals in the classified biosolids set
#' get_terminal_labels(data = biosolids_class[1:10, ])
#'
get_terminal_labels <- function(data,
                                entity_id_col = NULL,
                       tax_level_labels = chemont_tax_levels){
  labels <- add_terminal_label(data = data,
                     entity_id_col = entity_id_col,
                     tax_level_labels = tax_level_labels)[["terminal_label"]]
  return(labels)
}

#' Count entities per label
#'
#' Given a `data.frame` of classified entities, determine the number of entities
#' with each label at each level of the taxonomy.
#'
#' @param data A `data.frame` of classified entities.
#' @param entity_id_col `NULL` (default) or character vector giving the name(s)
#'   of variables in `data` that specify unique entities. If `NULL`, each row
#'   will be treated as a unique entity.
#' @param tax_level_labels A character vector of taxonomy levels. Default is
#'   \code{\link{chemont_tax_levels}}, the levels of the ClassyFire taxonomy.
#' @return A list the same length as , and named after, `tax_level_labels`,
#'   where each element is a `data.frame` with first variable `label`, giving
#'   the labels, and second variable `n`, giving the count of entities with each
#'   label in the data.
#' @examples
#' #count entities per label in the first ten BIOSOLIDS2021 chemicals
#' count_entities_per_label(data = biosolids_class[1:10, ])
#'
count_entities_per_label <- function(data,
                              entity_id_col = NULL,
                                 tax_level_labels = chemont_tax_levels){

  #if no entity ID column specified,
  #create one with row numbers
  id_null <- FALSE
  if(is.null(entity_id_col)){
    id_null <- TRUE
    #add a new variable to the data
    #ensure it does not conflict with any of the existing variable names
    entity_id_col <- rev(
      make.names(
        names = c(names(data),
                  "id"
        ),
        unique = TRUE
      )
    )[1]
    data[[entity_id_col]] <- 1:nrow(data)
  }


label_counts <-  sapply(tax_level_labels,
         function(this_level){
           if(this_level %in% names(data)){
           data %>%
             dplyr::group_by(
               dplyr::pick(
                 dplyr::all_of(this_level)
               )
             ) %>%
            dplyr::summarise(
              n = dplyr::n_distinct(
                dplyr::pick(
                  dplyr::all_of(entity_id_col)
                )
              )
            ) %>%
             dplyr::ungroup() %>%
               as.data.frame() %>%
             setNames(c("label",
                        "n"))
           }else{
             #if no variable for this level in the data,
             #return data.frame with label NA and n = number of rows
             #this is the same behavior as if a variable were present but all-NA
             data.frame(label = NA_character_,
                        n = nrow(data))
           }
         },
        simplify = FALSE,
        USE.NAMES = TRUE)

  return(label_counts)
}

#' Count unique labels
#'
#' Count the number of unique labels per taxonomy level in a classified data set.
#'
#' @param data A `data.frame` of classified entities.
#' @param tax_level_labels A character vector of taxonomy levels. Default is
#'   \code{\link{chemont_tax_levels}}, the levels of the ClassyFire taxonomy.
#' @return The number of labels per taxonomy level.
#' @examples
#' #count unique labels in the first ten BIOSOLIDS2021 chemicals
#' count_labels(data = biosolids_class[1:10, ])
#'
count_labels <- function(data,
                             tax_level_labels = chemont_tax_levels){

  lengths <- sapply(get_labels(data = data,
                               tax_level_labels = tax_level_labels),
                    length)
  lengths
}

#' Similarity for two datasets
#'
#' Calculate similarity measures for two datasets
#'
#' @param data_1 A `data.frame` of classified entities
#' @param data_2 Another `data.frame` of classified entities
#' @param terminal_label The variable name in the two data.frames that denotes
#'   the terminal label of the classification. Default `"terminal_label"`.
#' @param tree The taxonomy tree to use as a \code{\link[ape]{phylo}}-class object (see
#'   [ape::read.tree()] for description of this class). Default \code{\link{chemont_tree}}.
#' @param tax_level_labels The set of taxonomy levels to use. Default
#'   \code{\link{chemont_tax_levels}}.
#' @param sim_matrix Optional: A pre-computed similarity matrix to use as a
#'   lookup table. Default `NULL`, which will compute similarity from scratch
#'   using the metric specified in `similarity`. If non-`NULL`, will override
#'   anything specified in `similarity` (with a warning). For example, if `tree
#'   = chemont_tree`, you could use `sim_matrix = chemont_jaccard` to lookup
#'   values in the pre-computed pairwise Jaccard similarity matrix for all nodes
#'   in the ChemOnt taxonomy tree. Providing a pre-computed similarity matrix
#'   may be faster if `data_1` and `data_2` are large.
#' @param similarity The similarity metric to calculate. Used only if
#'   `sim_matrix = NULL`. Default `NULL`, which will calculate Jaccard
#'   similarity. Options include "jaccard", "resnik", "lin", and
#'   "jiang_conrath". Note: Ignored if `sim_matrix` is non-`NULL`, with a
#'   warning.
#' @return A similarity matrix with rows and columns corresponding to labels
#'   from the `terminal_label` variable in `data_1` and `data_2`, respectively.
#' @author Caroline Ring, Paul Kruse
#' @export
#'
#' @examples
#' #computing similarity from scratch
#' calc_similarity_data(data_1 = biosolids_class[1:10,], data_2 = usgs_class[1:10, ], sim_matrix = NULL, similarity = "jaccard")
#'
#' #using pre-computed similarity matrix as lookup table
#' calc_similarity_data(data_1 = biosolids_class[1:10,], data_2 = usgs_class[1:10, ], sim_matrix = chemont_jaccard)
#'
#' #providing both sim_matrix and similarity throws a warning
#' calc_similarity_data(data_1 = biosolids_class[1:10, ], data_2 = usgs_class[1:20, ], sim_matrix = chemont_jaccard, similarity = "resnik")
#'
#'
#' @seealso \code{\link{jaccard_similarity}}, \code{\link{resnik_similarity}},
#'   \code{\link{lin_similarity}}, \code{\link{jiang_conrath_similarity}},
#'   \code{\link{similarity_matrix}}
calc_similarity_data <- function(data_1,
                                 data_2,
                                 terminal_label = "terminal_label",
                                 tree = chemont_tree,
                                 tax_level_labels = chemont_tax_levels,
                                 sim_matrix = NULL,
                                 similarity = NULL){

  #calculate pairwise similarity of ancestry of terminal labels in two data sets

  #check for terminal_label
  if(terminal_label == "terminal_label"){
    if(!(terminal_label %in% names(data_1))){
      data_1 <- add_terminal_label(data = data_1,
                                   tax_level_labels = tax_level_labels)
    }

    if(!(terminal_label %in% names(data_2))){
      data_2 <- add_terminal_label(data = data_2,
                                   tax_level_labels = tax_level_labels)
    }
  }

  #Keep only data with terminal labels in the tree
  data_1 <- data_1[data_1[[terminal_label]] %in%
                     c(tree$tip.label, tree$node.label), ]
  data_2 <- data_2[data_2[[terminal_label]] %in%
                     c(tree$tip.label, tree$node.label), ]

  if(is.null(sim_matrix)){
    if(is.null(similarity)){
      warning("Defaulting to Jaccard similarity")
      similarity = "jaccard"
    }
    m <- similarity_matrix(nodes1 = data_1[[terminal_label]],
                           nodes2 = data_2[[terminal_label]],
                           tree = tree,
                           metric = similarity,
                           upper_tri = FALSE)

    rownames(m) <- data_1[[terminal_label]]
    colnames(m) <- data_2[[terminal_label]]
  }else{ #if sim.matrix provided, just use it as lookup table
    if(!is.null(similarity)){
      warning(
        paste0(
          "Both `sim_matrix` and `similarity` were provided.",
          " Ignoring `similarity = ",
          similarity, "` and using provided `sim_matrix` as lookup table.")
      )
    }
    m <- sim_matrix[data_1[[terminal_label]],
                    data_2[[terminal_label]]]
  }

  return(m)
}
