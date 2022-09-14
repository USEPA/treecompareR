#' Label bars
#'
#' This function takes in a (list of) data.table(s) and returns figures
#' illustrating the label numbers by taxonomy level.
#'
#' @param data A data.table or list of data.tables with chemicals and their
#'   classification data.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @return A list of two plots, the first showing number of labels by taxonomy
#'   level for each data.table and the second showing the number of labels by
#'   data.table for each taxonomy level.
#' @export
#' @import data.table
#' @importFrom tidyr pivot_longer
#' @import ggplot2
label_bars <- function(data = NULL,
                       tax_level_labels = chemont_tax_levels){
  tax_levels <- NULL
  count_sums <- NULL
  dataset <- NULL

  if (is.null(data))
    stop('Please input data!')
  if (data.table::is.data.table(data)){
    data <- list(data)
  } else {
    if (!is.list(data) | !all(sapply(data, data.table::is.data.table)))
      stop('Please input a single data.table or list of data.tables!')
  }

  if (is.null(tax_level_labels)){
    tax_level_labels <- c('kingdom', 'superclass', 'class', 'subclass',
                          'level5', 'level6', 'level7', 'level8',
                          'level9', 'level10', 'level11')
  }

  number <- length(data)

  if (number == 1){
    data_names <- c('Set_1')
  } else {
    if (is.null(names(data)) | length(names(data)) != number | any(is.na(names(data)))){
      data_names <- paste0('Set_', seq_len(number))
    } else {
      data_names <- names(data)
    }
  }

  df <- data.frame(tax_level_labels, unname(sapply(data, function(t) get_label_length(get_terminal_labels(t, tax_level_labels)))))
  names(df) <- c('tax_levels', data_names)
  transformed_df <- df %>%
    tidyr::pivot_longer(!tax_levels, names_to = 'dataset', values_to = 'count_sums')
  transformed_df["count_sums"] <- as.numeric(transformed_df$count_sums)
  transformed_df["tax_levels"] <- factor(transformed_df$tax_levels, levels = tax_level_labels)

  plot_1 <- ggplot(transformed_df) +
    facet_wrap(~dataset, scales = 'free') +
    geom_bar(aes(x = tax_levels, y = count_sums, fill = tax_levels), stat = 'identity') +
    scale_color_manual(name = 'Taxonomy level') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = 'Taxonomy levels',
         y = 'Number of unique labels',
         fill = 'Taxonomy levels')

  plot_2 <- ggplot(transformed_df) +
    facet_wrap(~tax_levels) +
    geom_bar(aes(x = dataset, y = count_sums, fill = dataset), stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = 'Data set',
         y = 'Number of unique labels',
         fill = 'Data set')

  return(list(plot_1, plot_2))
}


#'Display subtree
#'
#'This function takes in a base tree, and (optionally) one or two data sets
#'consisting of entities with classification labels. It returns a tree diagram
#'consisting of the base tree; if data sets were provided, branches of the base
#'tree will be highlighted according to their membership in the provided data
#'sets.
#'
#'# How to specify data sets
#'
#'## As `data.frame` objects containing classified entities
#'
#'If `data_1` or `data_2` are provided as `data.frame` objects containing
#'classified entities, they must be of the following format: Rows correspond to
#'individual entities. There must be one column corresponding to each of the
#'taxonomy levels specified in `tax_level_labels`, containing the label at that
#'level for each entity. These columns must have names corresponding to the
#'members of `tax_level_labels`. For example, to use the ChemOnt taxonomy, the
#'`data.frame` must have columns `kingdom`, `superclass`, `class`, `subclass`,
#'`level5`, `level6`, ... `level11`. Whenever an entity does not have a label at
#'a given level, the corresponding element in the data table for that row and
#'column should be `NA_character_`.
#'
#'The `data.frame` must have at least one additional column that uniquely
#'identifies individual entities; the name of that additional column does not
#'matter, as long as it is not the same as one of the taxonomy levels.
#'
#'For an example of a properly-formatted input `data.frame`, see the built-in
#'dataset \code{\link{biosolids_class}}. The only required columns in
#'\code{\link{biosolids_class}} are \code{kingdom, superclass, class, subclass,
#'level5, level6, level7, level8, level9, level10, level11} and any one of the
#'chemical identifier columns, for example \code{DTXSID}.
#'
#'## As vectors of node labels
#'
#'If `data_1` or `data_2` are provided as vectors of node labels, they refer to
#'node labels in \code{base_tree}, i.e. labels that appear in
#'\code{base_tree$tip.label} or \code{base_tree$node.label}. Any labels in
#'`data_1` or `data_2` that do not appear in the base tree will be ignored.
#'
#'The labels in `data_1` and/or `data_2`, plus all of their ancestors in the
#'base tree, will be highlighted. Their descendants will not be highlighted
#'unless those descendants themselves appear in `data_1` or `data_2`.
#'
#'## As vectors of node numbers
#'
#'If `data_1` or `data_2` are provided as vectors of node numbers, they refer to
#'node and tip index numbers in the base tree. Tip numbers run from 1 to
#'\code{ape::Ntip(base_tree)}, and node numbers run from
#'\code{ape::Ntip(base_tree) + 1} to \code{ape::Ntip(base_tree) +
#'ape::Nnode(base_tree) - 1}.
#'
#'You might specify `data_1` or `data_2` as vectors of node numbers if, for
#'example, you have used \code{\link[phangorn]{Ancestors}} or
#'\code{\link[phangorn]{Descendants}} to trace ancestors or descendants of a
#'specified node, and now you want to visualize those ancestors or descendants.
#'\code{\link[phangorn]{Ancestors}} and \code{\link[phangorn]{Descendants}} both
#'return lists or vectors of node numbers, so it would be convenient to be able
#'to pass those vectors of node numbers directly to \code{display_subtree}.
#'
#''The nodes/tips in `data_1` and/or `data_2`, plus all of their ancestors in
#'the base tree, will be highlighted. Their descendants will not be highlighted
#'unless those descendants themselves appear in `data_1` or `data_2`.
#'
#'## As tree objects
#'
#'If `data_1` or `data_2` are provided as \code{\link[ape]{phylo}}-class tree
#'objects, any labels in their \code{tip.label} or \code{node.label} elements
#'that do not appear in the tip or node labels of the base tree will be ignored.
#'In other words, only the portions of `data_1` and `data_2` that are subtrees
#'of `base_tree` wil be used.
#'
#'@param base_tree The "base tree" to plot, as a \code{\link[ape]{phylo}}-class
#'  object.  Default is the full ChemOnt taxonomy tree,
#'  \code{\link{chemont_tree}}.
#'@param base_name Will be used as the plot title. Usually this should name or
#'  describe the base tree. Default NULL will result in no plot title.
#'@param data_1 Optional: Highlight branches of the base tree according to their
#'  membership in this list. Default is \code{NULL}, to do no highlighting. If
#'  not \code{NULL}, must be one of the following options: A `data.frame`
#'  consisting of a list of entities with classification data; a vector of node
#'  numbers in the base tree; a vector of node labels in the base tree; or a
#'  subtree of `base_tree` as a \code{\link[ape]{phylo}}-class object. See Details.
#'@param data_2 Optional: Highlight branches of the base tree to compare
#'  membership in this list and \code{data_1}. Default is \code{NULL}, to do no
#'  highlighting. If not \code{NULL}, one of the following options: A
#'  `data.frame` consisting of a list of chemicals with classification data; a
#'  list of node numbers in the base tree; a list of node labels in the base
#'  tree; or a subtree of `base_tree` as a \code{\link[ape]{phylo}}-class
#'  object. See Details.
#'@param name_1 Optional: A string giving the name of the list in \code{data_1},
#'  for plot labeling. Default is "Set 1".
#'@param name_2 Optional: A string giving the name of the list in \code{data_2}
#'  (if any), for plot labeling. Default is "Set 2".
#'@param tax_level_labels Optional: a vector of the taxonomy level labels to be
#'  used. Default value: \code{\link{chemont_tax_levels}}, i.e., the levels of
#'  the ClassyFire taxonomy: \code{c("kingdom", "superclass", "class",
#'  "subclass", paste0("level", 5:11))}.
#'@param layout \code{\link{ggtree}} layout option. Default "circular."
#'@param base_color Color for base tree (branches not in any input data set).
#'  Default "gray80" for a light gray.
#'@param subtree_colors Vector of colors to use to highlight branches according
#'  to subtree membership. If \code{data_2} is NULL, only the first element of
#'  \code{subtree_colors} will be used (to signify that a branch is in
#'  \code{data_1}). If \code{data_2} is not NULL, then the colors will be used
#'  in the following order: branch in \code{data_1} only; branch in
#'  \code{data_2} only; branch in both \code{data_1} and \code{data_2}. Default
#'  is \code{c("#66C2A5", "#FC8D62","#8DA0CB")}, the first 3 colors of the
#'  ColorBrewer2 "Set2" set.
#'@param base_size Line width for base tree. If NULL, default line width will be
#'  used.
#'@param subtree_sizes Vector of line widths to use to highlight branches
#'  according to subtree membership. If \code{data_2} is NULL, only the first
#'  element of \code{subtree_sizes} will be used (to signify that a branch is in
#'  \code{data_1}). If \code{data_2} is not NULL, then the sizes will be used in
#'  the following order: branch in \code{data_1} only; branch in \code{data_2}
#'  only; branch in both \code{data_1} and \code{data_2}. If NULL, then size
#'  will not be used to highlight branches; all branches will be plotted with
#'  line width given by \code{base_size}. Default is NULL.
#'@param clade_level The taxonomy level at which to draw clade labels, if
#'  any. Root is level 0. Default is level 2 (superclass, in ChemOnt taxonomy).
#'  Set to NULL to suppress clade labels altogether.
#'@param clade_label_fontsize Font size for clade labels. Default 3.
#'@param clade_label_wrap Number of characters at which to wrap clade labels to
#'  a second line. Default 20.
#'@param clade_label_lineheight Line height for multi-line clade labels, as a
#'  multiple of the size of text. Controls vertical space between lines on
#'  multi-line clade labels. Default 0.7.
#'@return A \code{\link[ggtree]{ggtree}} object visualizing the full base tree,
#'  with branches highlighted to indicate presence in
#'  \code{data_1}, in \code{data_2} if supplied, neither, or both.
#'@export
#'@import ggtree
display_subtree <- function(base_tree = chemont_tree,
                            base_name = NULL,
                            data_1 = NULL,
                            data_2 = NULL,
                            name_1 = "Set 1",
                            name_2 = "Set 2",
                            tax_level_labels = chemont_tax_levels,
                            layout = "circular",
                            base_opts = list("color" = "black",
                                             "size" = 0.5,
                                             "linetype" = 1),
                            subtree_mapping = list(color = "default"),
                            clade_level = "auto",
                            clade_opts = list(wrap = 20,
                                              barsize = "alternate",
                                              fontsize = 3,
                                              lineheight = 0.7,
                                              default_to_tip = TRUE)
){

  #Keep defaults for any base_opts not otherwise specified
  base_opts_default <- list("color" = "black",
                                        "size" = 0.5,
                                        "linetype" = 1)
  base_opts <- c(base_opts,
                 base_opts_default[setdiff(names(base_opts_default),
                                           names(base_opts))])

  #Keep defaults for any subtree_mapping not otherwise specified
  subtree_mapping_default <- list(color = "default")
  subtree_mapping <- c(subtree_mapping,
                       subtree_mapping_default[
                         setdiff(names(subtree_mapping_default),
                                 names(subtree_mapping)
                         )
                       ]
  )


  if(is.null(data_1)){
    #ignore any default subtree mapping
    subtree_mapping <- NULL
  }else{
    if(is.null(subtree_mapping)){
      warning(paste("Subtree data was supplied, but subtree_mapping is NULL.",
                    "No subtrees will be highlighted."))
    }
    #Check subtree mapping names
    #They need to be valid aesthetics for ggtree
    bad_subtree_map <- any(c("color", "colour", "size", "linetype") %in%
                             names(subtree_mapping))

    if(!is.list(subtree_mapping) |
       is.null(names(subtree_mapping)) |
       !isTRUE(bad_subtree_map)){
      stop(paste("subtree_mapping should be a list with named elements.",
      "Names must be one or more of 'color' (or 'colour'), 'size', and/or 'linetype'. "))
    }else{
    #if something reasonable was provided for subtree_map,
      #check if any other aesthetics were provided and will be ignored
    bad_aes <- setdiff(names(subtree_mapping),
                       c("color", "colour", "size", "linetype"))
    good_aes <- intersect(names(subtree_mapping),
                       c("color", "colour", "size", "linetype"))
    if(length(bad_aes)>0){
      message(paste("In subtree_mapping, only aesthetics",
                    paste(good_aes, collapse = ", "),
                    "will be used. Aesthetics",
                    paste(bad_aes, collapse = ", "),
                    "were provided but will be ignored,",
                    "since ggtree::geom_tree() does not understand them."))
    }

    if("color" %in% names(subtree_mapping)){
      if (all(subtree_mapping$color %in% "default")){
        subtree_mapping$color <- c("gray70",
                                   "#66C2A5",
                                   "#8DA0CB",
                                   "#FC8D62")
      }
    }
  else if ("colour" %in% names(subtree_mapping)){
    if (all(subtree_mapping$colour %in% "default")){
      subtree_mapping$colour <- c("gray70",
                                 "#66C2A5",
                                 "#8DA0CB",
                                 "#FC8D62")
    }
  }

    #check lengths of subtree_mapping vs. number of datasets provided
    sm_length <- sapply(subtree_mapping,
                        length,
                        USE.NAMES = TRUE)

    if(!is.null(data_1) & !is.null(data_2)){
      if(any(sm_length < 4)){
        short_el <- names(subtree_mapping)[sm_length < 4]
        message(paste("Both data_1 and data_2 were provided, but subtree_mapping elements",
                paste(short_el, collapse = "; "),
                "have fewer than 4 elements. They will be recycled to length 4."))
        subtree_mapping[short_el] <- sapply(subtree_mapping[short_el],
                                            function(x) rep(x, length.out = 4),
                                            simplify = FALSE,
                                            USE.NAMES =  TRUE)
      }
    }else if(!is.null(data_1) & is.null(data_2))
      if(any(sm_length > 2)){
        long_el <- names(subtree_mapping)[sm_length > 2]
        message(paste("data_1 was provided, but subtree_mapping elements",
                paste(long_el, collapse = "; "),
                "have more than 2 elements. Only the first 2 elements will be used."))
        subtree_mapping[long_el] <- sapply(subtree_mapping[long_el],
                                            function(x) x[1:2],
                                            simplify = FALSE,
                                            USE.NAMES =  TRUE)
      }

    if(any(sm_length < 2)){
      short_el <- names(subtree_mapping)[sm_length < 2]
      message(paste("data_1 was provided, but subtree_mapping elements",
              paste(short_el, collapse = "; "),
              "have fewer than 2 elements. They will be recycled to length 2."))
      subtree_mapping[short_el] <- sapply(subtree_mapping[short_el],
                                          function(x) rep(x, length.out = 2),
                                          simplify = FALSE,
                                          USE.NAMES =  TRUE)
    }
    } #end if(!is.list(subtree_mapping) |
    # is.null(names(subtree_mapping)) |
    #   !isTRUE(bad_subtree_map))



  } #end if(is.null(data_1))



#Base tree
  #Use base options, unless they will be mapped to list presence later
  #(aes() doesn't seem to overwrite them as expected)
  #e.g. if subtree_mapping has a "color" element, don't use base_opts$color
  tree_plot <- do.call(ggtree,
                       c(list(tr = base_tree,
                      layout = layout),
                      base_opts[setdiff(names(base_opts),
                                        names(subtree_mapping))
                                ]
                      )
                      )

  if(!is.null(data_1)){

    #If it is a list, look at the data type of the list elements,
    #and concatenate them as appropriate
    if(class(data_1) %in% "list"){
      #use rapply to handle possible nested list
      el_class <- unique(rapply(data_1, class))
      if(length(unique(el_class))>1){
        stop("data_1 is a list, but not all list elements are of the same class.")
      }
      if(all(el_class %in% "data.frame")){
        #rbind the data.frames
        data_1 <- as.data.frame(dplyr::bind_rows(data_1))
      }else if(all(el_class %in% c("character", "numeric"))){
        #concatenate, recursively if necessary
        data_1 <- unlist(data_1, recursive = TRUE)
      }else if(all(el_class %in% "phylo")){
        #pull labels from each of the trees and combine
        data_1 <- unlist(lapply(data_1,
                         function(x){
                           c(x$tip.labels,
                             x$node.labels)
                         }))
      }else{
        stop("data_1 is a list, but one or more elements are not one of the recognized classes: data.frame, numeric, character, or phylo. ")
      }
      #and keep only the unique combined elements
      data_1 <- unique(data_1)
    } #end if(class(data_1) %in% "list")

    if(is.data.frame(data_1)){
      #Check formatting
      #Check that all taxonomy levels have a column
      if(!(all(tax_level_labels %in% names(data_1)))){
        stop(paste("When data_1 is a data.frame",
                   "it must have one column named for each taxonomy level",
                   "as defined in tax_level_labels, here:\n",
                   paste0(paste(tax_level_labels, collapse = ", "), "\n"),
                   "data_1 does not have columns for the following levels:\n",
                   paste(setdiff(tax_level_labels,
                                 names(data_1)),
                         collapse = ", ")))
      }else{ #if all taxonomy levels have a column,
        #Check for blank labels
        #If any, replace with NAs
        data_1 <- data_1 %>%
          dplyr::mutate(dplyr::across(
            dplyr::all_of(tax_level_labels),
            function(x) {
              x[!nzchar(trimws(x))] <- NA_character_
              return(x)
            }
          )
          )

        #If there is no entity ID column,
        #throw a warning and add one
        if(length(setdiff(names(data_1),
                          tax_level_labels))<1){
          warning(paste("data_1 is a data.frame,",
          "but it has no column to identify entities.",
          "Each row will be assumed to be one entity."))
          data_1$id <- 1:nrow(data_1)
        }

      } #end if(!(all(tax_level_labels %in% names(data_1))))

  #Get node numbers of data_1 subtree
  data_1_nodes <- get_subtree_nodes(data = data_1,
                                  base_tree = base_tree,
                                  tax_level_labels = tax_level_labels)
    }else if(is.numeric(data_1)){
      #interpret as node numbers

    data_1_nodes <- data_1

    }else if(is.character(data_1)){
    #interpret as node labels
      data_1_nodes <- get_node_from_label(label = data_1,
                                          tree = base_tree)

    }else if("phylo" %in% class(data_1)){
      #take node and tip labels
      data_1_labels <- c(data_1$tip.label,
                         data_1$node.label)
      #get nodes
      data_1_nodes <- get_node_from_label(label = data_1_labels,
                                          tree = base_tree)
    }else{
      stop("data_1 is not one of the recognized classes: data.frame, numeric, character, or phylo.")
    }

    #get ancestors of data_1_nodes
    data_1_all <- unique(
      c(data_1_nodes,
        unlist(
          phangorn::Ancestors(x = base_tree,
                              node = data_1_nodes)
        )
      )
    )

  #Data frame with all node numbers in taxonomy tree
   cohort_data <- get_tree_df(tree = base_tree)
   #Categorical column: is each node in Data Set 1?
   #0 = no, 1 = yes
   cohort_data$inSet1 <- ifelse(cohort_data$node %in% data_1_all,
                                1L,
                                0L)


  if (!is.null(data_2)) {

    #If it is a list, look at the data type of the list elements,
    #and concatenate them as appropriate
    if(class(data_2) %in% "list"){
      #use rapply to handle possible nested list
      el_class <- unique(rapply(data_2, class))
      if(length(unique(el_class))>1){
        stop("data_2 is a list, but not all list elements are of the same class.")
      }
      if(all(el_class %in% "data.frame")){
        #rbind the data.frames
        data_2 <- as.data.frame(dplyr::bind_rows(data_2))
      }else if(all(el_class %in% c("character", "numeric"))){
        #concatenate, recursively if necessary
        data_2 <- unlist(data_2, recursive = TRUE)
      }else if(all(el_class %in% "phylo")){
        #pull labels from each of the trees and combine
        data_2 <- unlist(lapply(data_2,
                                function(x){
                                  c(x$tip.labels,
                                    x$node.labels)
                                }))
      }else{
        stop("data_2 is a list, but one or more elements are not one of the recognized classes: data.frame, numeric, character, or phylo. ")
      }
      #and keep only the unique combined elements
      data_2 <- unique(data_2)
    } #end if(class(data_2) %in% "list")

    if(is.data.frame(data_2)){
      #Check formatting
      #Check that all taxonomy levels have a column
      if(!(all(tax_level_labels %in% names(data_2)))){
        stop(paste("When data_2 is a data.frame",
                   "it must have one column named for each taxonomy level",
                   "as defined in tax_level_labels, here:\n",
                   paste0(paste(tax_level_labels, collapse = ", "), "\n"),
                   "data_2 does not have columns for the following levels:\n",
                   paste(setdiff(tax_level_labels,
                                 names(data_2)),
                         collapse = ", ")))
      }else{ #if all taxonomy levels have a column,
        #Check for blank labels and replace with NAs
        data_2 <- data_2 %>%
          dplyr::mutate(dplyr::across(
            dplyr::all_of(tax_level_labels),
            function(x) {
              x[!nzchar(trimws(x))] <- NA_character_
              return(x)
            }
          )
          )

        #If there is no entity ID column,
        #throw a warning and add one
        if(length(setdiff(names(data_2),
                          tax_level_labels))<1){
          warning(paste("data_2 is a data.frame,",
                        "but it has no column to identify entities.",
                        "Each row will be assumed to be one entity."))
          data_2$id <- 1:nrow(data_2)
        }

      } #end if(!(all(tax_level_labels %in% names(data_2))))
      #Get node numbers of data_2 subtree
      data_2_nodes <- get_subtree_nodes(data = data_2,
                                      base_tree = base_tree,
                                      tax_level_labels = tax_level_labels)
    }else if(is.numeric(data_2)){
      #interpret as node numbers
      #get ancestors of these nodes
      data_2_nodes <- data_2
    }else if(is.character(data_2)){
      #interpret as node labels
      data_2_nodes <- get_node_from_label(label = data_2,
                                          tree = base_tree)
    }else if("phylo" %in% class(data_2)){
      #take node and tip labels
      data_2_labels <- c(data_2$tip.label,
                         data_2$node.label)
      #get nodes
      data_2_nodes <- get_node_from_label(label = data_2_labels,
                                          tree = base_tree)
    }else{
      stop("data_2 is not one of the recognized classes: data.frame, numeric, character, or phylo.")
    }

    #get ancestors of these nodes
    data_2_all <- unique(c(data_2_nodes,
                           unlist(phangorn::Ancestors(x = base_tree,
                                                      node = data_2_nodes))))

    #Categorical column: is each node in Data Set 2?
    #0 = no, 2 = yes
    cohort_data$inSet2 <- ifelse(cohort_data$node %in% data_2_all,
                                 2L,
                                 0L)

    #Categorical column: Is each node in Set 1, Set 2, neither, or both?
    #0 = neither
    #1 = Set 1 only
    #2 = Set 2 only
    #3 = both
    cohort_data$list_presence <- factor(cohort_data$inSet1 + cohort_data$inSet2,
                                 levels = 0:3,
                                 labels = c("Neither set",
                                            paste(name_1, "only"),
                                            paste(name_2, "only"),
                                            paste(name_1, "and", name_2)))

  }else{
    #Categorical column: Is each node in Data Set 1 or not?
    cohort_data$list_presence <- factor(cohort_data$inSet1,
                               levels = 0:1,
                               labels = c(paste("Not in", name_1),
                                          paste("In", name_1)))
  }

   #Name the subplot_mapping items after the categories in cohort_data$list_presence
   subtree_mapping <- sapply(subtree_mapping,
                             function(x) setNames(x, levels(cohort_data$list_presence)),
                             simplify = FALSE,
                             USE.NAMES = TRUE)
   #set up for aes call -- list of aesthetic mappings,
   #all applied to "list_presence" in cohort_data
   subtree_aes <- replicate(n= length(subtree_mapping),
                             expr = quote(list_presence))
   #name the list elements after the aesthetics in subtree_mapping
   subtree_aes <- setNames(subtree_aes, names(subtree_mapping))
   #you end up with something like:
   #`subtree_aes <- list(color = quote(list_presence),
   #                   size = quote(list_presence))``
   #`do.call(aes, subtree_aes)` is then equivalent to:
   #`aes(color = list_presence, size = list_presence)`

   #Prepare a list of manual scales as provided in subtree_mapping
   scale_list <- lapply(names(subtree_mapping),
                        function(aesthetic) {
                          ggplot2::scale_discrete_manual(
                            aesthetics = aesthetic,
                            name = "List presence",
                            values =  subtree_mapping[[aesthetic]],
                            breaks = levels(cohort_data$list_presence),
                            limits = levels(cohort_data$list_presence)
                          )
                        }
   )
   #The result of the above is something like
#    scale_list <- list(
#      color = ggplot2::scale_discrete_manual(
#        aesthetics = "color",
#    name = "List presence",
#    values = c("gray70",
#               "#66C2A5",
#               "#8DA0CB",
#               "#FC8D62"),
#    breaks = c("Neither set",
#               "In Set1",
#               "In Set2",
#               "Both sets"),
#    limits = c("Neither set",
#               "In Set1",
#               "In Set2",
#               "Both sets")
#    ),
#    size = ggplot2::scale_discrete_manual(
#      aesthetics = "size",
#      name = "List presence",
#      values = c(0.5,
#                 1,
#                 1,
#                 1),
#      breaks = c("Neither set",
#                 "In Set1",
#                 "In Set2",
#                 "Both sets"),
#      limits = c("Neither set",
#                 "In Set1",
#                 "In Set2",
#                 "Both sets")
#    )
# )
#But the idea is to generate it programatically from argument `subtree_mapping`

#add list presence highlighting to tree plot
   tree_plot <- tree_plot %<+% cohort_data +
     do.call(aes, #this maps list presence to all aesthetics specified in `subtree_mapping`
             subtree_aes) +
      scale_list + #this applies all of the scales in scale_list to the aesthetics
     theme(legend.position = "top")
  } #end if(!is.null(data_1))

   #if clade labels have been selected, add them to the plot
  if(!is.null(clade_level)){
  tree_plot <- add_cladelab(tree_plot = tree_plot,
                                          tree = base_tree,
                                          clade_level = clade_level,
                                          clade_opts = clade_opts)
  }

  #add title with base tree name, if provided
  if(!is.null(base_name)){
  tree_plot <- tree_plot + ggtitle(base_name)
  }

  return(tree_plot)

}


#'Prune and display subtree
#'
#'This function takes in a base tree, a pruning specification as for
#'\code{\link{prune_tree}}, and other arguments as for
#'\code{\link{display_subtree}}. It first prunes the base tree using
#'\code{\link{prune_tree}}, then calls \code{\link{display_subtree}}
#'substituting the pruned tree as the base tree for plotting.  It returns a tree
#'diagram consisting of the pruned base tree; if data sets were provided,
#'branches of the base tree will be highlighted according to their membership in
#'the provided data sets.
#'
#'
#'@param base_tree The "base tree" to prune and then plot, as a
#'  \code{\link[ape]{phylo}}-class object.  Default is the full ChemOnt taxonomy
#'  tree, \code{\link{chemont_tree}}.
#'@param prune_to As for \code{\link{prune_tree}}: What to *keep* from the base
#'  tree (everything else will be pruned away). May be a \code{data.frame} of
#'  classified data; one or more labels in the tree (tip or internal node
#'  labels); one or more node numbers in the tree (tip or internal nodes); or
#'  the name of a taxonomy level (one of the items in \code{tax_level_labels}).
#'  Default is NULL, which results in no pruning being done (i.e., the base tree
#'  is returned as-is). See Details.
#'@param prune_name Name of pruned tree to use as plot title. Default NULL
#'  results in no plot title.
#'@param adjust_branch_length Whether to adjust branch length so that all
#'  newly-pruned terminal nodes appear at the same length as tips, even if they
#'  were originally internal nodes. Default FALSE.
#'@param tax_level_labels Optional: a vector of the taxonomy level labels to be
#'  used. Default value: \code{\link{chemont_tax_levels}}, i.e., the levels of
#'  the ClassyFire taxonomy: \code{c("kingdom", "superclass", "class",
#'  "subclass", paste0("level", 5:11))}.
#'@param ... Other arguments to \code{\link{display_subtree}}.
#'@inheritDotParams display_subtree -base_tree -tax_level_labels
#'@return A ggtree object visualizing the pruned base tree. If\code{data_1}
#'  and/or \code{data_2} are supplied, branches will be highlighted to indicate
#'  whether they are present in each set, neither, or both.
#'@export
prune_and_display_subtree <- function(base_tree = chemont_tree,
                                      prune_to = NULL,
                                      prune_name = NULL,
                                      adjust_branch_length = FALSE,
                                      tax_level_labels = chemont_tax_levels,
                                      ...) { #args as for display_subtree()

  args <- list(...)

  pruned_tree <- prune_tree(tree = base_tree,
                            prune_to = prune_to,
                            adjust_branch_length = adjust_branch_length,
                            tax_level_labels = tax_level_labels)


  tree_plot <- do.call(display_subtree,
                       args = c(list(base_tree = pruned_tree,
                               base_name = prune_name,
                               tax_level_labels = tax_level_labels),
                               args))

  return(tree_plot)
}



#' Display taxonomy tree, annotated with numbers of chemicals and fraction
#' overlap between chemicals for two data sets.
#'
#' @param base_tree A \code{\link[ape]{phylo}}-class tree object to plot.
#' @param base_name A name for the base tree (used for the plot title)
#' @param data_1 A \code{data.frame} containing a classified list of entities
#' @param name_1 A name for the list in \code{data_1}, used for the plot legend
#' @param data_2 Another \code{data.frame} containing a classified list of
#'   entities
#' @param name_2 A name for the list in \code{data_2}, used for the plot legend
#' @param entity_id_col The name of the variable identifying unique entities in
#'   both \code{data_1} and \code{data_2}. Must be the same in both data sets.
#' @param group_level Taxonomy level at which to aggregate entities. Default
#'   \code{"terminal"}: Calculate number of entities and overlap for each tip
#'   label. Can also be any element of \code{tax_level_labels}, or an integer
#'   between 1 and \code{length{tax_level_labels}}. In this case, entities will
#'   be grouped by unique labels at the specified taxonomic level, rather than
#'   by terminal (tip) labels, for calculation and plotting of number of
#'   chemicals and overlap.
#' @param tax_level_labels Taxonomy level labels. By default, the
#'   ClassyFire/ChemOnt levels: kingdom, superclass, class, subclass, level5,
#'   ... level11.
#' @param annot_angle Angle at which to plot the annotation text (names of
#'   datasets).
#' @param ... Additional arguments as for \code{\link{display_subtree}}.
#' @return A \code{\link[ggtree]{ggtree}} plot object, branches highlighted by
#'   list membership as in \code{\link{display_subtree}}, with three layers of
#'   heatmap annotation at the tree tips. The innermost layer represents the
#'   number of entities in each tip label (or group of tip labels) in list
#'   \code{data_1}. The second (middle) layer represents the number of entities
#'   in each tip label (or group of tip labels) for the list \code{data_2}. The
#'   outermost layer represents the fraction of overlap between the entities at
#'   each tip label (or group of tip labels), calculated as (size of
#'   intersection)/(size of union).
#' @export
display_overlap <- function(base_tree,
                            base_name,
                            data_1,
                            name_1,
                            data_2,
                            name_2,
                            entity_id_col,
                            group_level = "terminal",
                            tax_level_labels = chemont_tax_levels,
                            annot_angle = "auto",
                            ...){

  args <- list(...)

  if(is.numeric(group_level)){
    if(group_level > length(tax_level_labels)){
      stop(
        paste("Cannot find overlap at level  'group_level' =",
              group_level,
              "because it is greater than the max level",
              "defined by the length of 'tax_level_labels' =",
              paste(tax_level_labels, collapse = ", ")
        )
      )
    }else{
      #pull the corresponding taxonomy level label
      group_level <- tax_level_labels[group_level]
    }
  }

overlap <- calc_number_overlap(data_1 = data_1,
                               data_2 = data_2,
                               entity_id_col = entity_id_col,
                               at_level = group_level,
                               tax_level_labels = tax_level_labels)

overlap$n_1[overlap$n_1==0] <- NA_real_
overlap$n_2[overlap$n_2==0] <- NA_real_

if(group_level %in% "terminal"){
  group_level <- "terminal_label"
}else{
#get tip labels associated with each label in overlap,
#if group_level is not already terminal
#geom_fruit only works on tip labels apparently

  overlap$at_node <- get_node_from_label(label = overlap[[group_level]],
                                         tree = base_tree)

  overlap <- overlap[!is.na(overlap$at_node), ] #NA for any label not in the base tree
  #get children nodes for each
  overlap_tip_nodes <- phangorn::Descendants(x = base_tree,
                                                      node = overlap$at_node,
                                                      type = "tips")

  #now repeat the rest of the columns for each one
  df_list <- lapply(seq_along(overlap_tip_nodes),
                    function(i){
                      data.frame(tip_nodes = overlap_tip_nodes[[i]],
                                 n_1 = overlap[i, "n_1"],
                                 n_2 = overlap[i, "n_2"],
                                 n_intersect = overlap[i, "n_intersect"],
                                 n_union = overlap[i, "n_union"],
                                 simil = overlap[i, "simil"])
                    })

  overlap <- dplyr::bind_rows(df_list)
  overlap$terminal_label <- get_label_from_node(node = overlap$tip_nodes,
                                               tree = base_tree)

}

#unless otherwise specified, use fan layout with open angle, to allow labels
if(!("layout" %in% names(args))){
  args$layout <- "fan"
  args$base_opts$open.angle <- 30
}

if(annot_angle %in% "auto"){
if(args$layout %in% "fan"){
annot_angle <- 360 - 2*args$base_opts$open.angle
}else if(args$layout %in% c("circular",
                                        "equal_angle",
                                        "daylight")){
  annot_angle <- 270
}else if(args$layout %in% c("rectangular",
                            "roundrect",
                            "slanted",
                            "ellipse")){
annot_angle <- 0
}
}

out_obj <- do.call(display_subtree,
                   args = c(list(base_tree = base_tree,
                base_name = base_name,
                data_1 = data_1,
                name_1 = name_1,
                data_2 = data_2,
                name_2 = name_2,
                tax_level_labels = tax_level_labels,
                clade_level = NULL),
                args[setdiff(names(args),
                             "clade_level")
                     ]
                )
                ) +
  ggtreeExtra::geom_fruit(data = overlap,
             geom = geom_tile,
             mapping = aes(y = terminal_label,
                           x = 5,
                           fill = n_1,
                           height =1,
                           width = 10),
             pwidth = 0.1,
             offset = 0.01) +
  ggtreeExtra::geom_fruit(data = overlap,
             geom = geom_tile,
             mapping = aes(y = terminal_label,
                           x = 5,
                           fill = n_2,
                           height =1,
                           width = 10),
             pwidth = 0.1,
             offset = 0.01) +
  ggplot2::scale_fill_viridis_c(trans = "log10",
                       name = "# entities",
                       na.value = "white") +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(data = overlap,
             geom = geom_tile,
             mapping = aes(y = terminal_label,
                           x = 5,
                           fill = simil,
                           height = 1,
                           width = 10),
             pwidth = 0.1,
             offset = 0.01) +
  ggplot2::scale_fill_viridis_c(option = "magma") +
  ggplot2::annotate(geom = "text",
                    x = 110,
                    y = 0,
                    label = name_1,
                    angle = annot_angle,
                    hjust = "left") +
  ggplot2::annotate(geom = "text",
                    x = 120,
                    y = 0,
                    label = name_2,
                    angle = annot_angle,
                    hjust = "left") +
  ggplot2::annotate(geom = "text",
                    x = 130,
                    y = 0,
                    label = "overlap",
                    angle = annot_angle,
                    hjust = "left") +
  ggplot2::theme(legend.position = "left")

if(!is.null(args$clade_level)){
  #first plot arcs without labels, no offset
  clade_opts_tmp <- args$clade_opts
  clade_opts_tmp$offset <- 0
  clade_opts_tmp$textcolour <- NA
  out_obj <- add_cladelab(tree_plot = out_obj,
                          tree = base_tree,
                          clade_level = args$clade_level,
                          clade_opts = clade_opts_tmp)

  #now add arcs with labels, offset of 40 to be outside of geom_fruit
  clade_opts_tmp <- args$clade_opts
  clade_opts_tmp$offset <- 40
 out_obj <- add_cladelab(tree_plot = out_obj,
                         tree = base_tree,
                         clade_level = args$clade_level,
                         clade_opts = clade_opts_tmp)
}

return(out_obj)
}


#'Add clade labels to an existing \code{\link[ggtree]{ggtree}} plot.
#'

#'Wrapper for \code{\link[ggtree]{geom_cladelab}} to add clade labels to an
#'existing \code{\link[ggtree]{ggtree}} plot. Handles some fiddly data wrangling
#'and allows the option to plot clade labels with alternating thick/thin bar
#'widths, which is useful to distinguish clades when they are separated by
#'little or no space on the plot.
#'
#'Currently, \code{add_cladelab} cannot be chained using the
#'\code{\link{[ggplot2]{`+.gg`}}} operator. This is because \code{add_cladelab}
#'is not currently defined as an S3 class with an associated \code{ggplot_add}
#'method.
#'
#'#'For example, the following code will *not* work:
#'
#'```R ggtree(tree) + add_cladelab() ```
#'
#'\code{add_cladelab} also cannot be chained using the
#'\code{\link[magrittr]{%>%}} operator if the preceding chain involves
#'\code{`+`}. This is because of operator precedence: R evaluates
#'\code{\link[magrittr]{%>%}} before \code{\link{[ggplot2]{`+.gg`}}}.
#'
#'
#'For example, the following code also will *not* work:
#'
#'```R ggtree(tree) + layout_circular() %>% add_cladelab() ```
#'
#'However, the following code *will* work:
#'
#'```R ggtree(tree) %>% add_cladelab() ```

#'
#'# Clade options
#'
#'The \code{cladeopt} parameter may include any of the "additional parameters"
#'understood by \code{\link[ggtree]{geom_cladelab}}, including
#'
#'* offset
#'* offset.text
#'* align
#'* extend
#'* angle
#' * horizontal
#' * barsize
#' * barcolour
#' * fontsize
#'* textcolour
#'* imagesize
#' * imagecolour
#'
#'Additionally, in \code{add_cladelab}, the following parameters are
#'understood:
#'
#' * \code{wrap}: the number of characters to wrap text labels. Default 20.
#'    Use NULL for no wrapping.
#' * \code{barsize = "alternate"}: Alternates thick and thin bars. This is the
#'    default. It is useful to distinguish tightly-packed clades.
#' * \code{default_to_tip} Logical: Whether to default to showing tip labels if
#'    no clade label exists at the specified level. For example, if a tip
#'     terminates at level 3, but \code{clade_level = 4}, should that tip be
#'     labeled with its tip label (\code{default_to_tip = TRUE}), or should
#'     that tip be unlabeled (\code{default_to_tip = FALSE})? The default
#'     value is \code{TRUE}.
#'
#'
#' @param tree_plot A \code{\link[ggtree]{ggtree}} plot object, e.g. the output
#'  of \code{\link{display_subtree}} or of \code{\link[ggtree]{ggtree}}.
#' @param tree Optional: the \code{\link[ape]{phylo}} tree object plotted in
#'  \code{tree_plot}. Function will run faster if \code{tree} is provided;
#'  otherwise \code{tree} will be reconstructed from \code{tree_plot$data} using
#'  \code{\link{generate_taxonomy_tree}}.
#' @param clade_level Numeric or \code{"auto"}: the taxonomy level at which to
#'  label clades, where level 0 is the root. Default value \code{"auto"} selects
#'  the taxonomy level of the MRCA of the tips plus one, or level 2, whichever
#'  is the greater number. If \code{clade_level = NULL}, no clade labels will be
#'  added and \code{tree_plot} will be returned unchanged.
#'@param clade_opts A named list of parameters controlling the appearance of
#'  clade labels. See "Details."
#'@return \code{tree_plot} with clade labels added.
add_cladelab <- function(tree_plot,
                          tree = NULL,
                          clade_level = "auto",
                          clade_opts = list(wrap = 20,
                                            barsize = "alternate",
                                            fontsize = 3,
                                            lineheight = 0.7,
                                            default_to_tip = TRUE))
{
  if(!is.null(clade_level)){
  clade_opts_default <- list(wrap = 20,
                             barsize = "alternate",
                             fontsize = 3,
                             lineheight = 0.7,
                             default_to_tip = TRUE)
 #keep defaults for any clade_opts not specified
  clade_opts <- c(clade_opts,
                      clade_opts_default[setdiff(names(clade_opts_default),
                                                 names(clade_opts))])

  if(is.null(tree)){
    #get tree from tree_plot$data
    tmp_df <- tree_plot$data[, c("label", "node", "parent")]
    tmp_df <- setNames(tmp_df, c("Name", "ID", "Parent_ID"))
    #root node has itself as parent -- fix that
    tmp_df[tmp_df$ID == tmp_df$Parent_ID, "Parent_ID"] <- NA_real_
    tmp_df <- as.data.frame(tmp_df)
    tree <- generate_taxonomy_tree(tmp_df)
    rm(tmp_df)
  }


    dat <- get_tree_df(tree)
    if(clade_level %in% "auto"){
      #automatically set clade level to be level of tree's MRCA plus one,
      #or at least 2, but not more than the tree's max level
      mrca_tree <- ape::getMRCA(tree, 1:ape::Ntip(tree))
      mrca_level <- dat[dat$node %in% mrca_tree, "level"]
      max_tree_level <- max(dat$level)
      clade_level <- min(max(mrca_level + 1, 2), max_tree_level)
    }

    dat2 <- dat[dat$level == clade_level, c("node", "Name")]

    #sensible default for clade label angles
    #for circular-ish layouts, use angle = "auto"
    #for rectangular-ish layouts, use angle = 0

    #get layout from tree_plot
    layout <- eval(expression(layout),
                   envir = tree_plot$plot_env)

    if(!("angle" %in% names(clade_opts))){
      if (layout %in% c("rectangular",
                        "roundrect",
                        "slanted",
                        "ellipse")){
        clade_opts$angle <- 0
      }else if(layout %in% c("circular",
                             "fan",
                             "equal_angle",
                             "daylight")){
        clade_opts$angle <- "auto"
      }
    }


  #plot clade bars with alternating widths
  #to do this:
  #first need to get order in which clades are plotted
  #start with order in which *tips* are plotted
  #ggtree:get_taxa_name() gives us tips in plotting order
  tips_plot <- ggtree::get_taxa_name(tree_view = tree_plot)
  #get clade label corresponding to each of these tips at the specified level
  clade_plot <- get_clade(node = get_node_from_label(label = tips_plot,
                                                     tree = tree),
                          tree = tree,
                          level = clade_level)

  if(isTRUE(clade_opts$default_to_tip)){
  #if there is no clade at the specified level, i.e. branch terminates before that level,
  #then label as tip
  clade_plot[is.na(clade_plot)] <- get_node_from_label(label = tips_plot[is.na(clade_plot)],
                                                       tree = tree)
  }else{
    #remove any NA values from clade_plot -- these will not be labeled
    clade_plot <- clade_plot[!is.na(clade_plot)]
  }
  #clade_plot <- clade_plot[!is.na(clade_plot)]
  #tips_clade gives the clades in order of plotting
  #assign alternating bar widths in plotting order
  clade_dat <- data.frame(phylo_node = unique(clade_plot),
                          clade_name = get_label_from_node(node = unique(clade_plot),
                                                           tree = tree),
                          barsize = rep(1:2,
                                        length.out = length(
                                          unique(clade_plot)
                                        )
                          )
  )

  #Alternating bar thickness only works if the dataset is in base tree order
  #This is because geom_cladelab() can't use aes() mapping for barsize
  #So we have to supply it as a non-mapped argument outside aes()
  #and ggtree automatically reorders it to match the base tree
  #which means, if it's already in plotting order,
  #it gets reordered wrongly.
  #so it needs to match the order of the base tree to begin with,
  #so that the auto-reordering will be correct.

  dat3 <- clade_dat[match(intersect(dat$node,
                                    clade_dat$phylo_node),
                          clade_dat$phylo_node), ]


  if(clade_opts$barsize %in% "alternate"){
    #Re-sort the alternating bar sizes in clade plotting order
    #to correspond to base tree order of clades
    clade_opts$barsize <- dat3$barsize
  }


  if("wrap" %in% names(clade_opts)){
    #wrap clade names to have width clade_label_wrap characters
    dat3$clade_name2 <- stringr::str_wrap(dat3$clade_name,
                                          clade_opts$wrap)
  }


  tree_plot +
    do.call(geom_cladelab,
            args = c(list(data = dat3,
                          mapping = aes(node = phylo_node,
                                        label = clade_name2,
                                        group = clade_name2)),
                     clade_opts))
  }else{
  tree_plot
  }
}


#' Side by side trees
#'
#' This function takes two data sets and construsts both trees and internal
#' layers for data visualization. The left and right out layer of the diagram
#' consists of subtrees corresponding to the left and right data sets. The left
#' and right inner layers display label statistics for the left and right data
#' sets. The center displays label statistics shared by the left and right data
#' sets.
#'
#' @param data_left The left data.table of chemical classifications.
#' @param data_right The right data.table of chemical classifications.
#' @param name_left Alternate parameter for left data set.
#' @param name_right Alternate parameter for right data set.
#' @param log_trans Alternate parameter for log-transforming data.
#' @return An `aplot` consisting of two outer layer `ggtree` objects and three
#'   inner layer `ggplot2` objects.
#' @export
side_by_side_trees <- function(data_left, data_right, name_left = 'Left tree', name_right = 'Right tree', log_trans = FALSE){
  if(TRUE) {

    if (!all('terminal_label' %in% names(data_left))){
      warning('The column `terminal_label` is missing from one of the first input data.table! Attaching column...')
      data_left <- add_terminal_label(data_left)
      }
    if (!all('terminal_label' %in% names(data_right))){
      warning('The column `terminal_label` is missing from one of the second input data.table! Attaching column...')
      data_right <- add_terminal_label(data_right)
      }

    terminal_label <- NULL
    INCHIKEY <- NULL
    tree <- NULL
    tip.label <- NULL
  terminal_labels_left <- data_left[!is.na(terminal_label), unique(terminal_label)]
  terminal_labels_right <- data_right[!is.na(terminal_label), unique(terminal_label)]

  terminal_labels <- union(terminal_labels_left, terminal_labels_right)

  left_initial_tree <- prune_and_display_subtree(data = data_left, no_plot = TRUE)
  right_initial_tree <- prune_and_display_subtree(data = data_right, no_plot = TRUE)

  left_labels <- c(left_initial_tree$tip.label, left_initial_tree$node.label)
  right_labels <- c(right_initial_tree$tip.label, right_initial_tree$node.label)

  all_labels <- union(left_labels, right_labels)

  union_tree <- drop_tips_nodes(tree = chemont_tree, labels = all_labels)
  union_tree$edge.length <- adjust_branch_lengths(union_tree)
  terminal_labels <- intersect(terminal_labels, union_tree$tip.label)

  left_ancestors <- lapply(left_labels, function(t) {get_ancestors(chemont_tree, t)})
  right_ancestors <- lapply(right_labels, function(t) {get_ancestors(union_tree, t)})

  all_left_tree <- union(unlist(left_ancestors), left_labels)
  all_right_tree <- union(unlist(right_ancestors), right_labels)

  #print(all_left_tree)
  #print(all_right_tree)

  # Adjust the tip label size depending on the number of tips of the tree
  #if (length(union_tree$tip.label) <= 200){
  #  tip_size = 3
  #} else if (length(tree_1$tip.label) <= 500){
  #  tip_size = 1.5
  #} else {
  #  tip_size = .5
  #}

  tip_size <- 2/length(union_tree$tip.label)

  left_tree <- ggtree(union_tree,
                      aes(color= (c(union_tree$tip.label, union_tree$node.label) %in% all_left_tree)),
                      branch.length = FALSE)+
    scale_color_manual(values = c('black', 'blue'),
                       labels = c('', name_left),
                       name = 'Left Tree') + geom_tippoint(size = tip_size)
  right_tree <- ggtree(union_tree,
                       aes(color= c(union_tree$tip.label, union_tree$node.label) %in% all_right_tree),
                       branch.length = FALSE) +
    scale_color_manual(values = c('black', 'red'),
                       labels = c('', name_right),
                       name = 'Right Tree') + geom_tippoint(size = tip_size) +
    geom_tiplab(as_ylab=TRUE, size = 3)
  right_tree <- right_tree + ggplot2::scale_x_continuous(trans = "reverse")

  if(FALSE){# Add clade labels to each tree.
    left_tree <- add_cladelab(left_tree)
    right_tree <- add_cladelab(right_tree)
  }

  nTip <- length(union_tree$tip.label)

  tree_data <- data.frame(tip.label = rep(union_tree$tip.label, 3),
                          tree = rep(c('left', 'right', 'center'), each = nTip))#,
                          #left = logical(length(union_tree$tip.label)),
                          #right = logical(length(union_tree$tip.label)))
  #tree_data$left <- tree_data$tip.label %in% terminal_labels_left
  #tree_data$right <- tree_data$tip.label %in% terminal_labels_right

  tree_data$tree <- factor(tree_data$tree, levels = c('left', 'center', 'right'))

  tree_data_leftval <- double(nTip)
  tree_data_centerval <- double(nTip)
  tree_data_rightval <- double(nTip)

  for (i in seq_along(union_tree$tip.label)){
    t <- tree_data$tip.label[[i]]
    data_left_chemicals <- data_left[terminal_label == t, unique(INCHIKEY)]
    data_right_chemicals <- data_right[terminal_label == t, unique(INCHIKEY)]
    shared_chemicals <- intersect(data_left_chemicals, data_right_chemicals)
    tree_data_leftval[[i]] <- length(data_left_chemicals)
    tree_data_rightval[[i]] <- length(data_right_chemicals)
    tree_data_centerval[[i]] <- length(shared_chemicals)
  }

  #print(tree_data_leftval)
  #print(tree_data_rightval)
  #print(tree_data_centerval)
  #print(nTip)

  value <- c(tree_data_leftval, tree_data_rightval, tree_data_centerval)
  tree_data <- cbind(tree_data, value)
  names(tree_data)[[3]] <- "value"

  #print(names(tree_data))

  #print(str(tree_data))
  #left_data_plot <- ggplot(tree_data, aes(y = tip.label, x = leftval)) +
  #  geom_tile()#Or some other data viz
  #return(left_data_plot)
  #center_data_plot <- ggplot(tree_data, aes(y = tip.label, x = centerval)) +
  #  geom_tile()#Or some other data viz
  #right_data_plot <- ggplot(tree_data, aes(y = tip.label, x = rightval)) +
  #  geom_tile()#Or some other data viz

  #return(tree_data)

  trans <- ifelse(log_trans, 'log1p', 'identity')

  data_plot <- ggplot(tree_data, aes(x = tree, y = tip.label)) +
    geom_tile(aes(fill = value)) + scale_fill_viridis(trans=trans) +
    theme_minimal() + ylab(NULL)  +
      theme(axis.text.y = element_text(size = 3),
            axis.title.y = NULL)
  #return(data_plot)

  #data_plot <- center_data_plot %>% aplot::insert_left(left_data_plot)
  #data_plot <- data_plot %>% aplot::insert_right(right_data_plot)

  data_plot <- data_plot %>% aplot::insert_left(left_tree)

  #return(data_plot)

  data_plot$n <- 3
  data_plot_new_col <- matrix(3, nrow = 1)

  data_plot$width <- c(data_plot$width, 1)
  data_plot$layout <- cbind(data_plot$layout, data_plot_new_col)
  data_plot_axis <- list(ylab(data_plot$plotlist[[2]]$labels$y))
  data_plot$plotlist[[2]] <- data_plot$plotlist[[2]] + data_plot_axis
  data_plot$plotlist[[3]] = right_tree




  return(data_plot)
  }


}


