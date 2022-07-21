



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
    dplyr::rename(terminal_label = label)

  data <- merge(data, labels, by = setdiff(names(data),
                                           tax_level_labels))

  return(data)
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

#' Get subtree node numbers
#'
#' This is a helper function that takes a classified set and a taxonomy tree and
#' provides all node numbers in the subtree corresponding to the data set.
#'
#' @param data A classified data set.
#' @param tree A taxonomy tree object, as returned by
#'   \code{\link{generate_taxonomy_tree}}.
#' @param tax_level_labels A vector of levels for the taxonomy.
#' @return A vector of node numbers in the base tree that are represented in the
#'   subtree corresponding to the data set.

get_subtree_nodes <- function(data,
                              base_tree = chemont_tree,
                              tax_level_labels = c('kingdom', 'superclass',
                                                'class', 'subclass',
                                                'level5', 'level6',
                                                'level7', 'level8',
                                                'level9', 'level10',
                                                'level11')){
  #get labels represented by the classified dataset
  data_labels <- get_labels(data = data,
                            tax_level_labels = tax_level_labels)

  #get node numbers, levels, & names of base tree
  tree_df <- get_tree_df(base_tree)

  #get node numbers represented in dataset
  data_nodes <- tree_df[tree_df$Name %in% data_labels, "node"]
  #get all ancestors of nodes in input data set,
  #plus the nodes themselves
  data_all_nodes <- unique(c(unlist(phangorn::Ancestors(x = base_tree,
                                                    node = data_nodes)),
                         data_nodes))
  #return the vector of node numbers in the subtree
  return(data_all_nodes)

}

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
label_bars <- function(data = NULL, tax_level_labels = NULL){
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

  df <- data.frame(tax_level_labels, unname(sapply(data, function(t) get_label_length(get_labels(t, tax_level_labels)))))
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
#'This function takes in one or two data.tables consisting of chemicals with
#'classification data and returns a tree diagram indiacting the subtree(s)
#'induced by the data.table(s).
#'
#'@param base_tree The "base tree" to plot, as a \code{\link[ape]{phylo}}-class
#'  object. Branches of this base tree will be color-coded to indicate subtree
#'  membership. Default is the full ChemOnt taxonomy,
#'  \code{\link{chemont_tree}}.
#'@param base_name Will be used as the plot title. Usually this should name or
#'  describe the base tree. Default is NULL, in which case the *name* of the
#'  variable passed to \code{base_tree} will be used as the title (with "tree" appended).
#'@param data_1 A data.frame consisting of a list of chemicals with
#'  classification data.
#'@param data_2 Optional: A second data.frame consisting of a list of chemicals
#'  with classification data.
#'@param name_1 A string giving the name of the first data.frame, for plot
#'  labeling. Default is "Set 1".
#'@param name_2 A string giving the name of the second data.frame (if there is
#'  one), for plot labeling. Default is "Set 2".
#'@param tax_level_labels An alternate parameter giving the taxonomy levels if
#'  not using ClassyFire taxonomy.
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
#'@param clade_label_level The taxonomy level at which to draw clade labels, if
#'  any. Root is level 0. Default is level 2 (superclass, in ChemOnt taxonomy).
#'  Set to NULL to suppress clade labels altogether.
#'@param clade_label_fontsize Font size for clade labels. Default 3.
#'@param clade_label_wrap Number of characters at which to wrap clade labels to
#'  a second line. Default 20.
#'@param clade_label_lineheight Line height for multi-line clade labels, as a
#'  multiple of the size of text. Controls vertical space between lines on
#'  multi-line clade labels. Default 0.7.
#'@return A ggtree object visualizing the full base tree, with branches
#'  color-coded to indicate whether they are present in \code{data_1}, in
#'  \code{data_2} if supplied, neither, or both.
#'@export
#'@import phangorn
#'@import ggtree
display_subtree <- function(base_tree = chemont_tree,
                            base_name = NULL,
                            data_1 = NULL,
                            data_2 = NULL,
                            name_1 = NULL,
                            name_2 = NULL,
                            tax_level_labels = c('kingdom', 'superclass',
                                                 'class', 'subclass',
                                                 'level5', 'level6',
                                                 'level7', 'level8',
                                                 'level9', 'level10',
                                                 'level11'),
                            layout = "circular",
                            base_opts = list("color" = "gray80",
                                             "size" = 0.5,
                                             "linetype" = 1),
                            subtree_mapping = list("color" = c("gray80",
                                                               "#66C2A5",
                                                               "#8DA0CB",
                                                               "#FC8D62")),
                            clade_level = 2,
                            clade_opts = list(wrap = 20,
                                              barsize = "alternate",
                                              fontsize = 3,
                                              lineheight = 0.7)
){

  if(is.null(base_name)){
    #get the *name* of the variable that was passed to base_tree and append "tree"
    base_name <- as.character(substitute(base_tree))
  }

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
    }





  }



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
    if(is.null(name_1)){
      #get the name of the variable passed to data_1
      name_1 <- paste(as.character(substitute(data_1)))
    }
  #Get node numbers of data_1 subtree
  data_1_all <- get_subtree_nodes(data = data_1,
                                  base_tree = base_tree,
                                  tax_level_labels = tax_level_labels)

  #Data frame with all node numbers in taxonomy tree
   cohort_data <- get_tree_df(tree = base_tree)
   #Categorical column: is each node in Data Set 1?
   #0 = no, 1 = yes
   cohort_data$inSet1 <- ifelse(cohort_data$node %in% data_1_all,
                                1L,
                                0L)


  if (!is.null(data_2)) {
    if(is.null(name_2)){
      #get the name of the variable passed to data_1
      name_2 <- paste(as.character(substitute(data_2)))
    }
#Get node numbers of Data Set 2 subtree
    data_2_all <- get_subtree_nodes(data = data_2,
                                    base_tree = base_tree,
                                    tax_level_labels = tax_level_labels)
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
    cohort_data$List <- factor(cohort_data$inSet1 + cohort_data$inSet2,
                                 levels = 0:3,
                                 labels = c("Neither set",
                                            paste(name_1, "only"),
                                            paste(name_2, "only"),
                                            paste(name_1, "and", name_2)))

  }else{
    #Categorical column: Is each node in Data Set 1 or not?
    cohort_data$List <- factor(cohort_data$inSet1,
                               levels = 0:1,
                               labels = c(paste("Not in", name_1),
                                          paste("In", name_1)))
  }

   #Name the subplot_mapping items after the categories in cohort_data$List
   subtree_mapping <- sapply(subtree_mapping,
                             function(x) setNames(x, levels(cohort_data$List)),
                             simplify = FALSE,
                             USE.NAMES = TRUE)
   #set up for aes call -- list of aesthetic mappings, all applied to "List" in cohort_data
   subtree_aes <- replicate(n= length(subtree_mapping),
                             expr = quote(List))
   #name the list elements after the aesthetics in subtree_mapping
   subtree_aes <- setNames(subtree_aes, names(subtree_mapping))
   #you end up with something like: list(color = quote(List), size = quote(List))
   #with do.call(aes, subtree_aes), it's equivalent to calling
   #aes(color = List, size = List)

   #prepare a list of manual scales as provided in subtree_mapping
   #do this using scale_discrete_manual() as follows:
   #specify aesthetic, like "color" or "size"
   #then specify the values for the manual scale,
   #like c("gray80", "#66C2A5", "#8DA0CB","#FC8D62")
   scale_list <- lapply(names(subtree_mapping),
                        function(aesthetic) {
                          ggplot2::scale_discrete_manual(aesthetics = aesthetic,
                                                         name = "List presence",
                                                         values =  subtree_mapping[[aesthetic]],
                                                         breaks = levels(cohort_data$List),
                                                         limits = levels(cohort_data$List))
                        }
   )

#add list presence data
   tree_plot <- tree_plot %<+% cohort_data +
     do.call(aes,
             subtree_aes) +
      scale_list +
     theme(legend.position = "top")
  } #end if(!is.null(data_1))

   #if clade labels have been selected
    if(!is.null(clade_level)){
      dat <- get_tree_df(base_tree)
      dat2 <- dat[dat$level == clade_level,]
      dat2 <- setNames(dat2, c("phylo_node", "level", "clade_name"))

      if("wrap" %in% names(clade_opts)){
      #wrap clade names to have width clade_label_wrap characters
      dat2$clade_name2 <- stringr::str_wrap(dat2$clade_name,
                                            clade_opts$wrap)
      }

      #sensible default for clade label angles
      #for circular-ish layouts, use angle = "auto"
      #for rectangular-ish layouts, use angle = 0
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
      if(clade_opts$barsize %in% "alternate"){
        #plot clade bars with alternating widths
        #to do this:
        #first need to get order in which clades are plotted
        #start with order in which *tips* are plotted
        #ggtree:get_taxa_name() gives us tips in plotting order
        tips_plot <- ggtree::get_taxa_name(tree_view = tree_plot)
        #get clade label corresponding to each of these tips at the specified level
        clade_plot <- get_clade(node = get_node_from_label(label = tips_plot,
                                                           tree = base_tree),
                                tree = base_tree,
                                level = clade_level)
        #if there is no clade at the specified level, i.e. branch terminates before that level,
        #then label as tip
        clade_plot[is.na(clade_plot)] <- get_node_from_label(label = tips_plot[is.na(clade_plot)],
                                                             tree = base_tree)
        clade_plot <- clade_plot[!is.na(clade_plot)]
        #tips_clade gives the clades in order of plotting
        #assign alternating bar widths in plotting order
        clade_dat <- data.frame(phylo_node = unique(clade_plot),
                                clade_name = get_label_from_node(node = unique(clade_plot),
                                                                 tree = base_tree),
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

        #Re-sort the alternating bar sizes in clade plotting order
        #to correspond to base tree order of clades
        barsize_vect <- clade_dat$barsize[match(dat2$clade_name,
                                                clade_dat$clade_name)]
        clade_opts$barsize <- barsize_vect
      }

      tree_plot <- tree_plot +
        do.call(geom_cladelab,
                args = c(list(data = dat2,
                              mapping = aes(node = phylo_node,
                                            label = clade_name2,
                                            group = clade_name2)),
                         clade_opts)) +
        theme(legend.position = "top")
    }

  #add title of base tree
  tree_plot <- tree_plot + ggtitle(base_name)

  return(tree_plot)

}

#'Prune a tree
#'
#'Prune a tree to keep only a subtree specified as a classified data set, a
#'vector of labels, or a vector of node numbers
#'
#'@param tree The "base tree" to be pruned.
#'@param prune_to What to keep. May be a \code{data.frame} of classified data;
#'  one or more labels in the tree (tip or internal node labels); one or more
#'  node numbers in the tree (tip or internal nodes). See Details. Default is
#'  NULL, which results in no pruning being done.
#'@param adjust_branch_length Whether to adjust branch length so that all
#'  newly-pruned terminal nodes appear at the same length as tips, even if they
#'  were originally internal nodes. Default FALSE.
#'@param tax_level_labels Vector of the possible taxonomy levels that can appear
#'  as column names in \code{prune_to} if it is a \code{data.frame} of
#'  classified data.
#'@return A \code{\link[ape]{phylo}}-class object representing the pruned tree.
#'@export
prune_tree <- function(tree,
                       prune_to = NULL,
                       adjust_branch_length = FALSE,
                       tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                            'level5', 'level6', 'level7', 'level8',
                                            'level9', 'level10', 'level11')){
  if(!is.null(prune_to)){ #if user has specified something to prune to
    if(is.data.frame(prune_to)){ #if user has specified a dataset to prune to
      #Prune the tree according to the specified dataset
      pruned_tree <- drop_tips_nodes(tree = tree,
                                     data = prune_to,
                                     tax_level_labels = tax_level_labels)
    }else if(is.character(prune_to)){
      #interpret as node/tip labels
      #prune to only the subtree including the descendents of this internal node(s)
      pruned_tree <- drop_tips_nodes(tree = tree,
                                     labels = prune_to)
    }else if(is.numeric(prune_to)){
      #interpret as node numbers
      #prune to only the subtree including the descendents of this internal node(s)
      pruned_tree <- drop_tips_nodes(tree = tree,
                                     nodes = prune_to)
    }

    if (adjust_branch_length) {
      pruned_tree$edge.length <- adjust_branch_lengths(pruned_tree)
    }

  }
  return(pruned_tree)
}

#' Prune and display subtree
#'
#' This function takes a data.table of chemicals and their classification data,
#' and returns the data-induced subtree with an optional tree diagram of the
#' subtree.
#'
#' @param data A data.table with classification data for chemicals.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @param show_tips An alternate parameter determining whether tip labels are
#'   displayed.
#' @param no_plot An alternate parameter for returning only the pruned tree
#'   without the tree visual.
#' @param adjust_branch_length An alternate parameter determining whether to
#'   resize branches of subtree.
#' @param xlimit An alternate parameter determining the scale of the tree
#'   visualization.
#' @return A pruned tree or list consisting of a pruned tree and ggtree diagram
#'   of the pruned tree.
#' @export
#' @importFrom ape drop.tip
#' @import ggtree
prune_and_display_subtree <- function(base_tree = chemont_tree,
                                      base_name = NULL,
                                      prune_to = NULL,
                                      prune_name = NULL,
                                      adjust_branch_length = FALSE,
                                      tax_level_labels = c('kingdom',
                                                           'superclass',
                                                           'class',
                                                           'subclass',
                                                           'level5',
                                                           'level6',
                                                           'level7',
                                                           'level8',
                                                           'level9',
                                                           'level10',
                                                           'level11'),
                                      ...) { #args as for display_subtree()

  args <- list(...)

  pruned_tree <- prune_tree(tree = base_tree,
                            prune_to = prune_to,
                            adjust_branch_length = adjust_branch_length,
                            tax_level_labels = tax_level_labels)

  #get prune_name if NULL
  if(is.null(prune_name)){
    if(is.null(prune_to)){
      if(!is.null(base_name)){
        prune_name <- base_name
      }else{
        prune_name <-as.character(substitute(base_tree))
      }
    }else{
      if(is.data.frame(prune_to)){
        prune_name <- as.character(substitute(prune_to))
      }else if(is.character(prune_to) | is.numeric(prune_to)){
        prune_name <- paste(prune_to, collapse = ", ")
      }
    }
  }

  tree_plot <- do.call(display_subtree,
                       args = c(list(base_tree = pruned_tree,
                               base_name = prune_name,
                               tax_level_labels = tax_level_labels),
                               args))

  return(tree_plot)
}

#' Circular tree with boxplots
#'
#' This function takes in a data.table of chemicals with classification data and
#' additional numeric data, and displays selected numeric data grouped by tip
#' label on the data-induced subtree of the classification taxonomy.
#'
#' @param data A data.table consisting of classification data and additional
#'   numeric data.
#' @param col A specified numeric column for use in displaying boxplots.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @param title An alternate parameter for the title of the plot.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @param layers An alternate parameter giving which taxonomic layers to display
#'   outside of the boxplot layer. This can be either a string with a single
#'   column name, a vector of column names, or a list of column names. If the
#'   input is a vector or a list, it is fine for it to be length 1.
#' @param tippoint_boxplot Alternate parameter for determining whether to color
#'   the tippoints and the boxplots. If TRUE, they will match in colors, and if
#'   FALSE, the boxplots will be filled white with black outline and the
#'   tippoints will be balck.
#' @return A ggtree object consisting of subtree induced by data and boxplots
#'   corresponding to specified numeric column of data.
#' @export
#' @import ggtree
#' @import ggtreeExtra
#' @import viridis
#' @importFrom ggnewscale new_scale_fill
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
circ_tree_boxplot <- function(data, col, tax_level_labels = NULL, title = NULL, tree = chemont_tree, layers = NULL, tippoint_boxplot = FALSE){
  val <- NULL
  terminal_label <- NULL
  grp <- NULL
  Label <- NULL
  ID <- NULL

  if (!(col %in% names(data)))
    stop(paste('The column', col, 'is not in the input data!'))

  #print(col)
  if (!('terminal_label' %in% names(data))){
    data <- add_terminal_label(copy(data))
  }
  columns <- which(names(data) %in% c(col, 'terminal_label'))
  #print(columns)

  new_data <- copy(as.data.frame(data[, .SD, .SDcols = columns]))
  index <- which(names(new_data) %in% col)
  new_data[, index] <- as.numeric(new_data[, index])
  new_data <- new_data[!is.na(new_data[, index]),]
  new_data[["grp"]] <- new_data$terminal_label
  new_data[["val"]] <- new_data[, index]
  new_data[["node"]] <- new_data$terminal_label

  #summary(new_data)

  new_data_tree <- prune_and_display_subtree(data, tax_level_labels, tree, show_tips = FALSE, no_plot = TRUE)
  #print(new_data_tree)
  num_nodes_tips <- length(new_data_tree$tip.label) + new_data_tree$Nnode
  #print(num_nodes_tips)

  tip_node_data <- data.frame('ID' = c(new_data_tree$tip.label, new_data_tree$node.label),
                              'Label' = c(new_data_tree$tip.label, new_data_tree$node.label))
  tip_node_data$Label <- factor(tip_node_data$Label)

  circ_plot <- ggtree(new_data_tree,
                      layout = 'circular')
  if (tippoint_boxplot) {
    circ_plot <- circ_plot %<+% tip_node_data +
    geom_tippoint(aes(color = Label), show.legend = FALSE) +
      scale_color_viridis(name = 'Terminal label',
                          option = 'magma',
                          discrete = TRUE) + ggnewscale::new_scale_fill()

    circ_plot <- circ_plot + ggtreeExtra::geom_fruit(data = new_data, geom = geom_boxplot,
                                                     mapping = aes(x = val,
                                                                   y = terminal_label,
                                                                   fill = grp),
                                                     size = 0.2,
                                                     outlier.size = 0.5,
                                                     outlier.stroke = 0.08,
                                                     outlier.shape = 21,
                                                     axis.params = list(axis = 'x',
                                                                        text.size = 1.8,
                                                                        text.angle = 270,
                                                                        hjust = 0),
                                                     grid.params = list(),
                                                     show.legend = FALSE) +
      scale_fill_viridis(name = 'Terminal label',
                         option = 'magma',
                         discrete = TRUE) + new_scale_fill()
  } else {
    circ_plot <- circ_plot + geom_tippoint()
    circ_plot <- circ_plot + ggtreeExtra::geom_fruit(data = new_data, geom = geom_boxplot,
                                                     mapping = aes(x = val,
                                                                   y = terminal_label),
                                                     size = 0.2,
                                                     outlier.size = 0.5,
                                                     outlier.stroke = 0.08,
                                                     outlier.shape = 21,
                                                     axis.params = list(axis = 'x',
                                                                        text.size = 1.8,
                                                                        text.angle = 270,
                                                                        hjust = 0),
                                                     grid.params = list(),
                                                     show.legend = FALSE) + new_scale_fill()
  }

  #circ_plot <- circ_plot + ggtreeExtra::geom_fruit(data = new_data, geom = geom_boxplot,
  #                                                 mapping = aes(x = val,
  #                                                               y = terminal_label,
  #                                                               fill = grp),
  #                                                 size = 0.2,
  #                                                 outlier.size = 0.5,
  #                                                 outlier.stroke = 0.08,
  #                                                 outlier.shape = 21,
  #                                                 axis.params = list(axis = 'x',
  #                                                                    text.size = 1.8,
  #                                                                    text.angle = 270,
  #                                                                    hjust = 0),
  #                                                 grid.params = list(),
  #                                                 show.legend = FALSE) + new_scale_fill()
  #circ_plot <- circ_plot + scale_fill_discrete(guide = 'none')
  #circ_plot <- circ_plot + scale_fill_discrete(name = 'Tip label',
  #                                             guide = guide_legend(keywidth = 0.2,
  #                                                                  keyheight = 0.2,
  #                                                                  ncol = 2))
  #

  if (!is.null(layers)){
    label_levels <- get_labels(data, tax_level_labels = tax_level_labels)
    if (is.list(layers) | is.vector(layers)){
    level_names <- names(label_levels)[which(names(label_levels) %in% layers)]
    } else {
      warning('The `layers` parameter must be a list or a vector! Skipping extra layers for now...')
      level_names <- c()
    }
    #print(which(names(label_levels) %in% layers))
    #print(names(label_levels))
    #print(level_names)

    fruit_data <- data.frame('ID' = c(new_data_tree$tip.label, new_data_tree$node.label))
    palettes <- c('Blues', 'Oranges', 'BuGn', 'OrRd', 'BuPu', 'Reds','GnBu', 'RdPu','Greens', 'YlOrBr',
                  'PuBu', 'YlOrRd', 'PuBuGn', 'YlGnBu', 'PuRd', 'YlGn', 'Purples', 'Greys')

    for (i in rev(seq_along(level_names))){
      values <- unname(as.list(data[, unique(.SD), .SDcol = level_names[[i]]]))[[1]]
      values <- values[!is.na(values)]
      values <- values[-which(sapply(values, function(t) {t == ''}))]
      #print(which(sapply(values, function(t) {t == ''})))
      #print(values)
      tree_nodes <- lapply(c(new_data_tree$tip.label, new_data_tree$node.label), function(x) {x})
      for (j in seq_along(values)){
        #print(values[[j]])
        total_descendants <- c(tree$tip.label, tree$node.label)[c(phangorn::Descendants(tree, which(c(tree$tip.label, tree$node.label) %in% values[[j]]), type = 'all'), which(c(tree$tip.label, tree$node.label) %in% values[[j]]))]
        name_indices <- which(c(new_data_tree$tip.label, new_data_tree$node.label) %in% total_descendants)
        #print(name_indices)
        names(tree_nodes)[name_indices] <- values[[j]]
        #names(tree_nodes)[c(phangorn::Descendants(new_data_tree, which(tree_nodes %in% values[[j]]), type = 'all'), which(tree_nodes %in% values[[j]]))] <- values[[j]]
        #print(which(is.na(names(tree_nodes))))
      }
      level_number <- length(unique(names(tree_nodes)))
      level_labels <- unique(names(tree_nodes))
      names(tree_nodes)[which(is.na(names(tree_nodes)))] <- paste0('_', level_names[[i]])

      fruit_data[[level_names[[i]]]] <- factor(names(tree_nodes))


      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = palettes[[i]]))
      current_palette <- colors(level_number)

      circ_plot <- circ_plot + geom_fruit(data = fruit_data,
                                          geom = geom_tile,
                                          mapping = aes(y = ID, x = .data[[level_names[[i]]]], fill = .data[[level_names[[i]]]]),
                                          width = 3,
                                          pwidth = 0,
                                          color = 'white') +
        scale_fill_manual(values = current_palette,
                          labels = level_labels) +
        ggnewscale::new_scale_fill()

      #print(names(tree_nodes))
      #attr(new_data_tree, paste0('_', level_names[[i]])) <- factor(names(tree_nodes))
      #new_data_tree[[paste0('_', level_names[[i]])]] <- factor(new_data_tree[[paste0('_', level_names[[i]])]])
      #new_data_tree <- groupOTU(new_data_tree, tree_nodes, paste0('_', level_names[[i]]))
    }


  }

  if (is.character(title)){
    circ_plot <- circ_plot + ggtitle(title, subtitle = paste('Boxplots of the data from', col, 'column')) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
  return(circ_plot)
}


#' Leaf fraction subtree
#'
#' This function takes in two data.tables, plots the subtree induced by the
#' first data.table and colors the tips based on the proportion of chemicals from the second data.table
#' that make up the chemicals from the first, grouped by tip (or terminal) label.
#'
#' @param data_1 A data.table of chemicals, classifications, and column `terminal_label`.
#' @param data_2 A data.table of chemicals, classifications, and column `terminal_label`
#' @param name_1 Alternate parameter for name of first data.table.
#' @param name_2 Alternate parameter for name of second data.table.
#' @param show_labels Alternate parameter indicating whether to show tip labels.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @return A ggtree plot.
#' @export
#' @import data.table
#' @import ggtree
#' @import ggtreeExtra
#' @import ggplot2
#'
#' @seealso \code{\link{add_terminal_label}}
#'
leaf_fraction_subtree <- function(data_1, data_2, name_1 = 'data_1', name_2 = 'data_2', show_labels = FALSE, tax_level_labels = NULL, tree = NULL){
  if (!all('terminal_label' %in% names(data_1))){
    warning('The column `terminal_label` is missing from one of the first input data.table! Attaching column...')
    data_1 <- add_terminal_label(data_1)
  }
  if (!all('terminal_label' %in% names(data_2))){
    warning('The column `terminal_label` is missing from one of the second input data.table! Attaching column...')
    data_2 <- add_terminal_label(data_2)  }
  # Find all the terminal_label values from data_1.
  terminal_label <- NULL
  INCHIKEY <- NULL
  percentages <- NULL
  terminal_labels <- data_1[!is.na(terminal_label), unique(terminal_label)]

  # For each terminal_label value, determine the chemicals from data_2 that are
  # also in data_1. This checks using the INCHIKEY of each chemical.
  label_percentages <- sapply(terminal_labels, function(t) {
    data_1_chemicals <- data_1[terminal_label == t, unique(INCHIKEY)]
    data_2_chemicals <- data_2[terminal_label == t, unique(INCHIKEY)]
    shared_chemicals <- intersect(data_1_chemicals, data_2_chemicals)
    return(length(shared_chemicals)/length(data_1_chemicals))
  })

  label_data <- data.frame('label' = terminal_labels,
                           'percentages' = unname(label_percentages),
                           'data_1_numbers' <- unname(sapply(terminal_labels, function(t) {
                             length(data_1[terminal_label == t, unique(INCHIKEY)])
                           })),
                           'data_2_numbers' <- unname(sapply(terminal_labels, function(t) {
                             length(intersect(data_1[terminal_label == t, unique(INCHIKEY)],
                                              data_2[terminal_label == t, unique(INCHIKEY)]))
                           }
                           ))
  )
  names(label_data)[3:4] <- c(paste(name_1, 'label numbers'), paste(name_2, 'label numbers in', name_1))

  #print(label_data)

  data_1_tree <- prune_and_display_subtree(data_1, tax_level_labels = tax_level_labels, tree = tree, no_plot = TRUE)

  if (length(data_1_tree$tip.label) <= 200){
    tip_size = 3
  } else if (length(data_1_tree$tip.label) <= 500){
    tip_size = 1.5
  } else {
    tip_size = .5
  }

  #tree <- full_join(data_1_tree, label_data, by = 'label')

  tree_plot <- ggtree(data_1_tree, layout = 'circular') %<+% label_data
  tree_plot <- tree_plot + ggtitle(paste0(name_1, ' subtree')) + theme(plot.title = element_text(hjust = 0.5))

  tree_plot <- tree_plot + geom_tippoint(aes(color = percentages), size = tip_size)
  tree_plot <- tree_plot +
    scale_color_viridis_c(name = paste0('Percentage of ',name_1, ' chemicals that are ', name_2, ' chemicals'),
                          option = 'plasma')

  if (show_labels){
    tree_plot <- tree_plot + geom_tiplab(aes(color = percentages), size = 1)
  } else {
    tree_plot <- tree_plot + ggtreeExtra::geom_fruit(geom = geom_tile,
                                        mapping = aes(color = percentages),
                                        width = 10,
                                        height = .1,
                                        offset = 0.1)
  }

  return(list(tree_plot, label_data))

}


#' Compare data set subtrees
#'
#' This function takes in two data.tables of chemicals with classification data,
#' and creates two plots each of which show a visual comparison of the subtrees
#' induced by two data sets. The visual comparison consists of branches of one
#' data-induced subtree colored and labeled based on whether they are
#' represented by chemicals from the other data set.
#'
#' @param data_1 A data.table of chemicals and their classifications.
#' @param data_2 A data.table of chemicals and their classifications.
#' @param name_1 An alternate parameter giving the name of the first data set.
#' @param name_2 An alternate parameter giving the name of the second data set.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @param show_tips An alternate parameter determining whether to show tip
#'   labels.
#' @param adjust_branch_length An alternate parameter determining whether to
#'   resize branches of subtree.
#' @return A list of two ggtree objects.
#' @export
#' @import ggtree
data_set_subtrees <- function(data_1,
                              data_2,
                              name_1 = 'data_1',
                              name_2 = 'data_2',
                              tax_level_labels = NULL,
                              tree = NULL,
                              show_tips = TRUE,
                              adjust_branch_length = TRUE){
  membership <- NULL
  # Prune subtrees for each data set.
  tree_1 <- prune_and_display_subtree(data_1,
                                      tax_level_labels = tax_level_labels,
                                      tree = tree,
                                      show_tips = show_tips,
                                      no_plot = TRUE,
                                      adjust_branch_length = adjust_branch_length)
  tree_2 <- prune_and_display_subtree(data_2,
                                      tax_level_labels = tax_level_labels,
                                      tree = tree,
                                      show_tips = show_tips,
                                      no_plot = TRUE,
                                      adjust_branch_length = adjust_branch_length)

  # Get all taxonomy labels associated with the data
  data_labels_1 <- setNames(unlist(get_labels(data_1)), NULL)
  data_labels_2 <- setNames(unlist(get_labels(data_2)), NULL)

  # Create membership tables
  membership_1 <- data.frame('node' = 1:length(c(tree_1$tip.label, tree_1$node.label)),
                             'membership' = c(tree_1$tip.label, tree_1$node.label) %in% data_labels_2)
  membership_2 <- data.frame('node' = 1:length(c(tree_2$tip.label, tree_2$node.label)),
                             'membership' = c(tree_2$tip.label, tree_2$node.label) %in% data_labels_1)

  # Create tree plots
  plot_1 <- ggtree(tree_1) %<+% membership_1 +
    aes(color = membership) +
    layout_circular() +
    scale_color_manual(name = paste('Labels from', name_2, 'that are in', name_1),
                       labels = c('True', 'False'),
                       values = c('TRUE' = '#66c2a5', 'FALSE' = '#cccccc')) +
    ggtitle(paste0(name_1, ' subtree')) +
    theme(plot.title = element_text(hjust = 0.5))

  plot_2 <- ggtree(tree_2) %<+% membership_2 +
    aes(color = membership) +
    layout_circular() +
    scale_color_manual(name = paste('Labels from', name_1, 'that are in', name_2),
                       labels = c('True', 'False'),
                       values = c('TRUE' = '#66c2a5', 'FALSE' = '#cccccc')) +
    ggtitle(paste0(name_2, ' subtree'))  +
    theme(plot.title = element_text(hjust = 0.5))

  if (show_tips){
    # Adjust the tip label size depending on the number of tips of the tree
    if (length(tree_1$tip.label) <= 200){
      tip_size_1 = 3
    } else if (length(tree_1$tip.label) <= 500){
      tip_size_1 = 1.5
    } else {
      tip_size_1 = .5
    }

    # Adjust the tip label size depending on the number of tips of the tree
    if (length(tree_2$tip.label) <= 200){
      tip_size_2 = 3
    } else if (length(tree_1$tip.label) <= 500){
      tip_size_2 = 1.5
    } else {
      tip_size_2 = .5
    }

    plot_1 <- plot_1 + geom_tiplab(size = tip_size_1)
    plot_1 <- plot_1 + xlim(0, max(plot_1$data$x) + 50)

    plot_2 <- plot_2 + geom_tiplab(size = tip_size_2)
    plot_2 <- plot_2 + xlim(0, max(plot_2$data$x) + 50)
  }

  plot_1 <- plot_1 + theme(legend.justification = 'top')
  plot_2 <- plot_2 + theme(legend.justification = 'top')

  return(list(plot_1, plot_2))
}



