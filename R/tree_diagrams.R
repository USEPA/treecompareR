
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
                              tax_level_labels = chemont_tax_levels){
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
#'  variable passed to \code{base_tree} will be used as the title (with "tree"
#'  appended).
#'@param data_1 Highlight branches of the base tree according to their
#'  membership in this list. One of the following options: A data.frame
#'  consisting of a list of chemicals with classification data; a list of node
#'  numbers in the base tree; a list of node labels in the base tree; or another
#'  phylo object.
#'@param data_2 Optional: highlight branches of the base tree to compare
#'  membership in this list and \code{data_1}. One of the following options: A
#'  data.frame consisting of a list of chemicals with classification data; a
#'  list of node numbers in the base tree; a list of node labels in the base
#'  tree; or another phylo object.
#'@param name_1 Optional: A string giving the name of the list in \code{data_1},
#'  for plot labeling. Default is the variable name passed to \code{data_1}.
#'@param name_2 Optional: A string giving the name of the list in \code{data_2}
#'  (if any), for plot labeling. Default is the variable name passed to \code{data_2}.
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
                        subtree_mapping_default[setdiff(names(subtree_mapping_default),
                                                        names(subtree_mapping))])


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

    if(is.data.frame(data_1)){
  #Get node numbers of data_1 subtree
  data_1_all <- get_subtree_nodes(data = data_1,
                                  base_tree = base_tree,
                                  tax_level_labels = tax_level_labels)
    }else if(is.numeric(data_1)){
      #interpret as node numbers
      #get ancestors of these nodes
    data_1_all <- unique(c(data_1,
                    unlist(phangorn::Ancestors(x = base_tree,
                                               node = data_1)))
    )
    }else if(is.character(data_1)){
    #interpret as node labels
      data_1_nodes <- get_node_from_label(label = data_1,
                                          tree = base_tree)
      data_1_all <- unique(c(data_1_nodes,
                             unlist(phangorn::Ancestors(x = base_tree,
                                                        node = data_1_nodes))))
    }else if("phylo" %in% class(data_1)){
      #take node and tip labels
      data_1_labels <- c(data_1$tip.label,
                         data_1$node.label)
      #get nodes
      data_1_nodes <- get_node_from_label(label = data_1_labels,
                                          tree = base_tree)
      #get ancestors of these nodes
      data_1_all <- unique(c(data_1_nodes,
                             unlist(phangorn::Ancestors(x = base_tree,
                                                        node = data_1_nodes))))
    }

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

    if(is.data.frame(data_2)){
      #Get node numbers of data_2 subtree
      data_2_all <- get_subtree_nodes(data = data_2,
                                      base_tree = base_tree,
                                      tax_level_labels = tax_level_labels)
    }else if(is.numeric(data_2)){
      #interpret as node numbers
      #get ancestors of these nodes
      data_2_all <- unique(c(data_2,
                             unlist(phangorn::Ancestors(x = base_tree,
                                                        node = data_2)))
      )
    }else if(is.character(data_2)){
      #interpret as node labels
      data_2_nodes <- get_node_from_label(label = data_2,
                                          tree = base_tree)
      data_2_all <- unique(c(data_2_nodes,
                             unlist(phangorn::Ancestors(x = base_tree,
                                                        node = data_2_nodes))))
    }else if("phylo" %in% class(data_2)){
      #take node and tip labels
      data_2_labels <- c(data_2$tip.label,
                         data_2$node.label)
      #get nodes
      data_2_nodes <- get_node_from_label(label = data_2,
                                          tree = base_tree)
      #get ancestors of these nodes
      data_2_all <- unique(c(data_2_nodes,
                             unlist(phangorn::Ancestors(x = base_tree,
                                                        node = data_2_nodes))))
    }

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
  tree_plot <- add_cladelab(tree_plot = tree_plot,
                                          tree = base_tree,
                                          clade_level = clade_level,
                                          clade_opts = clade_opts)

  #add title of base tree
  tree_plot <- tree_plot + ggtitle(base_name)

  return(tree_plot)

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
                                      tax_level_labels = chemont_tax_levels,
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
#'   by terminal (tip) labels, for calculation and plotting of number of chemicals and overlap.
#' @param tax_level_labels
#' @param annot_angle
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
