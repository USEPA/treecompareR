



#' This function returns the terminal classification label for classified chemicals
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
add_terminal_label <- function(data, tip = FALSE, labels = NULL, tree = NULL){
  if (is.null(labels)){
    labels <- c('kingdom', 'superclass', 'class', 'subclass',
                         'level5', 'level6', 'level7', 'level8',
                         'level9', 'level10', 'level11')
  }

  # Check to see if the table already has a terminal_label column
  if('terminal_label' %in% names(data))
    return(data)
  indices <- which(names(data) %in% labels)
  if (length(indices) != length(labels))
    stop('Missing col(s) \n', paste(base::setdiff(labels, names(data)[indices], collapse = '\n', '!')))


  new_dt <- copy(data)
  new_dt[, c("terminal_label") := terminal_function(t = .SD, tip = tip, tax_level_labels = labels, tree = tree), by = 1:nrow(new_dt)]

  return(new_dt)
}

#' Helper function to return the terminal label for `add_terminal_label()`
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

#' Returns unique labels in a classified data for a given taxonomy level
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
    stop('Please input a valid label!')

  # Collect the labels
  labels <- data[, sapply(.SD, unique), .SDcols = c(level_label)]

  # Remove NA from list of labels
  labels <- labels[!is.na(labels)]

  # Remove '' from list of labels
  labels <- labels[which(labels != '')]

  return(labels)
}


#' Helper function for retrieving labels in a data.table of classified chemicals
#'
#' @param data A data.table consisting of classified chemicals.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @return A list of classification labels for each level of taxonomy.
get_labels <- function(data, tax_level_labels = NULL){
  if (is.null(tax_level_labels)){
    tax_level_labels <- c('kingdom', 'superclass', 'class', 'subclass',
                          'level5', 'level6', 'level7', 'level8',
                          'level9', 'level10', 'level11')
  }
  labels <- sapply(tax_level_labels, function(t) {get_label_level(data, t, tax_level_labels)})
  labels
}

#' A helper function that returns the number of labels per level.
#'
#' @param label_list A named list of labels corresponding to taxonomy levels.
#' @return The number of labels per taxonomy level.
get_label_length <- function(label_list){
  lengths <- sapply(label_list, length)
  lengths
}


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
#' @import tidyr
#' @import ggplot2
label_bars <- function(data = NULL, tax_level_labels = NULL){
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


#' This function takes in one or two data.tables consisting of chemicals with
#' classification data and returns a tree diagram indiacting the subtree(s)
#' induced by the data.table(s).
#'
#' @param data_1 A data.table consisting of classification data for chemicals.
#' @param data_2 An alternate parameter giving a second data.table with
#'   classification data for chemicals.
#' @param name_1 An alternate parameter giving the name of `data.table_1`.
#' @param name_2 An alternate parameter giving the name of `data.table_1`.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @return A ggtree object showing the subtree induced by `data_1` (and `data_2`
#'   if supplied, and their intersection.).
#' @export
#' @import phangorn
#' @import ggtree
display_subtree <- function(data_1, data_2 = NULL, name_1 = NULL, name_2 = NULL, tree = NULL, tax_level_labels = NULL){
  if (is.null(tax_level_labels)){
    tax_level_labels <- c('kingdom', 'superclass', 'class', 'subclass',
                          'level5', 'level6', 'level7', 'level8',
                          'level9', 'level10', 'level11')
  }

  if (is.null(tree)){
    tree <- chemont_tree
  }

  tree_labels <- c(tree$tip.label, tree$node.label)


  # Get all labels associated with the input data sets
  data_1_labels <- setNames(unlist(get_labels(data_1, tax_level_labels)), NULL)
  select_data_1 <- sapply(data_1_labels, function(t) {match(t, tree_labels)})
  branches_data_1 <- sapply(select_data_1, phangorn::Ancestors, x = tree)
  data_1_all <- unique(unlist(branches_data_1))

  # Create data.frame to store subtree information
  data_1_subtree <- rep(FALSE, length(tree_labels))
  data_1_subtree[select_data_1] <- TRUE
  data_1_subtree[data_1_all] <- TRUE

  analytes_data <- data.frame(node = 1:length(tree_labels),
                              dataset_1 = data_1_subtree)

  #print(summary(analytes_data))

  if (is.null(name_1)){
    name_1 <- 'Set 1'
  }


  if (!is.null(data_2)) {
    data_2_labels <- setNames(unlist(get_labels(data_2, tax_level_labels)), NULL)
    select_data_2 <- sapply(data_2_labels, function(t) {match(t, tree_labels)})
    branches_data_2 <- sapply(select_data_2, phangorn::Ancestors, x = tree)
    data_2_all <- unique(unlist(branches_data_2))

    data_2_subtree <- rep(FALSE, length(tree_labels))
    data_2_subtree[select_data_2] <- TRUE
    data_2_subtree[data_2_all] <- TRUE

    analytes_data["dataset_2"] <- data_2_subtree

    all_cohort <- ifelse(data_1_subtree, 0L, 2L) + ifelse(data_2_subtree, 0L, 1L)

    analytes_data["cohort"] <- as.character(all_cohort)

    print(summary(analytes_data))


    if (is.null(name_2)){
      name_2 <- 'Set 2'
    }

    tree_plot <- ggtree(tree) %<+% analytes_data +
      aes(color = cohort) +
      layout_circular() +
      geom_tiplab(size = .5) +
      scale_color_manual(name = 'Cohort',
                         labels = c(paste0(name_1, ' and ', name_2), paste0('Only ', name_1),
                                    paste0('Only ', name_2), 'Neither set'),
                         values = c('0' = '#800026', '1' = '#ffeda0', '2' = '#41ab5d', '3' = '#f0f0f0'))

    return(tree_plot)
  }

  tree_plot <- ggtree(tree) %<+% analytes_data +
    aes(color = dataset_1) +
    layout_circular() +
    geom_tiplab(size = 0.5) +
    scale_color_manual(name = paste0(name_1, ' subtree'),
                       labels = c('TRUE', 'FALSE'),
                       values = c('TRUE' = '#800026', 'FALSE' = '#f0f0f0'))+
    guides(color = guide_legend(ovveride.aes = list(size = 2)))

  return(tree_plot)

}

#' This function takes a data.table of classified chemicals and returns a tree
#' with optional tree diagram
#'
#' @param data A data.table with classification data for chemicals.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @param show_tips An alternate parameter determining whether tip labels are
#'   displayed.
#' @param no_plot An alternate parameter for returning only the pruned tree
#'   without the tree visual.
#' @return A pruned tree or list consisting of a pruned tree and ggtree diagram
#'   of the pruned tree.
#' @export
#' @import ape
#' @import ggtree
prune_and_display_subtree <- function(data, tax_level_labels = NULL, tree = NULL, show_tips = TRUE, no_plot = FALSE) {
  if (is.null(tax_level_labels)){
    tax_level_labels <- c('kingdom', 'superclass', 'class', 'subclass',
                          'level5', 'level6', 'level7', 'level8',
                          'level9', 'level10', 'level11')
  }

  if (is.null(tree)){
    tree <- chemont_tree
  }

  # Get all taxonomy labels associated with the data
  data_labels <- setNames(unlist(get_labels(data, tax_level_labels)), NULL)

  # Prune the labels not represented by the data from the full tree
  pruned_tree <- ape::drop.tip(tree,
                               setdiff(tree$tip.label, intersect(data_labels, tree$tip.label)))
  if (no_plot)
    return(pruned_tree)

  # Adjust the tip label size depending on the number of tips of the tree
  if (length(pruned_tree$tip.label) <= 200){
    tip_size = 3
  } else if (length(pruned_tree$tip.label) <= 500){
    tip_size = 1.5
  } else {
    tip_size = .5
  }

  # Create the tree plot
  tree_plot <- ggtree(pruned_tree) +
    layout_circular()

  # If displaying the tip labels, add them in with the correct size
  if (show_tips)
    tree_plot <- tree_plot + geom_tiplab(size = tip_size)
  return(list(pruned_tree, tree_plot))
}

# This function takes in a data.table with column `terminal_label`, and a column
# name and produces a circular plot with boxplots displaying values from the
# specified column grouped by each value in `terminal_label`

#' This function takes in a data.table of chemicals with classification data and
#' additional numeric data and displays selected numeric data grouped by tip
#' label on data-induced subtree.
#'
#' @param data A data.table consisting of classification data and additional
#'   numeric data.
#' @param col A specified numeric column for use in displaying boxplots.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @return A ggtree object consisting of subtree induced by data and boxplots
#'   corresponding to specified numeric column of data.
#' @export
#' @import ggtree
#' @import ggtreeExtra
circ_tree_boxplot <- function(data, col, tax_level_labels = NULL, tree = NULL){
  if (!(col %in% names(data)))
    stop(paste('The column', col, 'is not in the input data!'))

  #print(col)
  columns <- which(names(data) %in% c(col, 'terminal_label'))
  #print(columns)

  new_data <- copy(as.data.frame(data[, .SD, .SDcols = columns]))
  index <- which(names(new_data) %in% col)
  new_data[, index] <- as.numeric(new_data[, index])
  new_data <- new_data[!is.na(new_data[, index]),]
  new_data["grp"] <- new_data$terminal_label
  new_data["val"] <- new_data[, index]
  new_data["node"] <- new_data$terminal_label

  #summary(new_data)

  new_data_tree <- prune_and_display_subtree(data, tax_level_labels, tree, show_tips = FALSE, no_plot = TRUE)
  #print(new_data_tree)
  num_nodes_tips <- length(new_data_tree$tip.label) + new_data_tree$Nnode
  #print(num_nodes_tips)

  circ_plot <- ggtree(new_data_tree,
                      layout = 'circular')
  circ_plot <- circ_plot + geom_tippoint()

  circ_plot <- circ_plot + ggtreeExtra::geom_fruit(data = new_data, geom = geom_boxplot,
                                      mapping = aes(x = val,
                                                    y = terminal_label,
                                                    fill = grp),
                                      size = 0.2,
                                      outlier.size = 0.5,
                                      outlier.stroke = 0.08,
                                      outlier.shape = 21,
                                      axis.params = list(axis = 'x',
                                                         text.size = 1.8),
                                      grid.params = list())
  circ_plot <- circ_plot + scale_fill_discrete(name = 'Tip label',
                                               guide = guide_legend(keywidth = 0.2,
                                                                    keyheight = 0.2,
                                                                    ncol = 2))
  return(circ_plot)
}


# This function takes in two data sets, creates the subtree induced by data_1,
# and colors the tips based on the number of chemicals from data_2 are in the
# set of chemicals from data_1 grouped by tip.

#' This function takes in two data.tables, plots the subtree induced by the
#' first and colors the tips based on the proportion of chemicalsfrom the second
#' that make up the chemials from the first, grouped by tip label.
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
leaf_fraction_subtree <- function(data_1, data_2, name_1 = 'data_1', name_2 = 'data_2', show_labels = FALSE, tax_level_labels = NULL, tree = NULL){
  # Find all the terminal_label values from data_1.
  terminal_labels <- data_1[!is.na(terminal_label), unique(terminal_label)]

  # For each terminal_label value, determine the chemicals from data_2 that are
  # also in data_1. This checks using the INCHIKEY of each chemical.
  label_percentages <- sapply(terminal_labels, function(t) {
    data_1_chemicals <- data_1[terminal_label == t, unique(INCHIKEY)]
    data_2_chemicals <- data_2[terminal_label == t, unique(INCHIKEY)]
    shared_chemicals <- intersect(data_1_chemicals, data_2_chemicals)
    return(length(shared_chemicals)/length(data_1_chemicals))
  })

  label_data <- data.frame(label = terminal_labels,
                           percentages = unname(label_percentages),
                           data_1_numbers <- unname(sapply(terminal_labels, function(t) {
                             length(data_1[terminal_label == t, unique(INCHIKEY)])
                           })),
                           data_2_numbers <- unname(sapply(terminal_labels, function(t) {
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
  tree_plot <- tree_plot + ggtitle(paste0(name_1, ' subtree'))

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


#' This function creates two plots that show a visual comparison of the subtrees induced by two data sets
#'
#' @param data_1 A data.table of chemicals and their classifications.
#' @param data_2 A data.table of chemicals and their classifications.
#' @param name_1 An alternate parameter giving the name of the first data set.
#' @param name_2 An alternate parameter giving the name of the second data set.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not using ClassyFire taxonomy.
#' @param tree An alternate parameter giving a taxonomy if not using ChemOnt.
#' @param show_tips An alternate parameter determining whether to show tip labels.
#' @return A list of two ggtree objects.
#' @export
#' @import ggtree
data_set_subtrees <- function(data_1, data_2, name_1 = 'data_1', name_2 = 'data_2', tax_level_labels = NULL, tree = NULL, show_tips = TRUE){
  # Prune subtrees for each data set.
  tree_1 <- prune_and_display_subtree(data_1, tax_level_labels = tax_level_labels, tree = tree, show_tips = show_tips, no_plot = TRUE)
  tree_2 <- prune_and_display_subtree(data_2, tax_level_labels = tax_level_labels, tree = tree, show_tips = show_tips, no_plot = TRUE)

  # Get all taxonomy labels associated with the data
  data_labels_1 <- setNames(unlist(get_labels(data_1)), NULL)
  data_labels_2 <- setNames(unlist(get_labels(data_2)), NULL)

  # Create membership tables
  membership_1 <- data.frame(node = 1:length(c(tree_1$tip.label, tree_1$node.label)),
                             membership = c(tree_1$tip.label, tree_1$node.label) %in% data_labels_2)
  membership_2 <- data.frame(node = 1:length(c(tree_2$tip.label, tree_2$node.label)),
                             membership = c(tree_2$tip.label, tree_2$node.label) %in% data_labels_1)

  # Create tree plots
  plot_1 <- ggtree(tree_1) %<+% membership_1 +
    aes(color = membership) +
    layout_circular() +
    scale_color_manual(name = paste('Labels from', name_2, 'that are in', name_1),
                       labels = c('True', 'False'),
                       values = c('TRUE' = '#66c2a5', 'FALSE' = '#cccccc'))

  plot_2 <- ggtree(tree_2) %<+% membership_2 +
    aes(color = membership) +
    layout_circular() +
    scale_color_manual(name = paste('Labels from', name_1, 'that are in', name_2),
                       labels = c('True', 'False'),
                       values = c('TRUE' = '#66c2a5', 'FALSE' = '#cccccc'))

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
    plot_2 <- plot_2 + geom_tiplab(size = tip_size_2)
  }


  return(list(plot_1, plot_2))
}



