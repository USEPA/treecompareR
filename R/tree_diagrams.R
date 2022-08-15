



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
#' This is a helper function for retrieving labels in a data.table of classified
#' chemicals, grouped by taxonomy level.
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

#' Number of labels
#'
#' This is a helper function to determine number of occurrences for each label
#' in a data.table containing chemicals and their classification data.
#'
#' @param data A data.table of chemicals with classifications.
#' @param tax_level_labels An alternate parameter giving the taxonomy levels if
#'   not  using ClassyFire taxonomy
#' @return A named list of chemical labels and their number of occurrences.
get_number_of_labels <- function(data, tax_level_labels = NULL){
  if (is.null(tax_level_labels)){
    tax_level_labels <- c('kingdom', 'superclass', 'class', 'subclass',
                          'level5', 'level6', 'level7', 'level8',
                          'level9', 'level10', 'level11')
  }

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


#' Display subtree
#'
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
  cohort <- NULL
  dataset_1 <- NULL

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

  analytes_data <- data.frame('node' = 1:length(tree_labels))
  analytes_data[['dataset_1']] <- data_1_subtree

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

    analytes_data[["dataset_2"]] <- data_2_subtree

    all_cohort <- ifelse(data_1_subtree, 0L, 2L) + ifelse(data_2_subtree, 0L, 1L)

    analytes_data[["cohort"]] <- as.character(all_cohort)

    #print(summary(analytes_data))


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
prune_and_display_subtree <- function(data, tax_level_labels = NULL, tree = NULL, show_tips = TRUE, no_plot = FALSE, adjust_branch_length = TRUE, xlimit = c(0, 150)) {
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
  #pruned_tree <- ape::drop.tip(tree,
  #                             setdiff(tree$tip.label, intersect(data_labels, tree$tip.label)))
  pruned_tree <- drop_tips_nodes(tree = tree, data = data, tax_level_labels = tax_level_labels)
  if (adjust_branch_length) {
    pruned_tree$edge.length <- adjust_branch_lengths(pruned_tree)
  }

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
    layout_circular() + xlim(xlimit[[1]], xlimit[[2]])

  # If displaying the tip labels, add them in with the correct size
  if (show_tips)
    tree_plot <- tree_plot + geom_tiplab(size = tip_size)
  return(list(pruned_tree, tree_plot))
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
data_set_subtrees <- function(data_1, data_2, name_1 = 'data_1', name_2 = 'data_2', tax_level_labels = NULL, tree = NULL, show_tips = TRUE, adjust_branch_length = TRUE){
  membership <- NULL
  # Prune subtrees for each data set.
  tree_1 <- prune_and_display_subtree(data_1, tax_level_labels = tax_level_labels, tree = tree, show_tips = show_tips, no_plot = TRUE, adjust_branch_length = adjust_branch_length)
  tree_2 <- prune_and_display_subtree(data_2, tax_level_labels = tax_level_labels, tree = tree, show_tips = show_tips, no_plot = TRUE, adjust_branch_length = adjust_branch_length)

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
#' @return An `aplot` consisting of two outer layer `ggtree` objects and three
#'   inner layer `ggplot2` objects.
#' @export
side_by_side_trees <- function(data_left, data_right, name_left = 'Left tree', name_right = 'Right tree'){
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

  print(all_left_tree)
  print(all_right_tree)

  left_tree <- ggtree(union_tree,
                      aes(color= (c(union_tree$tip.label, union_tree$node.label) %in% all_left_tree)),
                      branch.length = FALSE)+
    scale_color_manual(values = c('black', 'blue'),
                       labels = c('FALSE', name_left),
                       name = 'Left Tree') + geom_nodepoint() + geom_tippoint()# + geom_tiplab(size = 3)
  right_tree <- ggtree(union_tree,
                       aes(color= c(union_tree$tip.label, union_tree$node.label) %in% all_right_tree),
                       branch.length = FALSE) +
    scale_color_manual(values = c('black', 'red'),
                       labels = c('FALSE', name_right),
                       name = 'Right Tree') + geom_nodepoint()#geom_tiplab(size = 3)
  right_tree <- right_tree + ggplot2::scale_x_continuous(trans = "reverse")

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

  print(tree_data_leftval)
  print(tree_data_rightval)
  print(tree_data_centerval)
  print(nTip)

  value <- c(tree_data_leftval, tree_data_rightval, tree_data_centerval)
  tree_data <- cbind(tree_data, value)
  names(tree_data)[[3]] <- "value"

  print(names(tree_data))

  print(str(tree_data))
  #left_data_plot <- ggplot(tree_data, aes(y = tip.label, x = leftval)) +
  #  geom_tile()#Or some other data viz
  #return(left_data_plot)
  #center_data_plot <- ggplot(tree_data, aes(y = tip.label, x = centerval)) +
  #  geom_tile()#Or some other data viz
  #right_data_plot <- ggplot(tree_data, aes(y = tip.label, x = rightval)) +
  #  geom_tile()#Or some other data viz

  #return(tree_data)

  data_plot <- ggplot(tree_data, aes(x = tree, y = tip.label)) +
    geom_tile(aes(fill = value)) + theme_minimal() + ylab(NULL)# + theme(axis.text.y = element_text(size = 5))

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

