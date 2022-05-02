

#' Helper function for determining number of labels in a given data set
#'
#' @param datatable A data.table object of chemical classifications
#' @return Number of occurrences of each label in the parameter 'datatable'.
#' @import data.table

label_numbers <- function(datatable, chemont = TRUE, log = TRUE) {
  if (!chemont){
    top_level <- names(datatable)[[1]]
    complete_labels <- 'TO BE DETERMINED'
    return('TBD')
  }

  complete_labels <- c(unlist(unname(datatable[kingdom != '', .(unique(kingdom),
                                                      unique(superclass),
                                                      unique(class),
                                                      unique(subclass),
                                                      unique(level5),
                                                      unique(level6),
                                                      unique(level7),
                                                      unique(level8),
                                                      unique(level9),
                                                      unique(level10),
                                                      unique(level11)),
                                        by = .(PREFERRED_NAME)][, PREFERRED_NAME := NULL])))
  print(length(complete_labels))
  empty_indices <- which(sapply(complete_labels, function(t) {t == ''}))
  if (length(empty_indices) > 0){
    complete_labels <- complete_labels[-which(sapply(complete_labels, function(t) {t == ''}))]
  }

  complete_labels <- complete_labels[-which(sapply(complete_labels, is.na))]

  print(length(complete_labels))
  unique_labels <- unique(complete_labels)

  if (log) {
    unique_complete_label_numbers <- sapply(unique_labels, function(t){
      log10(length(which(complete_labels %in% t)))
      })
    print(length(unique_complete_label_numbers))
  } else {
    unique_complete_label_numbers <- sapply(unique_labels, function(t){
      length(which(complete_labels %in% t))
      })
  }
  return (list(unique_complete_label_numbers, complete_labels))
}








#' Generate a heatmap from input matrix and indices.
#'
#' @param tree_object A phylo object representing a rooted tree.
#' @param matrix A matrix of similarity measure values derived from parameter 'tree_object'.
#' @param row_indices The row indices for the matrix.
#' @param column_indices The column indices for the matrix.
#' @param row_data A data.table object of chemical classifications.
#' @param column_data A data.table object of chemical classifications.
#' @param name Name of the heatmap
#' @param row_split Number of clusters for rows.
#' @param col_split Number of cluster for columns.
#' @param row_title Title for rows.
#' @param col_title Title for columns.
#' @return A Heatmap object.
#' @export
#' @import ComplexHeatmap
#' @import circlize
#' @import viridis
#' @import grid


generate_heatmap <- function(tree_object, matrix, row_indices, column_indices, row_data, column_data, name = 'Name', row_split = 1, column_split = 1, row_title = 'Row title', column_title = 'Column title') {
  if (!identical(unlist(names(row_data)), unlist(names(column_data))))
      stop('The classification levels for the row data and column data do not match!')
  taxonomy_names <- names(row_data)
  # COLLECT LABEL NUMBERS FOR ROW DATA AND FOR COLUMN DATA
  row_label_data <- label_numbers(row_data)
  row_label_numbers <- row_label_data[[1]]
  row_labels <- row_label_data[[2]]
  row_anno_indices <- match(dimnames(matrix)[[2]][row_indices], row_labels)
  row_na_indices <- which(sapply(row_anno_indices, is.na))
  if (length(row_na_indices) > 0){
    row_anno_indices <- row_anno_indices[-row_na_indices]
  }

  column_label_data <- label_numbers(column_data)
  column_label_numbers <- column_label_data[[1]]
  column_labels <- column_label_data[[2]]
  column_anno_indices <- match(dimnames(matrix)[[2]][column_indices], column_labels)
  column_na_indices <- which(sapply(column_anno_indices, is.na))
  if (length(column_na_indices) > 0){
    column_anno_indices <- column_anno_indices[-column_na_indices]
  }


  matrix_row_indices <- intersect(which(dimnames(matrix)[[2]] %in% row_labels), row_indices)
  matrix_column_indices <- intersect(which(dimnames(matrix)[[2]] %in% column_labels), column_indices)



  heatmap <- ComplexHeatmap::Heatmap(#matrix = matrix[row_indices, column_indices],
                                     matrix = matrix[matrix_row_indices, matrix_column_indices],
                                     name = name,
                                     col = circlize::colorRamp2(seq(0, 1, len = 20), viridis::viridis(20, option = 'C')),

                                     # NEED TO ADD HELPER FUNCTIONS FOR THIS
                                     top_annotation = HeatmapAnnotation(#col_log_count_bar = anno_barplot(column_label_numbers[match(dimnames(matrix)[[2]][column_indices], column_labels)]),
                                                                        col_log_count_bar = anno_barplot(column_label_numbers[column_anno_indices]),
                                                                        annotation_name_rot = 45,
                                                                        annotation_label = c('log(col count) bars'),
                                                                        annotation_name_gp = grid::gpar(fontsize = 8)
                                     ),
                                     left_annotation = rowAnnotation(#row_log_count_bar = anno_barplot(row_label_numbers[match(dimnames(matrix)[[1]][row_indices], row_labels)],
                                                                      row_log_count_bar = anno_barplot(row_label_numbers[row_anno_indices],
                                                                      axis_param = list(direction = 'reverse')),
                                                                      annotation_name_rot = 45,
                                                                      annotation_label = c('log(row count) bars'),
                                                                      annotation_name_gp = grid::gpar(fontsize = 8)
                                                                      ),
                                     show_row_names = FALSE,
                                     show_column_names = FALSE,
                                     row_split = row_split,
                                     column_split = column_split
                                     )
  heatmap <- draw(heatmap,
                  row_title = row_title,
                  row_title_gp = grid::gpar(fontsize = 10, fontface = 'bold'),
                  column_title = column_title,
                  column_title_gp = grid::gpar(fontsize = 10, fontface = 'bold'))
}





















#############################################################
#NEED TO SPECIFY get_tip_level() FOR THE FUNCTION BELOW!!!!!#
#############################################################

#' Helper function that returns superclasses or classes for specified clusters in heatmap
#'
#' @param htmap A ComplexHeatmap object with hierarchical clustering.
#' @param row_cluster Index for the row cluster.
#' @param column_cluster Index for the column cluster.
#' @param level Alternate parameter indicating the level of depth for labels.
#' @param tree_object A phylo object representing a root tree, the taxonomy being investigated.
#' @param tree Alternate parameter, a phylo object representing a rooted tree, for restricting the labels.
#' @return A list of labels for the row and column clusters, based on the level specified.
#' @import stats
#' @import ComplexHeatmap
cluster_analysis <- function(htmap, row_cluster, column_cluster, level = 2, tree_object, tree = NULL){
  if (!is.null(tree)){
    tree_labels <- c(tree$tip.label, tree$node.label)
  }
  # get row levels for row cluster (restrict to tree if tree is given)
  row_names <- dimnames(htmap@ht_list[[1]]@matrix)[[1]][stats::order.dendrogram(ComplexHeatmap::row_dend(htmap)[[row_cluster]])]
  if (is.null(tree)){
    row_levels <- sort(unique(sapply(row_names, get_tip_level, tree = tree_object)))
  } else {
    row_levels <- sort(unique(sapply(intersect(tree_labels, row_names), get_tip_level, tree = tree_object)))
  }

  # get column levels for column cluster
  column_names <- dimnames(htmap@ht_list[[1]]@matrix)[[2]][stats::order.dendrogram(ComplexHeatmap::column_dend(htmap)[[column_cluster]])]
  if (is.null(tree)){
    column_levels <- sort(unique(sapply(column_names, get_tip_level, tree = tree_object)))
  } else {
    column_levels <- sort(unique(sapply(intersect(tree_labels, column_names), get_tip_level, tree = tree_object)))
  }

  if (level == 2){
    # get superclass labels per level for row cluster
    if (is.null(tree)) {
      row_superclass <- sapply(row_levels, function(q) {
        unique(unlist((sapply(row_names[which(sapply(row_names, function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-2]
        }
        )
        )
        )
        )
      }
      )
    } else {
      row_superclass <- sapply(row_levels, function(q) {
        unique(unlist((sapply(intersect(tree_labels, row_names)[which(sapply(intersect(tree_labels, row_names), function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-2]
        }
        )
        )
        )
        )
      }
      )
    }

    # get superclass labels per level for column cluster
    if (is.null(tree)) {
      column_superclass <- sapply(column_levels, function(q) {
        unique(unlist((sapply(column_names[which(sapply(column_names, function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-2]
        }
        )
        )
        )
        )
      }
      )
    } else {
      column_superclass <- sapply(column_levels, function(q) {
        unique(unlist((sapply(intersect(tree_labels, column_names)[which(sapply(intersect(tree_labels, column_names), function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-2]
        }
        )
        )
        )
        )
      }
      )
    }
    return(list('row_superclass' = unique(unlist(row_superclass)),
                'column_superclass' = unique(unlist(column_superclass))))
  } else if (level == 3) {
    # get class labels per level for row cluster
    if (is.null(tree)){
      row_class <- sapply(row_levels, function(q) {
        unique(unlist((sapply(row_names[which(sapply(row_names, function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-3]
        }
        )
        )
        )
        )
      }
      )
    } else {
      row_class <- sapply(row_levels, function(q) {
        unique(unlist((sapply(intersect(tree_labels, row_names)[which(sapply(intersect(tree_labels, row_names), function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-3]
        }
        )
        )
        )
        )
      }
      )
    }

    # get class labels per level for column cluster
    if (is.null(tree)){
      column_class <- sapply(column_levels, function(q) {
        unique(unlist((sapply(column_names[which(sapply(column_names, function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-3]
        }
        )
        )
        )
        )
      }
      )
    } else {
      column_class <- sapply(column_levels, function(q) {
        unique(unlist((sapply(intersect(tree_labels, column_names)[which(sapply(intersect(tree_labels, column_names), function(t) {
          unname(get_tip_level(tree = tree_object, t))
        }
        ) == q)], function(s) {
          get_ancestors(tree_object, s)[q-3]
        }
        )
        )
        )
        )
      }
      )
    }
    return(list('row_class' = unique(unlist(row_class)),
                'column_class' = unique(unlist(column_class))))
  }

}



#############################################################
#NEED TO SPECIFY get_tip_level() FOR THE FUNCTION BELOW!!!!!#
#############################################################


#' Helper function to display clade label associated with missing node.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param tree_object A phylo object representing the entire taxonomy from which the parameter 'tree' is derived.
#' @param list_superclasses A list of superclasses to be labeled.
#' @param tree_visual A ggtree object that will have clades labeled.
#' @param i The index for which list_superclasses element will be labeled.
#' @param color A color for the clade label.
#' @return A ggtree object that will have the specified clade labeled.
#' @import phangorn
#' @import ggplot2
#' @import ggtree
handle_missing_node_show_clade <- function(tree, tree_object, list_superclasses, tree_visual, i, color){
  #print(paste('superclasses:', list_superclasses))
  tree_labels <- c(tree$tip.label, tree$node.label)
  #print(intersect(c(tree_object$tip.label, tree_object$node.label)[Descendants(tree_object, which(c(tree_object$tip.label, tree_object$node.label) %in% list_superclasses[[i]]), type = 'all')], tree_labels))

  temp_descendants <- intersect(c(tree_object$tip.label, tree_object$node.label)[phangorn::Descendants(tree_object, which(c(tree_object$tip.label, tree_object$node.label) %in% list_superclasses[[i]]), type = 'all')], tree_labels)
  #print(temp_descendants)
  shallow_level <- min(sapply(temp_descendants, get_tip_level, tree = tree_object))
  shallow_descendants <- temp_descendants[which(sapply(temp_descendants, get_tip_level, tree = tree_object) == shallow_level)]
  middle <- (length(shallow_descendants)+1)%/% 2
  #print(middle)
  for (j in seq_along(shallow_descendants)){
    #print(j)
    #print(middle == j)
    tree_visual <- tree_visual + ggtree::geom_cladelab(node = which(tree_labels %in% shallow_descendants[[j]]),
                                               label = ifelse(j == middle, list_superclasses[[i]],''),
                                               textcolor = color,
                                               barcolor = color,
                                               offset = .2,
                                               offset.text = 3,
                                               fontsize = 2.6,
                                               angle = 'auto')
  }
  return(tree_visual)

}

#############################################################
#NEED TO SPECIFY get_tip_level() FOR THE FUNCTION BELOW!!!!!#
#############################################################

#' Helper function to display clade highlight associated with missing node.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param tree_object A phylo object representing the entire taxonomy from which parameter 'tree' is derived.
#' @param list_superclasses A list of superclasses to be highlighted.
#' @param tree_visual A ggtree object that will have clades highlighted.
#' @param i The index for which list_superclasses element will be highlighted.
#' @param color A color for the clade highlight.
#' @return A ggtree object that will have the specified clade highlighted.
#' @import phangorn
#' @import ggplot2
#' @import ggtree
handle_missing_node_highlight_clade <- function(tree, tree_object, list_superclasses, tree_visual, i, color){
  tree_labels <- c(tree$tip.label, tree$node.label)
  temp_descendants <- intersect(c(tree_object$tip.label, tree_object$node.label)[phangorn::Descendants(tree_object, which(c(tree_object$tip.label, tree_object$node.label) %in% list_superclasses[[i]]), type = 'all')], tree_labels)
  #print(temp_descendants)
  shallow_level <- min(sapply(temp_descendants, get_tip_level, tree = tree_object))
  shallow_descendants <- temp_descendants[which(sapply(temp_descendants, get_tip_level, tree = tree_object) == shallow_level)]
  middle <- (length(shallow_descendants)+1)%/% 2
  if (length(shallow_descendants) == 0) {
    print('whomp')
    return(tree_visual)
  }
  for (j in seq_along(shallow_descendants)){
    tree_visual <- tree_visual + ggtree::geom_hilight(node = which(tree_labels %in% shallow_descendants[[j]]),
                                              fill = color,
                                              alpha = .6)
  }
  return(tree_visual)

}



################################################################################
## CHECK THAT THE drop.tip FUCNTION CALL ON LINE 454 IS USING ape OR tree.io !##
################################################################################

#############################################################
#NEED TO SPECIFY get_tip_level() FOR THE FUNCTION BELOW!!!!!#
#############################################################

#' Generate tree visuals highlighting specified row and column cluster from a heatmap.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param tree_object A phylo object representing the entire taxonomy from which the parameter 'tree' is derived.
#' @param htmap A ComplexHeatmap object.
#' @param row_cluster Index for row cluster of htmap to be illustrated.
#' @param column_cluster Index for column cluster of htmap to be illustrated.
#' @param isolate_subtree Alternate parameter for pruning tree diagram to create a second diagram.
#' @param show_labels Alternate parameter specifying whether to show tip labels.
#' @param show_clades Alternate parameter for labeling superclass clades.
#' @param highlight_clades Alternate parameter for labeling superclass clades.
#' @param point_size Alternate parameter for size of tip points of represented tips.
#' @param bar_size Alternate parameter for size of bars highlighting clades.
#' @return A ggtree object or list of ggtree objects.
#' @export
#' @import ComplexHeatmap
#' @import stats
#' @import ggplot2
#' @import ggtree
#' @import phangorn
#' @import ape
generate_tree_cluster <- function(tree, tree_object, htmap, row_cluster, column_cluster, row_name = 'Row data set', column_name = 'Column data set',  isolate_subtree = FALSE, show_labels = FALSE, show_clades = TRUE, highlight_clades = TRUE, point_size = 2, bar_size = 1){
  # get tree labels
  tree_labels <- c(tree$tip.label, tree$node.label)
  # get row labels
  row_labels <- dimnames(htmap@ht_list[[1]]@matrix)[[1]][stats::order.dendrogram(ComplexHeatmap::row_dend(htmap)[[row_cluster]])]
  # get column labels
  column_labels <- dimnames(htmap@ht_list[[1]]@matrix)[[2]][stats::order.dendrogram(column_dend(htmap)[[column_cluster]])]
  # get shared labels
  shared_labels <- intersect(row_labels, column_labels)
  # get row superclasses
  row_superclasses <- unique(unname(unlist(cluster_analysis(htmap = htmap, row_cluster = row_cluster, column_cluster = column_cluster, tree_object = tree_object, tree = tree)[1])))
  # get column superclasses
  column_superclasses <- unique(unname(unlist(cluster_analysis(htmap = htmap, row_cluster = row_cluster, column_cluster = column_cluster, tree_object = tree_object, tree = tree)[2])))
  # shared superclasses
  shared_superclasses <- intersect(row_superclasses, column_superclasses)
  # cut shared superclasses from row and column lists
  row_superclasses <- setdiff(row_superclasses, shared_superclasses)
  column_superclasses <- setdiff(column_superclasses, shared_superclasses)

  #print(row_superclasses)
  #print(which(tree_labels %in% row_superclasses))
  #print(column_superclasses)
  #print(which(tree_labels %in% column_superclasses))
  #print(shared_superclasses)
  #print(which(tree_labels %in% shared_superclasses))

  # select only colors that represent present tips and data sets
  ###superclass_lengths <- c(length(row_superclasses), length(column_superclasses), length(shared_superclasses))
  cluster_lengths <- c(length(intersect(setdiff(row_labels, shared_labels),tree_labels)), length(intersect(setdiff(column_labels, shared_labels),tree_labels)), length(intersect(shared_labels,tree_labels)))
  #print(superclass_lengths)
  color_selection <- which(cluster_lengths > 0)
  #color_selection <- which(superclass_lengths > 0)
  #print(color_selection)
  color_values <- c("row" = "#053061", "column" = "#d73027", "both" = "#2d004b")[color_selection]
  #print(color_values)
  color_labels <- c(row_name, column_name, paste(row_name, 'and', column_name))[color_selection]
  #print(color_labels)


  # build tree visual
  tree_visual <- ggtree(tree) +
    ggtree::layout_circular() +
    #geom_tiplab(size = 0.1) +
    ggtree::geom_point2(aes(subset = (label %in% intersect(setdiff(row_labels, shared_labels), c(tree$tip.label, tree$node.label))),
                    color = "row"),
                size = point_size) +
    ggtree::geom_point2(aes(subset = (label %in% intersect(setdiff(column_labels, shared_labels), c(tree$tip.label, tree$node.label))),
                    color = "column"),
                size = point_size) +
    ggtree::geom_point2(aes(subset = (label %in% intersect(shared_labels, c(tree$tip.label, tree$node.label))),
                    color = "both"),
                size = point_size) +
    ggtree::scale_color_manual(name = 'Data sets',
                       values= color_values,
                       labels= color_labels) +
    ggtree::theme(legend.position = c(1.2, 0.2))

  if (show_clades){
    if (length(row_superclasses) > 0){
      for (i in seq_along(row_superclasses)){
        if (row_superclasses[[i]] %in% tree_labels){
          tree_visual <- tree_visual + ggtree::geom_cladelab(node = which(tree_labels %in% row_superclasses[[i]]),
                                                     label = row_superclasses[[i]],
                                                     textcolor = '#2166ac',
                                                     barcolor = '#2166ac',
                                                     barsize = bar_size,
                                                     offset = .2,
                                                     offset.text = 3,
                                                     fontsize = 3.8,
                                                     angle = 'auto')
        } else {
          tree_visual <- handle_missing_node_show_clade(tree, tree_object, row_superclasses, tree_visual, i, '#2166ac')}
      }
    }

    if (length(column_superclasses) > 0){
      for (i in seq_along(column_superclasses)){
        if (column_superclasses[[i]] %in% tree_labels){
          tree_visual <- tree_visual + ggtree::geom_cladelab(node = which(tree_labels %in% column_superclasses[[i]]),
                                                     label = column_superclasses[[i]],
                                                     textcolor = '#b2182b',
                                                     barcolor = '#b2182b',
                                                     barsize = bar_size,
                                                     offset = .2,
                                                     offset.text = 3,
                                                     fontsize = 3.8,
                                                     angle = 'auto')
        } else {
          tree_visual <- handle_missing_node_show_clade(tree, tree_object, column_superclasses, tree_visual, i, '#b2182b')}
      }
    }

    if (length(shared_superclasses) > 0){
      for (i in seq_along(shared_superclasses)){
        if (shared_superclasses[[i]] %in% tree_labels){
          tree_visual <- tree_visual + ggtree::geom_cladelab(node = which(tree_labels %in% shared_superclasses[[i]]),
                                                     label = shared_superclasses[[i]],
                                                     textcolor = '#542788',
                                                     barcolor = '#542788',
                                                     barsize = bar_size,
                                                     offset = .2,
                                                     offset.text = 3,
                                                     fontsize = 3.8,
                                                     angle = 'auto')
        } else {
          tree_visual <- handle_missing_node_show_clade(tree, tree_object, shared_superclasses, tree_visual, i, '#542788')}
      }
    }
  }

  if (highlight_clades){
    if (length(row_superclasses) > 0){
      for (i in seq_along(row_superclasses)){
        if (row_superclasses[[i]] %in% tree_labels){
          tree_visual <- tree_visual + ggtree::geom_hilight(node = which(tree_labels %in% row_superclasses[[i]]),
                                                    fill = "#e6f5d0",
                                                    alpha = .6)
        } else {
          tree_visual <- handle_missing_node_highlight_clade(tree, tree_object, row_superclasses, tree_visual , i, "#e6f5d0")
        }
      }
    }
    if (length(column_superclasses) > 0){
      for (i in seq_along(column_superclasses)){
        if (column_superclasses[[i]] %in% tree_labels){
          tree_visual <- tree_visual + ggtree::geom_hilight(node = which(tree_labels %in% column_superclasses[[i]]),
                                                    fill = "#e0f3f8",
                                                    alpha = .6)
        } else {
          tree_visual <- handle_missing_node_highlight_clade(tree, tree_object, column_superclasses, tree_visual , i, "#e0f3f8")
        }
      }
    }
    if (length(shared_superclasses) > 0){
      for (i in seq_along(shared_superclasses)){
        if (shared_superclasses[[i]] %in% tree_labels){
          tree_visual <- tree_visual + ggtree::geom_hilight(node = which(tree_labels %in% shared_superclasses[[i]]),
                                                    fill = "#fde0ef",
                                                    alpha = .6)
        } else {
          tree_visual <- handle_missing_node_highlight_clade(tree, tree_object, shared_superclasses, tree_visual, i, "#fde0ef")
        }
      }
    }
  }

  # if isolating tree
  if (isolate_subtree) {
    superclasses <- unique(unname(unlist(cluster_analysis(htmap = htmap, row_cluster = row_cluster, column_cluster = column_cluster, tree_object = tree_object, tree = tree))))
    ancestors <- tree_labels[unique(unname(unlist(phangorn::Ancestors(tree, which(tree_labels %in% superclasses)))))]
    descendants <- tree_labels[unname(unlist(phangorn::Descendants(tree, which(tree_labels %in% superclasses),type = 'tips')))]
    #print(paste('Desc', length(descendants)))
    # handle cases when superclass nodes are not in the tree_label list
    missing_superclasses <- superclasses[-which(superclasses %in% tree_labels)]
    #print(superclasses)
    #print(superclasses[which(superclasses %in% tree_labels)])
    #print(missing_superclasses)
    for (l in seq_along(missing_superclasses)){
      temp_descendants <- intersect(c(tree_object$tip.label, tree_object$node.label)[phangorn::Descendants(tree_object, which(c(tree_object$tip.label, tree_object$node.label) %in% missing_superclasses[[l]]), type = 'all')], tree_labels)
      #print(paste('length of temp desc for l = ', l, ':', length(temp_descendants)))
      shallow_level <- min(sapply(temp_descendants, get_tip_level, tree = tree_object))
      #print(shallow_level)
      shallow_descendants <- temp_descendants[which(sapply(temp_descendants, get_tip_level, tree = tree_object) == shallow_level)]
      #print(paste('there are this many shallow descendants', length(shallow_descendants)))
      ancestors <- unique(c(ancestors, tree_labels[unique(unname(unlist(phangorn::Ancestors(tree, which(tree_labels %in% shallow_descendants)))))]))
      descendants <- unique(c(descendants, tree_labels[unname(unlist(phangorn::Descendants(tree, which(tree_labels %in% shallow_descendants),type = 'tips')))]))
    }
    #print(paste('Desc', length(descendants)))


    subtree <- ape::drop.tip(tree, setdiff(tree$tip.label, intersect(tree$tip.label, c(superclasses, ancestors, descendants))))
    # get labels
    subtree_labels <- c(subtree$tip.label, subtree$node.label)

    tree_visual_sub <- ggtree(subtree) +
      ggtree::layout_circular() +
      #geom_tiplab(size = 0.1) +
      ggtree::geom_point2(aes(subset = (label %in% intersect(setdiff(row_labels, shared_labels), c(subtree$tip.label, subtree$node.label))),
                      color = "row"),
                  size = point_size) +
      ggtree::geom_point2(aes(subset = (label %in% intersect(setdiff(column_labels, shared_labels), c(subtree$tip.label, subtree$node.label))),
                      color = "column"),
                  size = point_size) +
      ggtree::geom_point2(aes(subset = (label %in% intersect(shared_labels, c(subtree$tip.label, subtree$node.label))),
                      color = "both"),
                  size = point_size) +
      ggtree::scale_color_manual(name = 'Data sets',
                         values= color_values,
                         labels= color_labels)+
      ggtree::theme(legend.position = c(1.2, 0.2))

    if (show_clades){
      if (length(row_superclasses) > 0){
        for (i in seq_along(row_superclasses)){
          if (row_superclasses[[i]] %in% subtree_labels){
            tree_visual_sub <- tree_visual_sub + ggtree::geom_cladelab(node = which(subtree_labels %in% row_superclasses[[i]]),
                                                               label = row_superclasses[[i]],
                                                               textcolor = "#2166ac",
                                                               barcolor = "#2166ac",
                                                               barsize = bar_size,
                                                               offset = .2,
                                                               offset.text = 3,
                                                               fontsize = 3.8,
                                                               angle = 'auto')
          } else {
            tree_visual_sub <- handle_missing_node_show_clade(subtree, tree_object, row_superclasses, tree_visual_sub, i, "#2166ac")
          }
        }
      }

      if (length(column_superclasses) > 0){
        for (i in seq_along(column_superclasses)){
          if (column_superclasses[[i]] %in% subtree_labels){
            tree_visual_sub <- tree_visual_sub + ggtree::geom_cladelab(node = which(subtree_labels %in% column_superclasses[[i]]),
                                                               label = column_superclasses[[i]],
                                                               textcolor = '#b2182b',
                                                               barcolor = '#b2182b',
                                                               barsize = bar_size,
                                                               offset = .2,
                                                               offset.text = 3,
                                                               fontsize = 3.8,
                                                               angle = 'auto')
          } else {
            tree_visual_sub <- handle_missing_node_show_clade(subtree, tree_object, column_superclasses, tree_visual_sub, i, "#b2182b")
          }
        }
      }
      if (length(shared_superclasses) > 0){
        for (i in seq_along(shared_superclasses)){
          if (shared_superclasses[[i]] %in% subtree_labels){
            tree_visual_sub <- tree_visual_sub + ggtree::geom_cladelab(node = which(subtree_labels %in% shared_superclasses[[i]]),
                                                               label = shared_superclasses[[i]],
                                                               textcolor = '#542788',
                                                               barcolor = '#542788',
                                                               barsize = bar_size,
                                                               offset = .2,
                                                               offset.text = 3,
                                                               fontsize = 3.8,
                                                               angle = 'auto')
          } else {
            tree_visual_sub <- handle_missing_node_show_clade(subtree, tree_object, shared_superclasses, tree_visual_sub, i, "#542788")
          }
        }
      }
    }

    if (highlight_clades){
      if (length(row_superclasses) > 0){
        for (i in seq_along(row_superclasses)){
          if (row_superclasses[[i]] %in% subtree_labels) {
            tree_visual_sub <- tree_visual_sub + ggtree::geom_hilight(node = which(subtree_labels %in% row_superclasses[[i]]),
                                                              fill = '#e6f5d0',
                                                              alpha = .6)
          } else {
            tree_visual_sub <- handle_missing_node_highlight_clade(subtree, tree_object, row_superclasses, tree_visual_sub , i, "#e6f5d0")
          }
        }
      }
      if (length(column_superclasses) > 0){
        for (i in seq_along(column_superclasses)){
          if (column_superclasses[[i]] %in% subtree_labels){
            tree_visual_sub <- tree_visual_sub + ggtree::geom_hilight(node = which(subtree_labels %in% column_superclasses[[i]]),
                                                              fill = '#e0f3f8',
                                                              alpha = .6)
          } else {
            tree_visual_sub <- handle_missing_node_highlight_clade(subtree, tree_object, column_superclasses, tree_visual_sub , i, "#e0f3f8")
          }
        }
      }
      if (length(shared_superclasses) > 0){
        for (i in seq_along(shared_superclasses)){
          if (shared_superclasses[[i]] %in% subtree_labels){
            tree_visual_sub <- tree_visual_sub + ggtree::geom_hilight(node = which(subtree_labels %in% shared_superclasses[[i]]),
                                                              fill = '#fde0ef',
                                                              alpha = .6)
          } else {
            tree_visual_sub <- handle_missing_node_highlight_clade(subtree, tree_object, shared_superclasses, tree_visual_sub , i, "#fde0ef")
          }
        }
      }
    }
    if (show_labels & !show_clades) {
      tree_visual_sub <- tree_visual_sub + ggtree::geom_tiplab(size = 0.1) +
        ggtree::scale_color_manual(values=c(rep('#2166ac', length(row_superclasses)),
                                    rep('#b2182b', length(column_superclasses)),
                                    rep('#542788', length(shared_superclasses)), "#053061", "#d73027", "#2d004b"),
                           labels=c(row_superclasses, column_superclasses, shared_superclasses, 'TSCA', 'NSSS', 'TSCA and NSSS'))
    }

    return(list(tree_visual, tree_visual_sub))
  }

  return(tree_visual)

}


