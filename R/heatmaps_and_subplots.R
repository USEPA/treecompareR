is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#' Label numbers
#'
#' This is a helper function, use for determining number of labels in a given
#' data set per taxonomical level.
#'

#' @param datatable A data.table object of chemical classifications.
#' @param chemont `TRUE`/`FALSE`: whether to use ChemOnt taxonomy.
#' @param log `TRUE`/`FALSE`: whether numbers are reported as-is or as their log.
#' @return Number of occurrences of each label in the parameter `datatable`.
#' @import data.table

label_numbers <- function(datatable, chemont = TRUE, log = TRUE) {
  kingdom <- NULL
  superclass <- NULL
  class <- NULL
  subclass <- NULL
  level5 <- NULL
  level6 <- NULL
  level7 <- NULL
  level8 <- NULL
  level9 <- NULL
  level10 <- NULL
  level11 <- NULL
  DTXSID <- NULL
  . <- NULL
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
                                        by = .(DTXSID)][, c('DTXSID') := NULL])))
  #print(length(complete_labels))
  empty_indices <- which(sapply(complete_labels, function(t) {t == ''}))

  #print(length(empty_indices))
  if (length(empty_indices) > 0){
    complete_labels <- complete_labels[-which(sapply(complete_labels, function(t) {t == ''}))]
  }

  #print(length(complete_labels))
  na_indices <- which(sapply(complete_labels, is.na))
  #print(length(na_indices))
  if (length(na_indices > 0)){
    complete_labels <- complete_labels[-na_indices]
  }
  #complete_labels <- complete_labels[-which(sapply(complete_labels, is.na))]

  #print(length(complete_labels))
  unique_labels <- unique(complete_labels)

  if (log) {
    unique_complete_label_numbers <- sapply(unique_labels, function(t){
      log10(length(which(complete_labels %in% t)))
      })
    #print(length(unique_complete_label_numbers))
  } else {
    unique_complete_label_numbers <- sapply(unique_labels, function(t){
      length(which(complete_labels %in% t))
      })
  }
  return (list(unique_complete_label_numbers, complete_labels))
}


#' @title Generate similarity heatmap
#'
#' @description Generate a heatmap showing similarity of two classified data
#'   sets.
#'
#' @details This function takes in two `data.frame`s of classified entities,
#'   along with the taxonomy tree used to classify them, and a pre-computed
#'   matrix of pairwise similarities between all nodes of that tree. It produces
#'   a heatmap showing the pairwise similarities between the classifications of
#'   the two data sets, one data set on the rows and the other on the columns,
#'   annotated with bar graphs showing the number of occurrences of each label
#'   in each data set. The heatmap will be automatically clustered on rows and
#'   columns. Optionally, the user may specify that rows and columns should be
#'   split according to these clusters for plotting. See [the ComplexHeatmap
#'   reference](https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html)
#'   for more details.
#'
#'
#'
#' @param row_data A `data.frame` object of classified entities. Typically, each
#'   row represents one entity. Must include variables with names matching
#'   everything in `tax_level_labels`, containing the classification at each
#'   taxonomy level. May include a variable with unique identifiers for
#'   entities; if so, this must be the same variable name in both `row_data` and
#'   `column_data`, and it should be supplied in `entity_id_col`.
#' @param column_data A `data.frame` object of classified entities.  Typically,
#'   each row represents one entity. Must include variables with names matching
#'   everything in `tax_level_labels`, containing the classification at each
#'   taxonomy level. May include a variable with unique identifiers for
#'   entities; if so, this must be the same variable name in both `row_data` and
#'   `column_data`, and it should be supplied in `entity_id_col`.
#' @param terminal_only `TRUE`/`FALSE`: Whether to show only terminal
#'   classifications of row/column data in the heatmap (`TRUE`, default) or
#'   whether to include classifications at all levels (`FALSE`).
#' @param tree_object A `phylo`-class object representing a rooted tree. Default
#'   [chemont_tree].
#' @param matrix A matrix of pairwise similarity measure values derived from
#'   `tree_object`. Default [chemont_jaccard], to plot Jaccard similarity of
#'   taxonomic classification. Other pre-built option include
#'   [chemont_resnik_IC_SVH] (to plot Resnik similarity), [chemont_lin_IC_SVH]
#'   (to plot Lin similarity), and [chemont_jiangconrath_IC_SVH] (to plot
#'   Jiang-Conrath similarity).
#' @param row_indices Integer vector: A subset of row indices of the similarity
#'   matrix to consider. Default `NULL` to include all rows.
#' @param column_indices Integer vector: The column indices of the similarity
#'   matrix to consider. Default `NULL` to include all columns.
#' @param entity_id_col Character: the variable name(s) in `row_data` and
#'   `column_data` that uniquely identifies entities. Must be the same name(s)
#'   in each `data.frame`. Default `NULL` to assume each row is its own entity.
#' @param tax_level_labels Character: levels of the taxonomy in `tree_object`.
#'   Default [chemont_tax_levels].
#' @param name Name of the heatmap similarity measure for legend title. Default
#'   `'Similarity'`.
#' @param row_split Positive integer: Number of clusters to split rows into.
#'   Default `NULL`, to do no row splitting.
#' @param column_split Positive integer: Number of clusters to split columns
#'   into. Default `NULL`, to do no column splitting.
#' @param row_title Character: Title for rows. Default `'Row title'`.
#' @param column_title Character: Title for columns. Default `'Column title'`.
#' @param log_trans `TRUE`/`FALSE`: Whether to log-transform numbers of label
#'   occurrence for labels present in each data set before plotting them as bar
#'   chart annotations on the heatmap. Default `TRUE`.
#' @param colors A colormap used for the heatmap. Should be a function produced
#'   by [circlize::colorRamp2()], which accepts a vector of numeric values and
#'   returns interpolated colors. Default `circlize::colorRamp2(breaks = seq(0,
#'   1, len = 20), viridis::viridis(n=20, option = 'C'))` to use the
#'   [viridis::viridis()] colormap.
#' @return A [ComplexHeatmap::Heatmap()] object.
#' @examples
#' generate_heatmap(tree_object = chemont_tree,
#'  matrix = chemont_jaccard,
#'   row_data = biosolids_class,
#'    column_data = usgs_class,
#'    row_split = 5L,
#'    column_split = 5L,
#'     row_title = "Biosolids",
#'      column_title = "USGS Water",
#'    name = "Jaccard")
#'
#' @export
#' @import ComplexHeatmap
#' @importFrom magrittr `%>%`
#'
#' @seealso \code{\link{generate_tree_cluster}}
#'
generate_heatmap <- function(
    row_data,
    column_data,
    terminal_only = TRUE,
    tree_object = chemont_tree,
    matrix = chemont_jaccard,
    row_indices = NULL,
    column_indices = NULL,
    entity_id_col = NULL,
    tax_level_labels = chemont_tax_levels,
    name = 'Similarity',
    row_split = NULL,
    column_split = NULL,
    row_title = 'Row title',
    column_title = 'Column title',
    log_trans = TRUE,
    colors = circlize::colorRamp2(
      breaks = seq(0, 1, len = 20),
      colors = c('#440154FF',
                 '#481568FF',
                 '#482677FF',
                 '#453781FF',
                 '#3F4788FF',
                 '#39558CFF',
                 '#32648EFF',
                 '#2D718EFF',
                 '#287D8EFF',
                 '#238A8DFF',
                 '#1F968BFF',
                 '#20A386FF',
                 '#29AF7FFF',
                 '#3CBC75FF',
                 '#56C667FF',
                 '#74D055FF',
                 '#94D840FF',
                 '#B8DE29FF',
                 '#DCE318FF',
                 '#FDE725FF')
    )) {

  if (!identical(unlist(names(row_data)), unlist(names(column_data))))
      stop('The classification levels for the row data and column data do not match!')

  if (is.null(row_indices)){
    row_indices <- 1:dim(matrix)[[1]]
  }

  if (is.null(column_indices)){
    column_indices <- 1:dim(matrix)[[2]]
  }

  if(!is.null(row_split)){
    if(!is.numeric(row_split)){
      stop("row_split must be a positive integer")
    }
    if(length(row_split)>1){
      row_split <- row_split[1]
      warning("row_split of length > 1; taking the first element")
    }
    if(!is.finite(row_split)){
      stop("row_split is NA/NaN/Inf; it must be a positive integer")
    }
    if(row_split < 1){
      stop("row_split is zero or negative; it must be a positive integer")
    }
    if(!is.integer(row_split)){
    #attempt to coerce to integer
    #with a warning if necessary
      if(is.na(as.integer(row_split))){
        stop("row_split could not be coerced to integer")
      }
      if(!is.wholenumber(row_split)){
        warning(paste("non-integer row_split = ",
                      row_split,
                      "will be coerced to integer row_split = ",
                      as.integer(row_split)))
      }
    row_split <- as.integer(row_split)
    }
  }

  if(!is.null(column_split)){
    if(!is.numeric(column_split)){
      stop("column_split must be a positive integer")
    }
    if(length(column_split)>1){
      column_split <- column_split[1]
      warning("column_split of length > 1; taking the first element")
    }
    if(!is.finite(column_split)){
      stop("column_split is NA/NaN/Inf; it must be a positive integer")
    }
    if(column_split < 1){
      stop("column_split is zero or negative; it must be a positive integer")
    }
    if(!is.integer(column_split)){
      #attempt to coerce to integer
      #with a warning if necessary
      if(is.na(as.integer(column_split))){
        stop("column_split could not be coerced to integer")
      }
      if(!is.wholenumber(column_split)){
        warning(paste("non-integer column_split = ",
                      column_split,
                      "will be coerced to integer column_split = ",
                      as.integer(column_split)))
      }
      column_split <- as.integer(column_split)
    }
  }

  taxonomy_names <- names(row_data)
  # COLLECT LABEL NUMBERS FOR ROW DATA AND FOR COLUMN DATA

  if(terminal_only %in% TRUE){
    #check whether row data already has terminal labels
    if(!("terminal_label" %in% names(row_data))){
      row_data <- add_terminal_label(data = row_data,
                         entity_id_col = entity_id_col,
                         tax_level_names = tax_level_names)
    }
  row_label_data <-
    count_entities_per_label(data = row_data,
                             entity_id_col = entity_id_col,
                             tax_level_labels = c(tax_level_labels,
                                                  "terminal_label"))
  #keep only terminal labels
  row_label_data <- row_label_data[["terminal_label"]]


  }else{ #if terminal_only == FALSE
    row_label_data <-
      count_entities_per_label(data = row_data,
                               entity_id_col = entity_id_col,
                               tax_level_labels = tax_level_labels) %>%
      dplyr::bind_rows()
  }

  row_label_data <- row_label_data %>%
    dplyr::filter(!is.na(label))

  row_label_numbers <- row_label_data$n
  if(log_trans %in% TRUE){
    row_label_numbers <- log10(row_label_numbers)
  }

  row_labels <- row_label_data$label
  names(row_label_numbers) <- row_labels

  row_anno_indices <- match(dimnames(matrix)[[1]][row_indices], row_labels)
  row_na_indices <- which(sapply(row_anno_indices, is.na))
  if (length(row_na_indices) > 0){
    row_anno_indices <- row_anno_indices[-row_na_indices]
  }


  if(terminal_only %in% TRUE){
    #check whether column data already has terminal labels
    if(!("terminal_label" %in% names(column_data))){
      column_data <- add_terminal_label(data = column_data,
                                     entity_id_col = entity_id_col,
                                     tax_level_names = tax_level_names)
    }
    column_label_data <-
      count_entities_per_label(data = column_data,
                               entity_id_col = entity_id_col,
                               tax_level_labels = c(tax_level_labels,
                                                    "terminal_label"))
    #keep only terminal labels
    column_label_data <- column_label_data[["terminal_label"]]


  }else{ #if terminal_only == FALSE
    column_label_data <-
      count_entities_per_label(data = column_data,
                               entity_id_col = entity_id_col,
                               tax_level_labels = tax_level_labels) %>%
      dplyr::bind_rows()
  }

  column_label_data <- column_label_data %>%
    dplyr::filter(!is.na(label))
  column_label_numbers <- column_label_data$n
  if(log_trans %in% TRUE){
    column_label_numbers <- log10(column_label_numbers)
  }
  column_labels <- column_label_data$label
  names(column_label_numbers) <- column_labels

  column_anno_indices <- match(dimnames(matrix)[[2]][column_indices],
                               column_labels)
  column_na_indices <- which(sapply(column_anno_indices, is.na))
  if (length(column_na_indices) > 0){
    column_anno_indices <- column_anno_indices[-column_na_indices]
  }


  matrix_row_indices <- intersect(
    which(dimnames(matrix)[[1]] %in% row_labels),
    row_indices)
  matrix_column_indices <- intersect(
    which(dimnames(matrix)[[2]] %in% column_labels),
    column_indices)

  if(log_trans){
    row_anno_label <- 'log(row count) bars'
    col_anno_label <- 'log(col count) bars'
  } else {
    row_anno_label <- 'row count bars'
    col_anno_label <- 'col count bars'
  }

  heatmap <- ComplexHeatmap::Heatmap(
    matrix = matrix[
      matrix_row_indices,
      matrix_column_indices
    ],
    name = name,
    col = colors,

    # NEED TO ADD HELPER FUNCTIONS FOR THIS
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      col_log_count_bar = anno_barplot(
        column_label_numbers[
          column_labels[
            column_anno_indices
          ]
        ]
      ),
      annotation_name_rot = 45,
      annotation_label = col_anno_label,
      annotation_name_gp = grid::gpar(fontsize = 8)
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      row_log_count_bar = anno_barplot(
        row_label_numbers[
          row_labels[
            row_anno_indices
          ]
        ],
        axis_param = list(direction = 'reverse')),
      annotation_name_rot = 45,
      annotation_label = row_anno_label,
      annotation_name_gp = grid::gpar(fontsize = 8)
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = row_split,
    column_split = column_split
  )
  heatmap <- draw(heatmap,
                  row_title = row_title,
                  row_title_gp = grid::gpar(fontsize = 10,
                                            fontface = 'bold'),
                  column_title = column_title,
                  column_title_gp = grid::gpar(fontsize = 10,
                                               fontface = 'bold'))
}




#' Heatmap cluster analysis
#'
#' This is a helper function that returns superclasses or classes for specified
#' clusters in \code{\link{generate_tree_cluster}}. This is used for providing
#' visual identification of different clades based on the specified taxonomic
#' level.
#'
#' @param htmap A ComplexHeatmap object with hierarchical clustering, as created
#'   by [generate_heatmap()].
#' @param row_cluster Index for the row cluster.
#' @param column_cluster Index for the column cluster.
#' @param level Integer: the level of depth for labels. Default 2.
#' @param tree_object A `phylo` object representing a rooted tree, the taxonomy
#'   being investigated. Default [chemont_tree].
#' @param subtree Optional: Another `phylo` object representing a rooted tree, to
#'   restrict cluster analysis to a subtree. Default `NULL`.
#' @return A list of labels for the row and column clusters, based on the level
#'   specified.
#' @author Paul Kruse
#' @examples
#' #first generate a heatmap on which to do the cluster analysis
#' my_htmap <- generate_heatmap(
#'   row_data = biosolids_class,
#'    column_data = usgs_class,
#'    tree_object = chemont_tree,
#'  matrix = chemont_jaccard,
#'  row_split = 5L,
#'  column_split = 5L,
#'     row_title = "Biosolids",
#'      column_title = "USGS Water",
#'    name = "Jaccard")
#' #now do the cluster analysis
#' cluster_analysis(htmap = my_htmap,
#' row_cluster = 2L,
#' column_cluster = 2L,
#' level = 2,
#' tree_object = chemont_tree
#' )
#'
cluster_analysis <- function(htmap,
                             row_cluster,
                             column_cluster,
                             level = 2,
                             tree_object,
                             subtree = NULL){
  if (!is.null(subtree)){
    subtree_labels <- c(subtree$tip.label, subtree$node.label)
  }

  # get row levels for row cluster (restrict to subtree if subtree is given)
  row_names <- dimnames(htmap@ht_list[[1]]@matrix)[[1]][
    stats::order.dendrogram(
      ComplexHeatmap::row_dend(htmap)[[row_cluster]]
      )
    ]

  row_nodes <- get_node_from_label(label = row_names,
                                   tree = tree_object)

  #get ancestors at the specified level for the row indexes
  row_clades <- get_clade(node = row_nodes,
            tree = tree_object,
            level = level)

  row_clades <- unique(
    row_clades[!is.na(row_clades)]
  )

  row_clades <- get_label_from_node(node = row_clades,
                                    tree = tree_object)

  # get row levels for row cluster (restrict to tree if tree is given)
  column_names <- dimnames(htmap@ht_list[[1]]@matrix)[[2]][
    stats::order.dendrogram(
      ComplexHeatmap::column_dend(htmap)[[column_cluster]]
    )
  ]

  column_nodes <- get_node_from_label(label = column_names,
                                   tree = tree_object)

  #get ancestors at the specified level for the row indexes
  column_clades <- get_clade(node = column_nodes,
                          tree = tree_object,
                          level = level)

  column_clades <- unique(
    column_clades[!is.na(column_clades)]
  )

  column_clades <- get_label_from_node(node = column_clades,
                                    tree = tree_object)

  return(list('row_class' = row_clades,
              'column_class' = column_clades))

}



#' @title Handle missing node show clade
#'
#' @description
#' Display clade labels associated with missing
#' nodes after pruning a tree
#'
#' @details This is a helper function to display clade labels associated with missing
#' nodes after pruning a tree. Missing nodes occur when a node from an original
#' tree loses all but one child in the pruning process and is thus no longer
#' considered a node in the subtree.
#'
#' @param subtree A `phylo` object representing a rooted tree.
#' @param tree_object A `phylo` object representing the entire taxonomy from
#'   which the `subtree` is derived.
#' @param list_superclasses A list of superclasses to be labeled.
#' @param tree_visual A `ggtree` object that will have clades labeled.
#' @param i The index for which list_superclasses element will be labeled.
#' @param color A color for the clade label.
#' @return A ggtree object that will have the specified clade labeled.
#' @import ggplot2
#' @import ggtree
#' @author Paul Kruse
#'
#' @seealso \code{\link{handle_missing_node_highlight_clade}}
#'
handle_missing_node_show_clade <- function(subtree,
                                           tree_object,
                                           list_superclasses,
                                           tree_visual,
                                           i,
                                           color){
  subtree_labels <- c(subtree$tip.label,subtree$node.label)

  temp_descendants <- intersect(
    c(
      tree_object$tip.label,
      tree_object$node.label)[
        phangorn::Descendants(
          tree_object,
          which(
            c(tree_object$tip.label,
              tree_object$node.label) %in%
              list_superclasses[[i]]),
          type = 'all')],
    subtree_labels)

  shallow_level <- min(
    sapply(
      temp_descendants,
      get_node_level,
      tree = tree_object
    )
  )
  shallow_descendants <- temp_descendants[
    which(
      sapply(
        temp_descendants,
        get_node_level,
        tree = tree_object) == shallow_level
    )
  ]
  middle <- (length(shallow_descendants)+1)%/% 2
  #print(middle)
  for (j in seq_along(shallow_descendants)){

    tree_visual <- tree_visual +
      ggtree::geom_cladelab(node = which
                            (subtree_labels %in% shallow_descendants[[j]]
                            ),
                            label = ifelse(
                              j == middle,
                              list_superclasses[[i]],
                              ''),
                            textcolor = color,
                            barcolor = color,
                            offset = .2,
                            offset.text = 3,
                            fontsize = 2.6,
                            angle = 'auto')
  }
  return(tree_visual)

}

#' @title Handle missing node highlight clade
#'
#' @description
#' Highlight clades associated with missing nodes
#' after pruning a tree
#'
#' @details This is a helper function to highlight clades associated with missing nodes
#' after pruning a tree. Missing nodes occur when a node from an original tree
#' loses all but one child in the pruning process and is thus no longer
#' considered a node in the subtree.
#'
#' @param subtree A `phylo` object representing a rooted tree (pruned from a larger tree).
#' @param tree_object A `phylo` object representing the entire taxonomy from which
#'   parameter `subtree` is derived.
#' @param list_superclasses A list of superclasses to be highlighted.
#' @param tree_visual A `ggtree` object that will have clades highlighted.
#' @param i The index for which list_superclasses element will be highlighted.
#' @param color A color for the clade highlight.
#' @return A ggtree object that will have the specified clade highlighted.
#' @import ggplot2
#' @import ggtree
#' @author Paul Kruse
#'
#' @seealso \code{\link{handle_missing_node_show_clade}}
#'
handle_missing_node_highlight_clade <- function(subtree,
                                                tree_object,
                                                list_superclasses,
                                                tree_visual,
                                                i,
                                                color){
  subtree_labels <- c(subtree$tip.label, subtree$node.label)
  temp_descendants <- intersect(
    c(
      tree_object$tip.label,
      tree_object$node.label
    )
    [phangorn::Descendants(
      tree_object,
      which(
        c(
          tree_object$tip.label,
          tree_object$node.label
        ) %in% list_superclasses[[i]]
      ),
      type = 'all')
    ],
    subtree_labels)
  #print(temp_descendants)
  shallow_level <- min(
    sapply(
      temp_descendants,
      get_node_level,
      tree = tree_object))
  shallow_descendants <- temp_descendants[
    which(
      sapply(
        temp_descendants,
        get_node_level,
        tree = tree_object
      ) == shallow_level
    )
  ]
  middle <- (length(shallow_descendants)+1)%/% 2
  if (length(shallow_descendants) == 0) {
    print('whomp')
    return(tree_visual)
  }
  for (j in seq_along(shallow_descendants)){
    tree_visual <- tree_visual +
      ggtree::geom_hilight(node = which(
        subtree_labels %in% shallow_descendants[[j]]
      ),
      fill = color,
      alpha = .6)
  }
  return(tree_visual)

}


#' @title Generate tree from cluster
#'
#' @description Generates a tree diagram highlighting specified row and column
#' clusters from a similarity heatmap
#'
#' @details This function generates tree visuals highlighting specified row and column
#' clusters produced from \code{\link{generate_heatmap}}. The `tree` parameter
#' gives an underlying tree that will be used in the plots. A second plot may
#' also be produced that prunes away extraneous superclasses to highlight only
#' the portions of the tree with labels present.
#'
#' @param subtree A `phylo` object representing a rooted tree.
#' @param tree_object A `phylo` object representing the entire taxonomy from
#'   which the parameter `tree` is derived.
#' @param htmap A [ComplexHeatmap::Heatmap()] object produced by
#'   [generate_heatmap()].
#' @param row_cluster Index for row cluster of htmap to be illustrated.
#' @param column_cluster Index for column cluster of htmap to be illustrated.
#' @param row_name Character: name of row data set.
#' @param column_name Character: name of column data set.
#' @param isolate_subtree `TRUE`/`FALSE`: Whether to prune tree diagram to
#'   create a second diagram. Default `FALSE`.
#' @param show_labels `TRUE`/`FALSE`: whether to show tip labels. Default
#'   `FALSE`.
#' @param show_clades `TRUE`/`FALSE`: Whether to label superclass clades.
#'   Default `TRUE`.
#' @param highlight_clades `TRUE`/`FALSE`: Whether to highlight superclass
#'   clades. Default `TRUE`.
#' @param point_size Numeric: Size of tip points of represented
#'   tips. Default 2.
#' @param bar_size Numeric: Size of bars highlighting clades. Default 1.
#' @return A ggtree object or list of ggtree objects.
#' @examples
#' #first generate a heatmap to highlight
#' my_htmap <- generate_heatmap(
#'   row_data = biosolids_class,
#'    column_data = usgs_class,
#'    tree_object = chemont_tree,
#'  matrix = chemont_jaccard,
#'  row_split = 5L,
#'  column_split = 5L,
#'     row_title = "Biosolids",
#'      column_title = "USGS Water",
#'    name = "Jaccard")
#' #then generate a tree with highlights
#' generate_tree_cluster(htmap = my_htmap,
#' tree = chemont_tree,
#' )
#'
#' @export
#' @import ggplot2
#' @import ggtree
#'
#' @seealso \code{\link{generate_heatmap}}
#'
generate_tree_cluster <- function(subtree,
                                  tree_object,
                                  htmap,
                                  row_cluster,
                                  column_cluster,
                                  row_name = 'Row data set',
                                  column_name = 'Column data set',
                                  isolate_subtree = FALSE,
                                  show_labels = FALSE,
                                  show_clades = TRUE,
                                  highlight_clades = TRUE,
                                  point_size = 2,
                                  bar_size = 1){
  label <- NULL
  # get tree labels
  subtree_labels <- c(subtree$tip.label, subtree$node.label)
  # get row labels
  row_labels <- dimnames(
    htmap@ht_list[[1]]@matrix)[[1]][
      stats::order.dendrogram(
        ComplexHeatmap::row_dend(htmap)[[row_cluster]]
        )
      ]
  # get column labels
  column_labels <- dimnames(
    htmap@ht_list[[1]]@matrix
    )[[2]][
      stats::order.dendrogram(
        column_dend(htmap)[[column_cluster]]
        )
      ]
  # get shared labels
  shared_labels <- intersect(row_labels, column_labels)
  # get row superclasses
  row_superclasses <- unique(
    unname(
      unlist(
        cluster_analysis(htmap = htmap,
                         row_cluster = row_cluster,
                         column_cluster = column_cluster,
                         tree_object = tree_object,
                         subtree = subtree)[1]
      )
    )
  )
  # get column superclasses
  column_superclasses <- unique(
    unname(
      unlist(
        cluster_analysis(htmap = htmap,
                         row_cluster = row_cluster,
                         column_cluster = column_cluster,
                         tree_object = tree_object,
                         subtree = subtree)[2])
      )
    )
  # shared superclasses
  shared_superclasses <- intersect(row_superclasses, column_superclasses)
  # cut shared superclasses from row and column lists
  row_superclasses <- setdiff(row_superclasses, shared_superclasses)
  column_superclasses <- setdiff(column_superclasses, shared_superclasses)

  # select only colors that represent present tips and data sets
  cluster_lengths <- c(
    length(
      intersect(
        setdiff(row_labels,
                shared_labels),
        subtree_labels
        )
    ),
    length(
      intersect(
        setdiff(column_labels,
                shared_labels),
        subtree_labels
        )
    ),
    length(
      intersect(
        shared_labels,
        subtree_labels
        )
    )
  )
  color_selection <- which(cluster_lengths > 0)
  color_values <- c("row" = "#053061",
                    "column" = "#d73027",
                    "both" = "#2d004b")[color_selection]
  color_labels <- c(row_name,
                    column_name,
                    paste(row_name,
                          'and',
                          column_name)
                    )[color_selection]


  # build tree visual
  tree_visual <- ggtree(subtree) +
    ggtree::layout_circular() +
    ggtree::geom_point2(
      aes(
        subset = (
          label %in% intersect(setdiff(row_labels, shared_labels),
                               c(subtree$tip.label, subtree$node.label)
                               )
          ),
                    color = "row"),
                size = point_size) +
    ggtree::geom_point2(
      aes(
        subset = (
          label %in% intersect(setdiff(column_labels, shared_labels),
                               c(subtree$tip.label, subtree$node.label)
                               )
          ),
                    color = "column"),
                size = point_size) +
    ggtree::geom_point2(
      aes(
        subset = (
          label %in% intersect(shared_labels,
                               c(subtree$tip.label, subtree$node.label)
                               )
          ),
                    color = "both"),
                size = point_size) +
    ggtree::scale_color_manual(name = 'Data sets',
                       values= color_values,
                       labels= color_labels) +
    ggtree::theme(legend.position = c(1.2, 0.2))

  if (show_clades){
    if (length(row_superclasses) > 0){
      for (i in seq_along(row_superclasses)){
        if (row_superclasses[[i]] %in% subtree_labels){
          tree_visual <- tree_visual +
            ggtree::geom_cladelab(
              node = which(
                subtree_labels %in% row_superclasses[[i]]
                ),
                                                     label = row_superclasses[[i]],
                                                     textcolor = '#2166ac',
                                                     barcolor = '#2166ac',
                                                     barsize = bar_size,
                                                     offset = .2,
                                                     offset.text = 3,
                                                     fontsize = 3.8,
                                                     angle = 'auto'
              )
        } else {
          tree_visual <- handle_missing_node_show_clade(
            subtree,
            tree_object,
            row_superclasses,
            tree_visual,
            i, '#2166ac')}
      }
    }

    if (length(column_superclasses) > 0){
      for (i in seq_along(column_superclasses)){
        if (column_superclasses[[i]] %in% subtree_labels){
          tree_visual <- tree_visual +
            ggtree::geom_cladelab(node = which(
              subtree_labels %in% column_superclasses[[i]]
              ),
                                                     label = column_superclasses[[i]],
                                                     textcolor = '#b2182b',
                                                     barcolor = '#b2182b',
                                                     barsize = bar_size,
                                                     offset = .2,
                                                     offset.text = 3,
                                                     fontsize = 3.8,
                                                     angle = 'auto')
        } else {
          tree_visual <- handle_missing_node_show_clade(
            subtree,
            tree_object,
            column_superclasses,
            tree_visual,
            i,
            '#b2182b')
          }
      }
    }

    if (length(shared_superclasses) > 0){
      for (i in seq_along(shared_superclasses)){
        if (shared_superclasses[[i]] %in% subtree_labels){
          tree_visual <- tree_visual +
            ggtree::geom_cladelab(node = which(subtree_labels %in% shared_superclasses[[i]]),
                                                     label = shared_superclasses[[i]],
                                                     textcolor = '#542788',
                                                     barcolor = '#542788',
                                                     barsize = bar_size,
                                                     offset = .2,
                                                     offset.text = 3,
                                                     fontsize = 3.8,
                                                     angle = 'auto')
        } else {
          tree_visual <- handle_missing_node_show_clade(
            subtree,
            tree_object,
            shared_superclasses,
            tree_visual,
            i, '#542788')
          }
      }
    }
  }

  if (highlight_clades){
    if (length(row_superclasses) > 0){
      for (i in seq_along(row_superclasses)){
        if (row_superclasses[[i]] %in% subtree_labels){
          tree_visual <- tree_visual +
            ggtree::geom_hilight(node = which(subtree_labels %in% row_superclasses[[i]]),
                                                    fill = "#e6f5d0",
                                                    alpha = .6)
        } else {
          tree_visual <- handle_missing_node_highlight_clade(subtree,
                                                             tree_object,
                                                             row_superclasses,
                                                             tree_visual ,
                                                             i, "#e6f5d0")
        }
      }
    }
    if (length(column_superclasses) > 0){
      for (i in seq_along(column_superclasses)){
        if (column_superclasses[[i]] %in% subtree_labels){
          tree_visual <- tree_visual +
            ggtree::geom_hilight(node = which(subtree_labels %in% column_superclasses[[i]]),
                                                    fill = "#e0f3f8",
                                                    alpha = .6)
        } else {
          tree_visual <- handle_missing_node_highlight_clade(subtree,
                                                             tree_object,
                                                             column_superclasses,
                                                             tree_visual ,
                                                             i, "#e0f3f8")
        }
      }
    }
    if (length(shared_superclasses) > 0){
      for (i in seq_along(shared_superclasses)){
        if (shared_superclasses[[i]] %in% subtree_labels){
          tree_visual <- tree_visual +
            ggtree::geom_hilight(node = which(subtree_labels %in% shared_superclasses[[i]]),
                                                    fill = "#fde0ef",
                                                    alpha = .6)
        } else {
          tree_visual <- handle_missing_node_highlight_clade(subtree,
                                                             tree_object,
                                                             shared_superclasses,
                                                             tree_visual,
                                                             i,
                                                             "#fde0ef")
        }
      }
    }
  }

  tree_visual <- tree_visual + xlim(0, max(tree_visual$data$x) + 50)

  # if isolating tree
  if (isolate_subtree) {
    superclasses <- unique(
      unname(
        unlist(
          cluster_analysis(
            htmap = htmap,
            row_cluster = row_cluster,
            column_cluster = column_cluster,
            tree_object = tree_object,
            subtree = subtree)
          )
        )
      )
    ancestors <- subtree_labels[
      unique(
        unname(
          unlist(
            phangorn::Ancestors(
              subtree,
              which(subtree_labels %in% superclasses)
              )
            )
          )
        )
      ]
    descendants <- subtree_labels[
      unname(
        unlist(
          phangorn::Descendants(subtree,
                                which(subtree_labels %in% superclasses),
                                type = 'all')
          )
        )
      ]

    # handle cases when superclass nodes are not in the tree_label list
    missing_superclasses <- superclasses[-which(superclasses %in% subtree_labels)]
    for (l in seq_along(missing_superclasses)){
      temp_descendants <- intersect(
        c(tree_object$tip.label,
                                      tree_object$node.label)[
                                        phangorn::Descendants(
                                          tree_object,
                                          which(
                                            c(
                                              tree_object$tip.label,
                                              tree_object$node.label) %in%
                                              missing_superclasses[[l]]
                                            ),
                                          type = 'all'
                                          )
                                        ],
        subtree_labels
        )
      shallow_level <- min(sapply(temp_descendants,
                                  get_node_level,
                                  tree = tree_object))
      #print(shallow_level)
      shallow_descendants <- temp_descendants[
        which(
          sapply(
            temp_descendants,
            get_node_level,
            tree = tree_object) == shallow_level
          )
        ]
      ancestors <- unique(
        c(
          ancestors,
          tree_labels[
            unique(
              unname(
                unlist(
                  phangorn::Ancestors(
                    subtree,
                    which(subtree_labels %in% shallow_descendants)
                    )
                  )
                )
              )
            ]
          )
        )
      descendants <- unique(
        c(
          descendants,
          subtree_labels[
            unname(
              unlist(
                phangorn::Descendants(
                  tree,
                  which(subtree_labels %in% shallow_descendants),
                  type = 'all'
                )
              )
            )
          ]
        )
      )
    }

    subtree <- drop_tips_nodes(tree = subtree,
                               labels = c(superclasses,
                                          ancestors,
                                          descendants),
                               keep_descendants = FALSE)
    # get labels
    subtree_labels <- c(subtree$tip.label, subtree$node.label)

    tree_visual_sub <- ggtree(subtree) +
      ggtree::layout_circular() +
      ggtree::geom_point2(
        aes(subset = (
          label %in% intersect(setdiff(row_labels,
                                                                     shared_labels),
                                                             c(subtree$tip.label,
                                                               subtree$node.label)
                                                             )
                                        ),
                      color = "row"),
                  size = point_size) +
      ggtree::geom_point2(
        aes(
          subset = (
            label %in% intersect(
              setdiff(
                column_labels,
                shared_labels),
              c(subtree$tip.label, subtree$node.label)
              )
            ),
                      color = "column"),
                  size = point_size) +
      ggtree::geom_point2(
        aes(
          subset = (
            label %in% intersect(shared_labels,
                                 c(subtree$tip.label,
                                   subtree$node.label)
                                 )
            ),
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
            tree_visual_sub <- tree_visual_sub +
              ggtree::geom_cladelab(node = which(
                subtree_labels %in% row_superclasses[[i]]
              ),
              label = row_superclasses[[i]],
              textcolor = "#2166ac",
              barcolor = "#2166ac",
              barsize = bar_size,
              offset = .2,
              offset.text = 3,
              fontsize = 3.8,
              angle = 'auto')
          } else {
            tree_visual_sub <- handle_missing_node_show_clade(subtree,
                                                              tree_object,
                                                              row_superclasses,
                                                              tree_visual_sub,
                                                              i,
                                                              "#2166ac")
          }
        }
      }

      if (length(column_superclasses) > 0){
        for (i in seq_along(column_superclasses)){
          if (column_superclasses[[i]] %in% subtree_labels){
            tree_visual_sub <- tree_visual_sub +
              ggtree::geom_cladelab(node = which(
                subtree_labels %in% column_superclasses[[i]]
                ),
                                                               label = column_superclasses[[i]],
                                                               textcolor = '#b2182b',
                                                               barcolor = '#b2182b',
                                                               barsize = bar_size,
                                                               offset = .2,
                                                               offset.text = 3,
                                                               fontsize = 3.8,
                                                               angle = 'auto')
          } else {
            tree_visual_sub <- handle_missing_node_show_clade(
              subtree,
              tree_object,
              column_superclasses,
              tree_visual_sub,
              i,
              "#b2182b")
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

    tree_visual_sub <- tree_visual_sub + xlim(0, max(tree_visual_sub$data$x) + 50)

    return(list(tree_visual, tree_visual_sub))
  }

  return(tree_visual)

}


