#' @title Is whole number
#'
#' @description Checks whether a number is a whole number (i.e., an integer)
#'
#' @details Adapted from the example in the documentation for
#'   [base::is.integer()].
#'
#' @examples
#' is.wholenumber(5) #returns TRUE
#' is.wholenumber(1.3) #returns FALSE
#' is.wholenumber(NA_real_) #returns NA
#' is.wholenumber(numeric(0)) #returns logical(0)

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  if(!is.numeric(x)){
    message("x is non-numeric")
    return(FALSE)
  }else{
    abs(x - round(x)) < tol
  }

}

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


#'@title Generate similarity heatmap
#'
#'@description Generate a heatmap showing similarity of two classified data
#'  sets.
#'
#'@details This function takes in two `data.frame`s of classified entities,
#'  along with the taxonomy tree used to classify them, and a pre-computed
#'  matrix of pairwise similarities between all nodes of that tree. It produces
#'  a heatmap showing the pairwise similarities between the classifications of
#'  the two data sets, one data set on the rows and the other on the columns,
#'  annotated with bar graphs showing the number of occurrences of each label in
#'  each data set. The heatmap will be automatically clustered on rows and
#'  columns.
#'
#'   See [the ComplexHeatmap
#'  reference](https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html)
#'  for more details.
#'
#'
#'  # Row and column splitting
#'
#'  Arguments `row_split` and `column_split` control row and column splitting,
#'  respectively, of the clustered heatmap. They may be provided in one of
#'  several ways.
#'
#'   - Integer of length 1: As for [ComplexHeatmap::Heatmap()]; `cutree()` will be applied to the row/column dendrogram resulting from clustering the matrix, with `k` equal to the supplied value.
#'   - A vector of categorical variables, or a `data.frame` of categorical variables: As for [ComplexHeatmap::Heatmap()]. The length of the vector, or the number of rows of the `data.frame`, must match the corresponding dimension of `matrix` (i.e., `nrow(matrix)` for `row_split`, and `ncol(matrix)` for `column_split`).
#'   - One of the taxonomy levels (i.e., one of `tax_level_labels`): Create a a character vector containing the label of the ancestor at the specified taxonomy level for each label in the row/column of `matrix`. Then this vector will be used as a vector of categorical variables for splitting, as described in [ComplexHeatmap::Heatmap()].
#'   - The string `"level"`: Create a character vector the same length as the corresponding dimension of `matrix`, containing the taxonomy level of each label in the row/column of `matrix`. Then this vector will be used as a vector of categorical variables for splitting, as described in [ComplexHeatmap::Heatmap()].
#'
#'  # Extracting matrix
#'
#'  You can access the matrix that was ultimately plotted in the heatmap.
#'
#' ```
#' #generate a heatmap
#' my_htmap <- generate_heatmap(
#'    tree_object = chemont_tree,
#'    matrix = chemont_jaccard,
#'    row_data = biosolids_class,
#'    column_data = usgs_class,
#'    terminal_only = TRUE)
#'
#' #extract the heatmap matrix
#' my_matrix <- my_htmap@ht_list[[1]]@matrix
#' #compare this to chemont_jaccard!
#' ```
#'
#'  # How to add more annotations after the fact
#'
#'  If you would like to add your own annotations after the fact, you may want
#'  to use the option `draw = FALSE`. This returns a
#'  [ComplexHeatmap::HeatmapList()] object that has not had the
#'  `draw()` method applied to it yet. This means
#'  you can still concatenate it with other heatmaps and heatmap annotations
#'  using `+` and `%v%`. (If you use `draw =
#'  TRUE`, you won't be able to concatenate the result with any other heatmaps
#'  or annotations!)
#'
#' ```
#'   my_htmap <- generate_heatmap(
#'    tree_object = chemont_tree,
#'    matrix = chemont_jaccard,
#'    row_data = biosolids_class,
#'    column_data = usgs_class,
#'    terminal_only = TRUE,
#'    row_split = "superclass",
#'    column_split = "superclass",
#'     draw = FALSE,
#'      row_title_rot = 0,
#'      column_title_rot = 90,
#'       row_title_gp = grid::gpar(fontsize = 6),
#'        column_title_gp = grid::gpar(fontsize = 6),
#'         row_gap = unit(0.1, "mm"),
#'         column_gap = unit(0.1, "mm"))
#'
#'     #extract matrix
#'     #notice difference in syntax from above, now that draw = FALSE!
#'     my_matrix <- my_htmap@matrix
#'
#'     #construct a row annotation using information from my_matrix
#'     #this is just 1:nrow(my_matrix) in order
#'     my_row_ann <- ComplexHeatmap::rowAnnotation(foo = 1:nrow(my_matrix))
#'
#'     #add the new row annotation
#'     my_row_ann + my_htmap
#'     #note that the row annotation is automatically reordered
#'     # to match the clustering of the matrix!
#' ```
#'
#'@param row_data A `data.frame` object of classified entities. Typically, each
#'  row represents one entity. Must include variables with names matching
#'  everything in `tax_level_labels`, containing the classification at each
#'  taxonomy level. May include a variable with unique identifiers for entities;
#'  if so, this must be the same variable name in both `row_data` and
#'  `column_data`, and it should be supplied in `entity_id_col`.
#'@param column_data A `data.frame` object of classified entities.  Typically,
#'  each row represents one entity. Must include variables with names matching
#'  everything in `tax_level_labels`, containing the classification at each
#'  taxonomy level. May include a variable with unique identifiers for entities;
#'  if so, this must be the same variable name in both `row_data` and
#'  `column_data`, and it should be supplied in `entity_id_col`.
#'@param terminal_only `TRUE`/`FALSE`: Whether to show only terminal
#'  classifications of row/column data in the heatmap (`TRUE`, default) or
#'  whether to include classifications at all levels (`FALSE`).
#'@param tree_object A `phylo`-class object representing a rooted tree. Default
#'  [chemont_tree].
#'@param matrix A matrix of pairwise similarity measure values derived from
#'  `tree_object`. Default [chemont_jaccard], to plot Jaccard similarity of
#'  taxonomic classification. Other pre-built option include
#'  [chemont_resnik_IC_SVH] (to plot Resnik similarity), [chemont_lin_IC_SVH]
#'  (to plot Lin similarity), and [chemont_jiangconrath_IC_SVH] (to plot
#'  Jiang-Conrath similarity).
#'@param prune_to Optional: A pruning specification, as for [prune_tree()] (see
#'  the documentation for that function for acceptable pruning specifications).
#'  Default `NULL` to do no pruning. If `prune_to` is non-`NULL`, then the
#'  heatmap will include only labels that are in the pruned tree.
#'@param entity_id_col Character: the variable name(s) in `row_data` and
#'  `column_data` that uniquely identifies entities. Must be the same name(s) in
#'  each `data.frame`. Default `NULL` to assume each row is its own entity.
#'@param tax_level_labels Character: levels of the taxonomy in `tree_object`.
#'  Default [chemont_tax_levels].
#'@param name Name of the heatmap similarity measure for legend title. Default
#'  `'Similarity'`.
#'@param row_split A specification for row splitting; see Details. Default
#'  `NULL`, to do no row splitting.
#'@param column_split A specification for column splitting; see Details. Default
#'  `NULL`, to do no column splitting.
#'@param row_title Character: Title for rows. Default `'Row title'`.
#'@param column_title Character: Title for columns. Default `'Column title'`.
#'@param annotation_count `TRUE`/`FALSE`: Whether to add row and column
#'  annotations containing bar charts of the count of entities for each label.
#'  Default `TRUE`.
#'@param log_trans `TRUE`/`FALSE`: Whether to log-transform numbers of label
#'  occurrence for labels present in each data set before plotting them as bar
#'  chart annotations on the heatmap. Default `FALSE`.
#'@param colors A colormap used for the heatmap. Should be a function produced
#'  by [circlize::colorRamp2()], which accepts a vector of numeric values and
#'  returns interpolated colors. Default `circlize::colorRamp2(breaks = seq(0,
#'  1, len = 20), viridis::viridis(n=20, option = 'C'))` to use the
#'  [viridis::viridis()] colormap.
#'@param ... Other arguments to be passed to [ComplexHeatmap::Heatmap()].
#'@return A [ComplexHeatmap::Heatmap()] object.
#' @examples
#' #heatmap of biosolids vs. USGS water terminal labels
#' generate_heatmap(tree_object = chemont_tree,
#'  matrix = chemont_jaccard,
#'   row_data = biosolids_class,
#'    column_data = usgs_class,
#'    terminal_only = TRUE,
#'    #split rows/columns by superclasses
#'    row_split = "superclass",
#'    column_split = "superclass",
#'     row_title = "Biosolids",
#'      column_title = "USGS Water",
#'    name = "Jaccard",
#'    #ComplexHeatmap::Heatmap() args follow
#'     row_title_rot = 0,
#'      column_title_rot = 90,
#'       row_title_gp = grid::gpar(fontsize = 6),
#'        column_title_gp = grid::gpar(fontsize = 6),
#'         row_gap = unit(0.1, "mm"),
#'         column_gap = unit(0.1, "mm")
#'         )
#'
#' #an example of pruning:
#' #the following plots only the block of the previous heatmap
#' # labeled "Organoheterocyclic compounds" on both row and column
#' generate_heatmap(tree_object = chemont_tree,
#'  matrix = chemont_jaccard,
#'  prune_to = "Organoheterocyclic compounds",
#'   row_data = biosolids_class,
#'    column_data = usgs_class,
#'    terminal_only = TRUE,
#'    #since all are in the same superclass, split by class this time
#'    row_split = "class",
#'    column_split = "class",
#'    row_title = "Biosolids",
#'    column_title = "USGS Water",
#'    name = "Jaccard",
#'    #ComplexHeatmap::Heatmap() args follow
#'    row_title_rot = 0,
#'    column_title_rot = 90,
#'    row_title_gp = grid::gpar(fontsize = 6),
#'    column_title_gp = grid::gpar(fontsize = 6),
#'    row_gap = unit(0.1, "mm"),
#'    column_gap = unit(0.1, "mm"))
#'
#' #heatmap of similarity at all levels of classification
#' #with rows/column split by taxonomic level
#' generate_heatmap(tree_object = chemont_tree,
#'  matrix = chemont_jaccard,
#'  prune_to = "Organoheterocyclic compounds",
#'   row_data = biosolids_class,
#'    column_data = usgs_class,
#'    #include labels at all levels of classification
#'    terminal_only = FALSE,
#'    row_split = "level",
#'    column_split = "level",
#'    row_title = "Biosolids",
#'    column_title = "USGS Water",
#'    name = "Jaccard",
#'    #ComplexHeatmap::Heatmap() args follow
#'    row_title_rot = 0,
#'    column_title_rot = 90,
#'    row_title_gp = grid::gpar(fontsize = 10),
#'    column_title_gp = grid::gpar(fontsize = 10),
#'    row_gap = unit(0.1, "mm"),
#'    column_gap = unit(0.1, "mm"))
#'
#'@export
#'@import ComplexHeatmap
#'@importFrom magrittr `%>%`
#'
#'@seealso \code{\link{generate_tree_cluster}}
#'
generate_heatmap <- function(
    row_data,
    column_data,
    terminal_only = TRUE,
    tree_object = chemont_tree,
    matrix = chemont_jaccard,
    prune_to = NULL,
    entity_id_col = NULL,
    tax_level_labels = chemont_tax_levels,
    name = 'Similarity',
    row_split = NULL,
    column_split = NULL,
    row_title = 'Row title',
    column_title = 'Column title',
    annotation_count = TRUE,
    log_trans = FALSE,
    draw = TRUE,
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
    ),
    ...) {

  #row_split and column_split may be any of the following:
  #
  #a single integer -- in which case it will be treated as cuttree with that
  #argument
  #
  #the name of a taxonomy level -- in which case the dimension will be assigned
  #its values for that taxonomy level and split accordingly
  #
  #the special value "level" -- in which case each label in rows/cols of matrix
  #will be assigned its taxonomy level and split accordingly (e.g. all kingdoms
  #together, all superclasses together, etc.)
  #
  #some vector of categorical variable the same length and order as the
  #corresponding dimension of `matrix` -- in which case it will be subset the
  #same way as `matrix` and used accordingly.
  #
  #a data.frame with the same number of rows as the corresponding dimension of
  #`matrix`, where the variables are categorical and combinations of them
  #uniquely define groups for splitting.

    row_indices <- 1:dim(matrix)[[1]]
    column_indices <- 1:dim(matrix)[[2]]


  ##################################
  # Prune tree & matrix if requested
  ###################################
  if(!is.null(prune_to)){
    tree_object <- prune_tree(tree = tree_object,
                              prune_to = prune_to,
                              tax_level_labels = tax_level_labels,
                              entity_id_col = entity_id_col)
    #subset the matrix to keep only the labels in the pruned tree
    tree_labels <- c(tree_object$tip.label,
                     tree_object$node.label)
    #get row and column indices of matrix that are in the pruned tree
    #but don't actually subset matrix yet --
    #that will be done later
    row_indices <- which(rownames(matrix) %in% tree_labels)
    column_indices <- which(colnames(matrix) %in% tree_labels)
  }




  ######################
  # Set up row indexes
  #######################
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

  ######################
  # Set up column indexes
  #######################
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

  ######################
  # Get matrix row and column indexes
  #######################
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



  #handle row splitting
  row_split_err <- paste(
    "row_split should be a single integer,",
    "the name of a taxonomy level (one of tax_level_labels),",
    "a vector the same length as nrow(matrix)",
    "that can be coerced to a factor,",
    "or a data.frame with the same number of rows as",
    "nrow(matrix) containing categorical variables",
    "(that can be coerced to factor)",
    "whose combinations define groups for splitting the rows."
  )

row_split <- htmap_split_check(htmap_split = row_split,
                               matrix_dim = nrow(matrix),
                               matrix_inds = matrix_row_indices,
                               matrix_names = rownames(matrix),
                               tree_object = tree_object,
                               tax_level_labels = tax_level_labels,
                               err_msg = row_split_err)

column_split_err <- paste(
  "column_split should be a single integer,",
  "the name of a taxonomy level (one of tax_level_labels),",
  "a vector the same length as ncol(matrix)",
  "that can be coerced to a factor,",
  "or a data.frame with the same number of rows as",
  "ncol(matrix) containing categorical variables",
  "(that can be coerced to factor)",
  "whose combinations define groups for splitting the columns."
)
column_split <- htmap_split_check(htmap_split = column_split,
                               matrix_dim = ncol(matrix),
                               matrix_inds = matrix_column_indices,
                               matrix_names = colnames(matrix),
                               tree_object = tree_object,
                               tax_level_labels = tax_level_labels,
                               err_msg = column_split_err)

  ######################
  # Create Heatmap
  #######################


if(annotation_count %in% TRUE){
    ######################
    # Add top (column) annotation
    #######################
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      col_log_count_bar = ComplexHeatmap::anno_barplot(
        column_label_numbers[
          column_labels[
            column_anno_indices
          ]
        ]
      ),
      annotation_name_rot = 45,
      annotation_label = col_anno_label,
      annotation_name_gp = grid::gpar(fontsize = 8)
    )

    ######################
    # Add left (row) annotation
    #######################
    left_annotation <- ComplexHeatmap::rowAnnotation(
      row_log_count_bar = ComplexHeatmap::anno_barplot(
        row_label_numbers[
          row_labels[
            row_anno_indices
          ]
        ],
        axis_param = list(direction = 'reverse')),
      annotation_name_rot = 45,
      annotation_label = row_anno_label,
      annotation_name_gp = grid::gpar(fontsize = 8)
    )

    #heatmap with annotations
    heatmap <- ComplexHeatmap::Heatmap(
      matrix = matrix[
        matrix_row_indices,
        matrix_column_indices
      ],
      name = name,
      col = colors,
      show_row_names = FALSE,
      show_column_names = FALSE,
      row_split = row_split,
      column_split = column_split,
      top_annotation = top_annotation,
      left_annotation = left_annotation,
      ...
    )
}else{
  #heatmap without annotations
  heatmap <- ComplexHeatmap::Heatmap(
    matrix = matrix[
      matrix_row_indices,
      matrix_column_indices
    ],
    name = name,
    col = colors,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = row_split,
    column_split = column_split,
    ...
  )
}


  ######################
  # Draw heatmap to lock in row/column ordering
  # See https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#get-orders-and-dendrograms-from-heatmap
  #######################
if(draw %in% TRUE){
  heatmap <- draw(heatmap,
                  row_title = row_title,
                  row_title_gp = grid::gpar(fontsize = 10,
                                            fontface = 'bold'),
                  column_title = column_title,
                  column_title_gp = grid::gpar(fontsize = 10,
                                               fontface = 'bold'))
}

  return(heatmap)
}


#'@title Heatmap cluster tree
#'
#'@description Tree visualization showing membership in a heatmap cluster
#'
#'@details Given a heatmap produced using [generate_heatmap()], with clusters
#'  produced using some row/column splitting specifications (`row_split` and
#'  `column_split` arguments to [generate_heatmap()]), produce a tree diagram
#'  with branches highlighted by their membership in the heatmap row labels for
#'  that cluster, and the heatmap column labels for that cluster.
#'
#'@param htmap A [ComplexHeatmap::HeatmapList()] object produced using
#'  [generate_heatmap()] with some row/column split specifications.
#'@param cluster_row Integer: the row index for the cluster of interest.
#'@param cluster_column: Integer: the column index for the cluster of interest.
#'@param tree_obj The `phylo`-class tree object used to construct the heatmap
#'  (i.e., the same `tree_obj` that was supplied to [generate_heatmap()]).
#'@param prune_to A pruning specification to apply to the tree visualization
#'  (see [prune_tree()] for options). Default `NA`, which will prune the tree to
#'  show only the union of the row and column labels in the heatmap. Use
#'  `prune_to = NULL` (*not* the default!) if you do not want to prune the base
#'  tree.
#'@param prune_args A list of additional arguments to [prune_tree()]. Default
#'  `list(keep_descendants = FALSE, adjust_branch_length = TRUE)`.
#'@param ... Other arguments as for [display_subtree()], excluding `base_tree`
#'  `prune_to`, and `prune_args` (which are already supplied to this function as
#'  `tree_obj`, `prune_to`, and `prune_args`), `highlight_by` (which will always
#'  be `"set"` for this function), and `data_1`, `data_2`, `name_1`, and `name_2`
#'  (which will be automatically determined from the heatmap row/column labels
#'  and row/column titles).
#'@return A [ggtree::ggtree()] object produced by [display_subtree()].
#'
#'@export
#'
#'@examples
#'#generate a heatmap with clusters
#' my_htmap <- generate_heatmap(
#'tree_object = chemont_tree,
#'matrix = chemont_jaccard,
#'row_data = biosolids_class,
#'column_data = usgs_class,
#'terminal_only = TRUE,
#'row_split = 6L,
#'column_split = 6L,
#'row_title = "Biosolids",
#'column_title = "USGS Water",
#'row_gap = unit(0.5, "mm"),
#'column_gap = unit(0.5, "mm"))
#'
#'#now show a tree higlighted by membership in cluster 1/1 (high similarity)
#'
#'cluster_tree(htmap = my_htmap,
#'cluster_row = 1,
#'cluster_column = 1,
#'tree_obj = chemont_tree,
#'show_tiplabs = FALSE,
#'clade_level = 2,
#'bg_tree_scale = 0,
#'base_opts = list(size = 1))
#'
#'#compare to a tree highlighted by row/column labels in cluster 2/3 (some high and some low similarity)
#'
#'#'cluster_tree(htmap = my_htmap,
#'cluster_row = 2,
#'cluster_column = 3,
#'tree_obj = chemont_tree
#'show_tiplabs = FALSE,
#'clade_level = 2,
#'bg_tree_scale = 0,
#'base_opts = list(size = 1))
#'
cluster_tree <- function(htmap,
                         cluster_row,
                         cluster_column,
                         tree_obj,
                         prune_to = NA,
                         prune_args = list(keep_descendants = FALSE,
                                           adjust_branch_length = TRUE),
                         ...){
  htmap <- draw(htmap)
  #To get labels by cluster:
  row_inds <- row_order(htmap)

  #Convert row indexes into labels:
  #get matrix
  my_mat <- htmap@ht_list[[1]]@matrix
  #pull rownames by cluster
  row_labs <- lapply(row_inds,
                     function(this_cluster){
                       rownames(my_mat)[this_cluster]
                     })

  #similarly for columns
  column_inds <- column_order(htmap)

  #convert column indexes into labels
  column_labs <- lapply(column_inds,
                        function(this_cluster){
                          colnames(my_mat)[this_cluster]
                        })

  name_1 <- paste(htmap@row_title, "cluster",
                  paste0(cluster_row, "/", cluster_column))
  name_2 <- paste(htmap@column_title, "cluster",
                  paste0(cluster_row, "/", cluster_column))

  if(is.na(prune_to)){
    #by default, prune to labels in the heatmap
    prune_to <- do.call(union, dimnames(my_mat))
  }

  display_subtree(base_tree = tree_obj,
                  prune_to = prune_to,
                  prune_args = prune_args,
                  data_1 = row_labs[[cluster_row]],
                  name_1 = name_1,
                  data_2 = column_labs[[cluster_column]],
                  name_2 = name_2,
                  highlight_by = "set",
                  ...)

}

#' @title Heatmap split check
#'
#' @description A helper function for [generate_heatmap()] to check the
#'   row/column split specifications.
#'
#' @details See [generate_heatmap()] for possible heatmap split specifications.
#'
#' @param htmap_split The split specification (either `row_split` or
#'   `column_split` as provided to [generate_heatmap()]).
#' @param matrix_dim Integer: The size of the appropriate dimension of the full
#'   matrix to be heatmapped (i.e., argument `matrix` as supplied to
#'   [generate_heatmap()]): `nrow(matrix)` for `row_split`, `ncol(matrix)` for
#'   `column_split`.
#' @param matrix_inds Integer vector: The row or column indices of the full
#'   matrix to be heatmapped (i.e., argument `matrix` as supplied to
#'   [generate_heatmap()]).
#' @param matrix_names Character vector: The row or column names of the full
#'   matrix (i.e., argument `matrix` as supplied to [generate_heatmap()]).
#' @param tree_object The tree object from which the matrix was derived (i.e.,
#'   argument `tree_object` as supplied to [generate_heatmap()]). May be a
#'   pruned subset of the original tree (i.e., if [generate_heatmap()] was
#'   called with a non-`NULL` value for argument `prune_to`, then this will be
#'   `prune_tree(tree = tree_obj, prune_to = prune_to)`).
#' @param tax_level_labels Character vector of taxonomy level labels, as
#'   supplied to [generate_heatmap()].
#' @param err_msg Character: An error message to be thrown if the `htmap_split`
#'   specification is not one of the accepted types. Default `"Invalid heatmap
#'   split specification."`, but [generate_heatmap()] internally supplies its
#'   own, more-informative error messages.
#' @return If it does not stop with an error, returns the final row/column split
#'   specification. If split specification was `"level"` or the name of a
#'   taxonomy level, then the returned object will be a character vector the
#'   same length as `matrix_inds`. Otherwise, it will be the same class as
#'   `htmap_split`, but with the same length or number of rows as `matrix_inds`,
#'   unless `htmap_split` was a single integer, in which case the returned
#'   object is the same single integer.
#' @examples
#' #No examples. See ?generate_heatmap for acceptable split specifications.
#' @author Caroline Ring
#'
htmap_split_check <- function(htmap_split,
                              matrix_dim,
                              matrix_inds,
                              matrix_names,
                              tree_object,
                              tax_level_labels,
                              err_msg = "Invalid heatmap split specification."){

  if(!is.null(htmap_split)){
    if(is.data.frame(htmap_split)){
      if(!(nrow(htmap_split) %in% matrix_dim)){
        stop(err_msg)
      }else{
        #reorder/subset to match the matrix
        htmap_split <- htmap_split[matrix_inds, ]
      }
    }else{ #if it's not a data.frame
      if(length(htmap_split) %in% 1){ #if it's length 1
        if(is.numeric(htmap_split)){ #if numeric
          if(!is.wholenumber(htmap_split)){ #if not integer
            stop(err_msg)
          }
        }else{ #if length 1 but not numeric
          if(!is.character(htmap_split)){ #if not character either
            stop(err_msg) #then stop
          }else{ #if length 1 and character
            if(htmap_split %in% tax_level_labels){
              #get the corresponding ancestor for rows of matrix
              nodes_row <- get_node_from_label(
                label = matrix_names[matrix_inds],
                tree = tree_object)
              htmap_split <- get_label_from_node(
                node = get_clade(node = nodes_row,
                                 tree = tree_object,
                                 level = match(htmap_split, tax_level_labels)
                ),
                tree = tree_object)
            }else{ #if character, but not in tax_level_labels
              if(htmap_split %in% "level"){
                #get taxonomic level of each row label in matrix
                htmap_split <- tax_level_labels[
                  get_node_level(tree = tree_object,
                                 node = matrix_names[matrix_inds])
                ]
              }else{
                #if single character, but neither one of the taxonomy levels nor
                #the word "level"
                stop(err_msg)
              }
            } #end if character, but not in tax_level_labels
          } #end if length 1 and character
        } #end if length 1 but not numeric
      }else{ #if length(htmap_split) > 1
        #check and make sure it's the same length as matrix_dim
        if(!(length(htmap_split) %in% matrix_dim)){
          stop(err_msg)
        }else{
          htmap_split <- htmap_split[matrix_inds]
        }
      }
    }
  }

  return(htmap_split)
}



