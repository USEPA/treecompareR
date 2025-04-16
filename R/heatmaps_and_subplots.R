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
#'@param row_indices Integer vector: A subset of row indices of the similarity
#'  matrix to consider. Default `NULL` to include all rows.
#'@param column_indices Integer vector: The column indices of the similarity
#'  matrix to consider. Default `NULL` to include all columns.
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

  # get row levels for row cluster
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



