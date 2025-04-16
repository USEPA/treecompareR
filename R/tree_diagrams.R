#' Label bars
#'
#' Bar plot of numbers of labels by taxonomy level
#'
#' This function takes in a (list of) `data.frame`(s) and returns figures
#' illustrating the label numbers by taxonomy level.
#'
#' @param data A `data.frame` or list of `data.frames` with chemicals and their
#'   classification data.
#' @param tax_level_labels Character vector of taxonomy levels. Default [chemont_tax_levels].
#' @return A list of two plots, the first showing number of labels by taxonomy
#'   level for each `data.frame` and the second showing the number of labels by
#'   data.table for each taxonomy level.
#' @author Paul Kruse
#' @examples
#' label_bars(data = biosolids_class[1:10, ])
#'
#' @export
#' @import ggplot2
label_bars <- function(data = NULL,
                       tax_level_labels = chemont_tax_levels){
  tax_levels <- NULL
  count_sums <- NULL
  dataset <- NULL

  if (is.null(data))
    stop('Please input data!')

  if (is.data.frame(data)){
    data <- list(data)
  } else {
    if (!is.list(data) | !all(sapply(data, function(t) {
      is.data.frame(t)})))
      stop('Please input a single data.frame or list of data.frames!')
  }

  data <- lapply(data, as.data.frame)

  number <- length(data)

  if (number == 1){
    data_names <- c('Set_1')
  } else {
    if (is.null(names(data)) |
        length(names(data)) != number |
        any(is.na(names(data))) |
        any(names(data) == '')
    ){
      data_names <- names(data)
      missing_names <- which(names(data) == '')
      data_names[missing_names] <- paste0('Set_', missing_names)
      warning("There were missing names! Attaching substitute names...")
    } else {
      data_names <- names(data)
    }
  }

  names(data) <- data_names


  df <- data.frame(
    tax_level_labels,
    unname(
      sapply(
        data,
        function(t) count_labels(
          data = t,
          tax_level_labels = tax_level_labels
        )
      )
    )
  )
  names(df) <- c('tax_levels', data_names)
  transformed_df <- df %>%
    tidyr::pivot_longer(!tax_levels,
                        names_to = 'dataset',
                        values_to = 'count_sums')
  transformed_df["count_sums"] <- as.numeric(transformed_df$count_sums)
  transformed_df["tax_levels"] <- factor(transformed_df$tax_levels,
                                         levels = tax_level_labels)
  transformed_df['dataset'] <- factor(transformed_df$dataset,
                                      levels = data_names)

  plot_1 <- ggplot(transformed_df) +
    facet_wrap(~dataset, scales = 'free') +
    geom_bar(aes(x = tax_levels,
                 y = count_sums,
                 fill = tax_levels),
             stat = 'identity') +
    scale_color_manual(name = 'Taxonomy level') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = 'Taxonomy levels',
         y = 'Number of unique labels',
         fill = 'Taxonomy levels')

  plot_2 <- ggplot(transformed_df) +
    facet_wrap(~tax_levels) +
    geom_bar(aes(x = dataset,
                 y = count_sums,
                 fill = dataset),
             stat = 'identity') +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = 'Data set',
         y = 'Number of unique labels',
         fill = 'Data set')

  return(list(plot_1, plot_2))
}


#'@title Display subtree
#'
#'@description Plot a tree with branches highlighted according to their
#'  membership in the provided data sets.
#'
#'@details This function takes in a base tree; an optional pruning specification
#'  as for [prune_tree()]; and optionally, one or two data sets used to specify
#'  subtrees of interest to highlight. If a pruning specification is given, the
#'  base tree will be pruned, and the pruned tree will be used as a new base
#'  tree. The function returns a tree diagram that plots the base tree (or the
#'  pruned base tree); if data sets were provided, branches of the plotted tree
#'  will be highlighted according to their membership in the provided data sets.
#'
#'  If one data set is provided, branches will be highlighted in one of two
#'  ways: whether they are in the data set or not.
#'
#'  If two data sets are provided, branches will be highlighted in one of four
#'  ways: whether they are in neither data set; whether they are in data set 1
#'  only; whether they are in data set 2 only; or whether they are in both
#'  datasets.
#'
#'  More information on how to specify some of the function arguments:
#'
#'  # How to specify pruning
#'
#'  ## `prune_to` See [prune_tree()] for details about how to specify `prune_to`
#'  and the list of additional arguments in `prune_args`. `prune_to` can be
#'  specified as a `data.frame`, as an integer vector of node numbers, as a
#'  character vector of node labels, as a character scalar defining a taxonomy
#'  level. Note that default pruning behavior changes depending on how you
#'  specify `prune_to`! This is explained in [prune_tree()].
#'
#'  ## `prune_args`
#'
#'  This is a named list of additional arguments for [prune_tree()]; please see
#'  the documentation for that function for details. The default is `NULL`, to
#'  use the defaults for [prune_tree()]. Names in this list can include
#'  `keep_descendants`, `adjust_branch_length`, and/or `tax_level_labels`.
#'
#'  # How to specify data sets
#'
#'  Arguments `data_1` and `data_2` are data sets that define subtrees of
#'  interest to highlight on the tree plot. They can be specified in one of
#'  several different ways. You are allowed to mix and match (specify `data_1`
#'  in a different way from `data_2`).
#'
#'  ## As `data.frame` objects containing classified entities
#'
#'  If `data_1` or `data_2` are provided as `data.frame` objects containing
#'  classified entities, they must be of the following format: There must be one
#'  variable corresponding to, and named for, each of the taxonomy levels
#'  specified in `tax_level_labels`, containing the label at that level for each
#'  entity. For example, to use the ChemOnt taxonomy, the `data.frame` must have
#'  columns `kingdom`, `superclass`, `class`, `subclass`, `level5`, `level6`,
#'  `level7`, `level8`, `level9`, `level10`, `level11`. Whenever an entity does
#'  not have a label at a given level, the corresponding element in the data
#'  table for that row and column should be `NA_character_`.
#'
#'  The `data.frame` must have at least one additional column that uniquely
#'  identifies individual entities; the name of that additional column does not
#'  matter, as long as it is not the same as one of the taxonomy levels.
#'
#'  For an example of a properly-formatted input `data.frame`, see the built-in
#'  dataset [biosolids_class]. The only required columns in [biosolids_class]
#'  are the ChemOnt taxonomy levels (`kingdom`, `superclass`, `class`,
#'  `subclass`, `level5`, `level6`, `level7`, `level8`, `level9`, `level10`,
#'  `level11`), and at least one entity identifer column, for example `DTXSID`.
#'
#'  ## As character vectors of node labels
#'
#'  If `data_1` or `data_2` are provided as character vectors of node labels,
#'  they refer to node labels in `base_tree` (or in the pruned base tree),
#'  *i.e.*, labels that appear in `base_tree$tip.label` or
#'  `base_tree$node.label`. Any labels in ` data_1` or `data_2` that do not
#'  appear in the base tree will be ignored.
#'
#'  The labels in `data_1` and/or `data_2`, plus all of their ancestors in the
#'  base tree, will be highlighted. Their descendants will not be highlighted
#'  unless those descendants themselves appear in `data_1` or `data_2`.
#'
#'  ## As vectors of node numbers
#'
#'  If `data_1` or `data_2` are provided as vectors of node numbers, they refer
#'  to node and tip index numbers in the base tree (or the pruned base tree).
#'  Tip numbers run from 1 to `ape::Ntip(base_tree)`, and node numbers run from
#'  `ape::Ntip(base_tree) + 1` to `ape::Ntip(base_tree) + ape::Nnode(base_tree)
#'  - 1`.
#'
#'  You might specify `data_1` or `data_2` as vectors of node numbers if, for
#'  example, you have used [phangorn::Ancestors()] or [phangorn::Descendants()]
#'  to trace ancestors or descendants of a specified node, and now you want to
#'  visualize those ancestors or descendants. [phangorn::Ancestors()] and
#'  [phangorn::Descendants()] both return lists or vectors of node numbers, so
#'  it would be convenient to be able to pass those vectors of node numbers
#'  directly to [display_subtree()].
#'
#'  The nodes/tips in `data_1` and/or `data_2`, plus all of their ancestors in
#'  the base tree, will be highlighted. Their descendants will not be
#'  highlighted unless those descendants themselves appear in `data_1` or
#'  `data_2`.
#'
#'  ## As tree objects
#'
#'  If `data_1` or `data_2` are provided as `phylo`-class tree objects, any
#'  labels in their \code{tip.label} or \code{node.label} elements that do not
#'  appear in the tip or node labels of the base tree will be ignored. In other
#'  words, only the portions of `data_1` and `data_2` that are subtrees of
#'  `base_tree` wil be used.
#'
#'  # How to specify `subtree_mapping`
#'
#'  Each element of `subtree_mapping` should be a vector of values. `color`
#'  should be a character vector of color names or hex specifications; `size`
#'  should be a numeric vector indicating line width;  `linetype` should be an
#'  integer vector indicating line type (1 = solid, 2 = dashed, 3 = dotted, 4 =
#'  dot-dashed, etc.).
#'
#'  If you supplied only `data_1`, then each vector should be size 2,
#'  corresponding (in order) to the following categories: "not in set 1" and "in
#'  set 1." If you want to name the elements of the vectors, they should be
#'  named according to the following recipe: `paste("Not in", name_1)`,
#'  `paste("In", name_1))`. If (and only if) you provide those names, then
#'  aesthetics will be assigned based on name rather than order.
#'
#'
#'  If you supplied both `data_1` and `data_2`, each vector should be 4 elements
#'  long, corresponding (in order) to the following categories: "neither set",
#'  "set 1 only", "set 2 only", and "both set 1 and set 2." If you want to name
#'  the elements of the vectors, they should be named according to the following
#'  recipe: `"Neither set"`, `paste(name_1, "only")`, `paste(name_2, "only")`,
#'  `paste(name_1, "and", name_2)`. If (and only if) you provide those names,
#'  then aesthetics will be assigned based on name rather than order.
#'
#'  If you provide a vector that is too long; it will be truncated to the needed
#'  length (with a warning). If you provide a vector that is too short, it will
#'  be recycled to the needed length (with a warning). If any of `color`,
#'  `size`, or `linetype` are not provided, the value from `base_opts` will be
#'  used for all data sets. If you provide any aesthetics (any named list
#'  elements) other than `color`/`colour`, `size`, or `linetype`, they will be
#'  ignored (with a warning).
#'
#'  # Helpful hints
#'
#'  If you supply `show_tiplabs = TRUE` and some of the tip labels are cut off,
#'  you can enlarge the margin using [ggplot2::xlim()]. For example:
#'
#' ```
#' my_tree <- display_subtree(prune_to = biosolids_class,
#' base_name = "Biosolids",
#' data_1 = usgs_class,
#' name_1 = "USGS",
#' show_tiplabs = TRUE)
#'
#' my_tree #labels are cut off
#'
#' my_tree + xlim(c(0, max(my_tree$data$x) + 50) #add extra space
#' ```
#'
#'@param base_tree The "base tree" to plot, as a `phylo`-class object.  Default
#'  is the full ChemOnt taxonomy tree, `chemont_tree`.
#'@param prune_to Optional: Use this argument if you want to prune the base tree
#'  before plotting. Default is `NULL`, which results in no pruning being done
#'  (*i.e.*, the `base_tree` is plotted as-is). Options are as for
#'  [prune_tree()]: What to *keep* from the base tree (everything else will be
#'  pruned away). May be a `data.frame` of classified data; a character vector
#'  of one or more labels in the tree (tip or internal node labels); an integer
#'  vector of one or more node numbers in the tree (tip or internal nodes); or a
#'  character scalar giving the name of a taxonomy level (one of the items in
#'  `tax_level_labels`).  See the help for [prune_tree()].
#'@param prune_args A named list of additional arguments to be passed to
#'  [prune_tree()] (used only if `prune_to` is non-`NULL`). Default:
#'  `list(adjust_branch_length = FALSE, tax_level_labels = chemont_tax_levels,
#'  keep_descendants = NULL)`. See the help for [prune_tree()] for details.
#'@param base_name Character: A string to use as the plot title. Usually this
#'  should name or describe the base tree (or the pruned base tree). Default
#'  `NULL` will result in no plot title.
#'@param data_1 Optional: Highlight branches of the base tree according to their
#'  membership in this data set. Default is `NULL`, to do no highlighting. If
#'  not `NULL`, must be one of the following options: A `data.frame` of
#'  classified entities; an integer vector of node numbers in the base tree; a
#'  character vector of node labels in the base tree; or a subtree of
#'  `base_tree` as a `phylo`-class object. `data_1` may also be a `list` of any
#'  of the preceding items, in which case all list elements must be of the same
#'  class, as they will be concatenated and the concatenation used as `data_1`.
#'  See Details.
#'@param data_2 Optional: Highlight branches of the base tree to compare
#'  membership in this data set and `data_1`. Default is `NULL`, to do no
#'  highlighting. If not `NULL`, one of the following options: A `data.frame` of
#'  classified entities; an integer vector of node numbers in the base tree; a
#'  character vector of node labels in the base tree; or a subtree of
#'  `base_tree` as a `phylo`-class object. `data_2` may also be a `list` of any
#'  of the preceding items, in which case all list elements must be of the same
#'  class, as they will be concatenated and the concatenation used as `data_2`.See Details.
#'@param name_1 Optional: Character string giving the name of `data_1` that
#'  should be used in the plot legend. Default is `"Set 1"`.
#'@param name_2 Optional: Character string giving the name of `data_2` (if
#'  `data_2` is provided) that should be used in the plot legend. Default is
#'  `"Set 2"`. If `data_2` is `NULL`, `name_2` will be ignored.
#'@param tax_level_labels Optional: a vector of the taxonomy level labels to be
#'  used. Default value:[chemont_tax_levels].
#'@param entity_id_col Character: the name of the variable in `prune_to`,
#'  `data_1`, and/or `data_2` that identifies entities, if those arguments are
#'  provided as `data.frame`s of classified entities. Default `NULL` to assume
#'  that each row is a unique entity.
#'@param layout [ggtree::ggtree()] tree layout option. Default `"circular."`
#'@param base_opts List of parameters to control color, size, and linetype of
#'  the base tree diagram. Default `list(color = "black", size = 0.5, linetype =
#'  1)`. Should be a named list with one or more of `color`, `size`, and
#'  `linetype`. Each should be a single value: `color` should be a character
#'  specifying a color name or hex code; `size` should be numeric specifying
#'  line weight; and `linetype` should be an integer specifying line type (1 =
#'  solid, 2 = dotted, 3 = dashed, etc.) If any of the three are not provided,
#'  their value from the default list will be substituted. For example, if you
#'  provide `base_opts = list(color = "red")`, it will be the same as if you
#'  provided `base_opts = list(color = "red", size = 0.5, linetype = 1)`.
#'@param bg_tree_scale Numeric: The factor by which to increase or decrease the
#'  size (linewidth) of the base tree (as specified in `base_opts$size`), on
#'  which the highlighted tree will be overlaid. Default `1.8`. In short, a tree
#'  is first drawn using the `base_opts`, but with a slightly thicker line
#'  weight. Then, the highlighted tree is drawn on top of it. This is a way of
#'  drawing an outline or border around the highlighted tree, using the base
#'  color. The larger `bg_scale` is, the wider the border will be.
#'@param highlight_by One of `"set"` (default), `"overlap"`, `"simil"`, or
#'  `"none"`. `"set"` highlights branches by their membership in `data_1` and/or
#'  `data_2` if provided. `"overlap"` highlights branches according to
#'  fractional overlap in number of entities, if two data sets were provided.
#'  `"sim"` highlights branches according to their maximum similarity of
#'  taxonomic ancestry between two data sets, if provided. `"none"` does not
#'  highlight branches. To use `highlight_by = "set"`, `data_1` must be provided
#'  (and `data_2` may be provided). To use `"overlap"`, either `prune_to` and
#'  `data_1` must both be provided as `data.frame`s of classified entities (in
#'  which case the overlap in number of entities will be computed between
#'  `prune_to` and `data_1`), or `data_1` and `data_2` must both be provided as
#'  `data.frame`s of classified entities (in which case the overlap in number of
#'  entities will be computed between `data_1` and `data_2`). If all three of
#'  `prune_to`, `data_1`, and `data_2` are provided as `data.frame`s of
#'  classified entities, then the overlap will be computed between `data_1` and
#'  `data_2`. To use `"simil"`, either `prune_to` and `data_1` must be provided,
#'  or `data_1` and `data_2` must be provided, although they may be provided in
#'  any of the allowable formats (unlike `highlight_by = "overlap"`, they need
#'  not be provided as `data.frame`s of classified entities).
#'@param subtree_mapping List of aesthetics used to control branch highlighting.
#'  Must be a named list with one or more of the following named elements:
#'  `color` (or `colour`), `size`, and/or `linetype`. Default `NULL` to choose a
#'  default color mapping based on `highlight_by`. If `highlight_by %in%
#'  c("overlap", "simil")`, then elements `size` or `linetype` in
#'  `subtree_mapping` will be ignored.
#'@param show_tiplabs `TRUE`/`FALSE`: Whether to display the tip labels. Default
#'  `TRUE`. Note that when `show_tiplabs = TRUE`, clade labels will be
#'  suppressed (*i.e.*, arguments `clade_opts` and `clade_level` will be
#'  ignored).
#'@param tiplab_opts A named list of additional arguments to
#'  [ggtree::geom_tiplab()]; see the documentation for that function for more
#'  options. Default `NULL` to use defaults for that function. For example, you
#'  could specify `tiplab_opts = list(offset = 3)` to increase the space between
#'  tips and labels. If `size` is not specified, this function will attempt to
#'  auto-determine size based on the number of tips (`size = 3` if fewer than
#'  200 tips; `size = 1.5` if between 200 and 500 tips; if more than 500 tips,
#'  `size = 0.5`).
#'@param show_tippoints `TRUE`/`FALSE`: Whether to display tip points. Default
#'  `FALSE`.
#'@param tippoint_opts A named list of additional arguments to
#'  [ggtree::geom_tippoint()]; see the documentation for that function for more
#'  options. Default `NULL` to use defaults for that function. If `size` is not
#'  specified, this function will attempt to auto-determine size based on the
#'  number of tips (`size = 1.5` if fewer than 200 tips; `size = 1` if between
#'  200 and 500 tips; if more than 500 tips, `size = 0.5`).
#'@param show_nodepoints `TRUE`/`FALSE`: Whether to display (internal) node
#'  points. Default `FALSE`.
#'@param nodepoint_opts A named list of additional arguments to
#'  [ggtree::geom_nodepoint()]; see the documentation for that function for more
#'  options. Default `NULL` to use defaults for that function. If `size` is not
#'  specified, this function will attempt to auto-determine size based on the
#'  number of tips (`size = 1.5` if fewer than 200 tips; `size = 1` if between
#'  200 and 500 tips; if more than 500 tips, `size = 0.5`).
#'@param sim_mat Optional: A numeric matrix of pairwise similarities including
#'  all the nodes in the base tree (or just the pruned tree, if `prune_to` is
#'  non-`NULL`), as produced by [similarity_matrix()]. Default
#'  `chemont_jaccard`. Used only if `highlight_by = "sim"`. If `NULL`, then
#'  similarity highlighting will not be performed.
#'@param clade_level Integer or character: The taxonomy level at which to draw
#'  clade labels, if any. Default is `NULL`, which will suppress clade labels
#'  altogether. `clade_level = "auto"` is a special value that automatically
#'  selects the taxonomy level of the MRCA of the tips plus one, or level 2,
#'  whichever is the greater number. See [add_cladelabels()]. If `show_tiplabs =
#'  TRUE`, `clade_level` will be ignored and clade labels will not be drawn.
#'@param clade_opts List of parameters defining aesthetic options for labeling
#'  clades. See [add_cladelabels()] for options.  Default is `list(wrap = 20,
#'  barsize = "alternate", fontsize = 3, lineheight = 0.7, default_to_tip =
#'  TRUE, draw_text = TRUE)`. You can provide any of the additional parameters
#'  of [ggtree::geom_cladelab()]. If `show_tiplabs = TRUE`, `clade_opts` will be
#'  ignored and clade labels will not be drawn.
#'@return A [ggtree::ggtree()] object visualizing the full base tree, with
#'  branches highlighted to indicate presence in `data_1`, in `data_2` (if
#'  supplied), neither, or both.
#'@author Caroline Ring, Paul Kruse
#' @examples
#'
#' #make a smaller base tree for visibility
#' oh_tree <- prune_tree(tree = chemont_tree,
#' prune_to = "Organohalogen compounds")
#'
#' #show this tree by itself
#' display_subtree(base_tree = oh_tree)
#'
#' #one data set
#' display_subtree(base_tree = oh_tree,
#'  data_1 = biosolids_class,
#'  entity_id_col = "DTXSID",
#'  highlight_by = "set")
#'
#' #two data sets
#' display_subtree(base_tree = oh_tree,
#'  data_1 = biosolids_class,
#' data_2 = usgs_class,
#' entity_id_col = "DTXSID",
#'  highlight_by = "set")
#'
#' #increasing line width
#' display_subtree(base_tree = oh_tree,
#'  data_1 = biosolids_class,
#' data_2 = usgs_class,
#' entity_id_col = "DTXSID",
#'  base_opts = list(size = 1),
#'  highlight_by = "set")
#'
#' #different subtree color mapping
#' display_subtree(base_tree = oh_tree,
#'  data_1 = biosolids_class,
#' data_2 = usgs_class,
#' entity_id_col = "DTXSID",
#' base_opts = list(size = 1),
#'  subtree_mapping = list(color = c("black", "red", "blue", "purple")),
#'  highlight_by = "set")
#'
#'  #clade labels
#'  display_subtree(base_tree = oh_tree,
#'   data_1 = biosolids_class,
#' data_2 = usgs_class,
#' entity_id_col = "DTXSID",
#'  base_opts = list(size = 1),
#'  show_tiplabs = FALSE, #they'll clash with clade labels
#' clade_level = 2,
#'  highlight_by = "set")
#'
#' #prune full ChemOnt tree to the union of BIOSOLIDS2021 and USGSWATER classes only
#' #and color by set membership
#' display_subtree(base_tree = chemont_tree,
#' prune_to = list(biosolids_class, usgs_class),
#' data_1 = usgs_class,
#' name_1 = "USGS Water",
#' data_2 = biosolids_class,
#' name_2 = "Biosolids",
#' entity_id_col = "DTXSID",
#' base_opts = list(size = 1),
#' bg_tree_scale = 0, #suppress border drawing
#' highlight_by = "set",
#' show_tiplabs = FALSE,
#' show_tippoints = TRUE,
#' tippoint_opts = list(size = 2),
#' clade_level = 2,
#' clade_opts = list(offset = 3,
#' offset.text = 1,
#' fontsize = 4))
#'
#' #prune to the union of BIOSOLIDS2021 and USGSWATER classes only
#' #and color by Jaccard similarity of classifications between BIOSOLIDS2021 and USGS Water.
#' #compare this plot to the previous one!
#' display_subtree(base_tree = chemont_tree,
#' prune_to = list(biosolids_class, usgs_class),
#' data_1 = usgs_class,
#' name_1 = "USGS Water",
#' data_2 = biosolids_class,
#' name_2 = "Biosolids",
#' entity_id_col = "DTXSID",
#' base_opts = list(size = 1),
#' bg_tree_scale = 0, #suppress border drawing
#' highlight_by = "sim",
#' sim_mat = chemont_jaccard,
#' show_tiplabs = FALSE,
#' show_tippoints = TRUE,
#' tippoint_opts = list(size = 2),
#' clade_level = 2,
#' clade_opts = list(offset = 3,
#' offset.text = 1,
#' fontsize = 4))
#'
#'@export
#'@import ggtree
#'@importFrom magrittr `%>%`
display_subtree <- function(base_tree = chemont_tree,
                            prune_to = NULL,
                            prune_args = list(adjust_branch_length = FALSE,
                                              tax_level_labels = chemont_tax_levels,
                                              keep_descendants = NULL),
                            base_name = NULL,
                            data_1 = NULL,
                            data_2 = NULL,
                            name_1 = "Set 1",
                            name_2 = "Set 2",
                            tax_level_labels = chemont_tax_levels,
                            entity_id_col = NULL,
                            layout = "circular",
                            base_opts = list("color" = "black",
                                             "size" = 0.5,
                                             "linetype" = 1),
                            bg_tree_scale = 1.8,
                            highlight_by = "set", #or "overlap" or "sim"
                            subtree_mapping = NULL,
                            show_tiplabs = TRUE,
                            tiplab_opts = NULL,
                            show_tippoints = FALSE,
                            tippoint_opts = NULL,
                            show_nodepoints = FALSE,
                            nodepoint_opts = NULL,
                            sim_mat = chemont_jaccard,
                            clade_level = NULL,
                            clade_opts = list(wrap = 20,
                                              barsize = "alternate",
                                              fontsize = 3,
                                              lineheight = 0.7,
                                              default_to_tip = TRUE,
                                              draw_text = TRUE)
){


  ############
  #Base options
  #############
  #Keep defaults for any base_opts not otherwise specified
  base_opts_default <- list("color" = "black",
                            "size" = 0.5,
                            "linetype" = 1)
  #harmonize color/colour
  if("colour" %in% names(base_opts)){
    base_opts$color <- base_opts$colour
    base_opts$colour <- NULL
  }

  base_opts <- c(base_opts,
                 base_opts_default[setdiff(names(base_opts_default),
                                           names(base_opts))])

  ############
  #Base tree
  #############

  # Prune base tree if user specified pruning
  if(!is.null(prune_to)){
    #concatenate prune_to here, since it will be used again later.
    prune_to <- prune_c(prune_to = prune_to)
    base_tree <- do.call(prune_tree,
                         args = c(list(tree = base_tree,
                                       prune_to = prune_to,
                                       entity_id_col = entity_id_col,
                                       tax_level_labels = tax_level_labels),
                                  prune_args[setdiff(names(prune_args),
                                                     c("tree",
                                                       "prune_to",
                                                       "entity_id_col",
                                                       "tax_level_labels")
                                                     )
                                             ]
                                  )
    )
  }

  ###############
  # Tip label, tip point, and node point options
  #################
  #(auto-size depends on the pruned tree, so it must be after pruning happens)
  if(is.null(tiplab_opts$size)){
    if (length(base_tree$tip.label) <= 200){
      tiplab_size = 3
    } else if (length(base_tree$tip.label) <= 500){
      tiplab_size = 1.5
    } else {
      tiplab_size = 0.5
    }

    tiplab_opts <- c(list(size = tiplab_size),
                     tiplab_opts)
  }

  if(is.null(tippoint_opts$size)){
    if (length(base_tree$tip.label) <= 200){
      tippoint_size = 1.5
    } else if (length(base_tree$tip.label) <= 500){
      tippoint_size = 1
    } else {
      tippoint_size = 0.5
    }

    tippoint_opts <- c(list(size = tippoint_size),
                       tippoint_opts)
  }

  if(is.null(nodepoint_opts$size)){
    if (length(base_tree$tip.label) <= 200){
      nodepoint_size = 1.5
    } else if (length(base_tree$tip.label) <= 500){
      nodepoint_size = 1
    } else {
      nodepoint_size = 0.5
    }

    nodepoint_opts <- c(list(size = nodepoint_size),
                        nodepoint_opts)
  }

  ########################################
  # Handle different formats for data_1
  ########################################
  if(!is.null(data_1)){

    #If it is a list, look at the data type of the list elements,
    #and concatenate them as appropriate
    if(is(data_1, "list")){ #this will appropriately return FALSE if it is a data.frame
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
        stop(paste("data_1 is a list, but one or more elements are not one of",
                   "the recognized classes: data.frame, numeric, character, or phylo."))
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
        if(is.null(entity_id_col)){
          message(paste("data_1 is a data.frame,",
                        "but no entity ID variable name has been specified.",
                        "Each row will be assumed to be one entity."))
          #add a new variable to the data
          #ensure it does not conflict with any of the existing variable names
          entity_id_col <- rev(
            make.names(
              names = c(names(data_1),
                        "id"
              ),
              unique = TRUE
            )
          )[1]
          data_1[[entity_id_col]] <- 1:nrow(data_1)
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
      stop(paste("data_1 is not one of the recognized classes:",
                 "data.frame, numeric, character, or phylo."))
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

    #convert ancestors to labels
    data_1_all_labs <- get_label_from_node(node = data_1_all,
                                           tree = base_tree)
  } #end if(!is.null(data_1))

  ########################################
  # Handle different formats for data_2
  ########################################
  if (!is.null(data_2)) {

    #If it is a list, look at the data type of the list elements,
    #and concatenate them as appropriate
    if(is(data_2, "list")){ #this will return FALSE if it is a data.frame
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
        stop(paste0("data_2 is a list, but one or more elements are not one of",
                    "the recognized classes: data.frame, numeric, character, or phylo."))
      }
      #and keep only the unique combined elements
      data_2 <- unique(data_2)
    } #end if(is.list(data_2))

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
        #add a new variable to the data
        #ensure it does not conflict with any of the existing variable names
        if(is.null(entity_id_col)){
          warning(paste("data_2 is a data.frame,",
                        "but no entity ID variable name has been specified.",
                        "Each row will be assumed to be one entity."))
          #add a new variable to the data
          #ensure it does not conflict with any of the existing variable names
          entity_id_col <- rev(
            make.names(
              names = c(names(data_2),
                        "id"
              ),
              unique = TRUE
            )
          )[1]
          data_2[[entity_id_col]] <- 1:nrow(data_2)
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
      stop(paste("data_2 is not one of the recognized classes: data.frame,",
                 "numeric, character, or phylo."))
    }

    #get ancestors of these nodes
    data_2_all <- unique(
      c(
        data_2_nodes,
        unlist(phangorn::Ancestors(x = base_tree,
                                   node = data_2_nodes)
        )
      )
    )

    #convert all ancestors to labels
    data_2_all_labs <- get_label_from_node(node = data_2_all,
                                           tree = base_tree)
  } #end if(!is.null(data_2))


  #############################
  # Set Up Branch Coloring Data
  ##############################

  #Data frame with all node numbers in taxonomy tree
  cohort_data <- get_tree_df(tree = base_tree)

  if(highlight_by %in% "set"){
    #Categorical column: is each node in Data Set 1?
    #0 = no, 1 = yes

    if(!is.null(data_1)){
      cohort_data$inSet1 <- ifelse(cohort_data$node %in% data_1_all,
                                   1L,
                                   0L)

      if(!is.null(data_2)){
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
                                                       "Both sets"))

      }else{
        #Categorical column: Is each node in Data Set 1 or not?
        cohort_data$list_presence <- factor(cohort_data$inSet1,
                                            levels = 0:1,
                                            labels = c(paste("Not in", name_1),
                                                       paste("In", name_1)))
      }
    } #end if(!is.null(data_1))


    #########################################################
    # Set up aesthetics for set presence coloring
    #########################################################

    if(is.null(data_1)){
      #ignore any subtree mapping
      subtree_mapping <- NULL
        highlight_by <- "none"
        tree_plot <- do.call(ggtree,
                             c(list(tr = base_tree,
                                    layout = layout),
                               base_opts
                             )
        )
    }else{ #if data_1 supplied, check subtree mapping

      if(is.null(subtree_mapping)){
        if(!is.null(data_2)){
          #If no aesthetic mapping for list presence specified,
          #then default to color only
          subtree_mapping <- list(color = c("gray70",
                                            "#66C2A5",
                                            "#8DA0CB",
                                            "#FC8D62"))
        }else{
          subtree_mapping <- list(color = c("gray70",
                                            "#66C2A5"))
        }
      }



      #Check subtree mapping names
      #They need to be valid aesthetics for ggtree
      good_subtree_map <- any(c("color", "colour", "size", "linetype") %in%
                                names(subtree_mapping))

      if(!is.list(subtree_mapping) |
         is.null(names(subtree_mapping)) |
         !isTRUE(good_subtree_map)){
        stop(paste("subtree_mapping should be a list with one or more named elements.",
                   "Names must be one or more of 'color' (or 'colour'),",
                   "'size', and/or 'linetype'. "))
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

        #harmonize color/colour
        if("colour" %in% names(subtree_mapping)){
          subtree_mapping$color <- subtree_mapping$colour
          subtree_mapping$colour <- NULL
        }

      } #end if(!is.list(subtree_mapping) |
      # is.null(names(subtree_mapping)) |
      #   !isTRUE(bad_subtree_map))


      ####################################
      # PLOT BACKGROUND TREE FOR SET HIGHLIGHTING
      #######################################
      #Use base options, unless they will be mapped to list presence later
      #(aes() doesn't seem to overwrite them as expected)
      #e.g. if subtree_mapping has a "color" element, don't use base_opts$color
      # tree_plot <- do.call(ggtree,
      #                      c(list(tr = base_tree,
      #                             layout = layout),
      #                        base_opts[setdiff(names(base_opts),
      #                                          names(subtree_mapping))
      #                        ]
      #                      )
      # )

      bg_opts <- base_opts
      bg_opts$size <- base_opts$size * bg_tree_scale
      tree_plot <- do.call(ggtree,
                           c(list(tr = base_tree,
                                  layout = layout),
                             bg_opts
                           )
      )

      #add background tip/node points if requested
      if(show_tippoints %in% TRUE){
        #scale size
        bg_tippoint_opts <- tippoint_opts
        bg_tippoint_opts$size <- bg_tippoint_opts$size * bg_tree_scale
        tree_plot <- tree_plot +
          do.call(geom_tippoint,
                  args = bg_tippoint_opts)
      }

      if(show_nodepoints %in% TRUE){
        bg_nodepoint_opts <- nodepoint_opts
        bg_nodepoint_opts$size <- bg_nodepoint_opts$size * bg_tree_scale
        tree_plot <- tree_plot +
          do.call(geom_nodepoint,
                  args = bg_nodepoint_opts
          )
      }


      #Name the subtree_mapping items after the categories in cohort_data$list_presence
      subtree_mapping <- sapply(subtree_mapping,
                                function(x) setNames(x,
                                                     levels(
                                                       cohort_data$list_presence
                                                     )
                                ),
                                simplify = FALSE,
                                USE.NAMES = TRUE)
      #set up for aes call -- list of aesthetic mappings,
      #all applied to "list_presence" in cohort_data
      subtree_aes <- replicate(n= length(subtree_mapping),
                               expr = quote(list_presence))
      #name the list elements after the aesthetics in subtree_mapping
      subtree_aes <- setNames(subtree_aes, names(subtree_mapping))
      #you end up with something like this:
      #`subtree_aes <- list(color = quote(list_presence),
      #                   size = quote(list_presence))``
      #`do.call(aes, subtree_aes)` is then equivalent to:
      #`aes(color = list_presence, size = list_presence)`

      #Prepare a list of manual scales as provided in subtree_mapping
      scale_list <- sapply(names(subtree_mapping),
                           function(aesthetic) {
                             ggplot2::scale_discrete_manual(
                               aesthetics = aesthetic,
                               name = "List presence",
                               values =  subtree_mapping[[aesthetic]],
                               breaks = levels(cohort_data$list_presence),
                               limits = levels(cohort_data$list_presence)
                             )
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE
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

      #any aesthetics not scaled to data, set to their base values.
      unscaled_aesthetics <- setdiff(names(base_opts),
                                     names(scale_list))

      #add list presence highlighting to tree plot
      tree_plot <- tree_plot %<+% cohort_data +
        do.call(geom_tree,
                args = c(
                  list(do.call(aes,
                               subtree_aes)),
                  base_opts[unscaled_aesthetics])) +
        scale_list
    } #end if(!is.null(data_1))
    ### end if highlight_by %in% "set"#####
  }else if(highlight_by %in% c("overlap", "sim")){
    ##############################
    # HIGHLIGHTING BY OVERLAP OR SIMILARITY
    ##############################

    #if neither data_1 nor data_2 were provided, then ignore with a warning
    if(is.null(data_1)){
      message(paste("To use highlight_by = ", paste0(highlight_by, ","),
                    "either data_1 and data_2 must be provided",
                    "or prune_to and data_1 must be provided.",
                    "Here, data_1 was not provided.",
                    "Therefore, highlight_by will be ignored,",
                    "and only the base tree will be plotted."))
      dat_A <- NULL
    }else{ #if data_1 is provided
      if(is.null(data_2)){
        #if data_2 not provided, then prune_to must be provided
        if(is.null(prune_to)){
          #if data_2 not provided and prune_to also not provided,
          #there is nothing to calculate overlap with
          message(paste("To use highlight_by = ", paste0(highlight_by, ","),
                        "either data_1 and data_2 must be provided",
                        "or prune_to and data_1 must be provided.",
                        "Here, only data_1 was provided (neither prune_to nor data_2 provided).",
                        "Therefore, highlight_by will be ignored,",
                        "and only the base tree will be plotted."))
          dat_A <- NULL
        }else{
          dat_A <- prune_to
          dat_B <- data_1
          name_A <- base_name
          name_B <- name_1
          dat_A_labs <- cohort_data$Name
          dat_B_labs <- data_1_all_labs
          dat_A_from <- "prune_to"
          dat_B_from < "data_1"
        }
      }else{
        dat_A <- data_1
        dat_B <- data_2
        name_A <- name_1
        name_B <- name_2
        dat_A_labs <- data_1_all_labs
        dat_B_labs <- data_2_all_labs
        dat_A_from <- "data_1"
        dat_B_from < "data_2"
      }
    }

    if(!is.null(dat_A)){ #i.e., if we can proceed
      ##################
      # Construct cohort_data for overlap or similarity
      ################
      if(highlight_by %in% "overlap"){
        #ensure that dat_A and dat_B are both data.frames of classified entities
        if(!is.data.frame(dat_A)){
          stop(paste("Error in treecompareR::display_subtree():",
                     "to use highlight_by = 'overlap'",
                     dat_A_from,
                     "must be provided as a data.frame of classified entities."))
        }else{
          if(!all(c(tax_level_labels,
                    entity_id_col) %in%
                  names(dat_A))){
            stop(paste("Error in treecompareR::display_subtree():",
                       "to use highlight_by = 'overlap'",
                       dat_A_from,
                       "must be provided as a data.frame of classified entities",
                       "whose variable names must include tax_level_labels = ",
                       paste(tax_level_labels, collapse = ", "),
                       "and entity_id_col = ",
                       paste(entity_id_col, collapse = ", ")))
          }
        }

        if(!is.data.frame(dat_B)){
          stop(paste("Error in treecompareR::display_subtree():",
                     "to use highlight_by = 'overlap'",
                     dat_B_from,
                     "must be provided as a data.frame of classified entities."))
        }else{
          if(!all(c(tax_level_labels,
                    entity_id_col) %in%
                  names(dat_B))){
            stop(paste("Error in treecompareR::display_subtree():",
                       "to use highlight_by = 'overlap'",
                       dat_B_from,
                       "must be provided as a data.frame of classified entities",
                       "whose variable names must include tax_level_labels = ",
                       paste(tax_level_labels, collapse = ", "),
                       "and entity_id_col = ",
                       paste(entity_id_col, collapse = ", ")))
          }
        }

        message(
          paste(
            "computing overlap in numbers of entities for each label of",
            name_A,
            "and",
            name_B)
        )
        overlap_dat <- lapply(tax_level_labels,
                              function(this_level) {
                                calc_number_overlap(data_1 = dat_A,
                                                    data_2 = dat_B,
                                                    entity_id_col = entity_id_col,
                                                    at_level = this_level,
                                                    tax_level_labels = tax_level_labels) %>%
                                  dplyr::rename(
                                    label = dplyr::all_of(this_level))
                              }
        )
        #calculate fractional overlaps of entities at each node
        cohort_data <- overlap_dat  %>%
          dplyr::bind_rows() %>%
          dplyr::mutate(node = get_node_from_label(
            label,
            base_tree))

      }else if(highlight_by %in% "sim"){
        if(is.null(sim_mat)){
          stop("Pre-computed similarity matrix required")
        }else{
          #check that sim_mat is a numeric matrix
          if(!is.matrix(sim_mat)){
            stop("sim_mat must be a numeric matrix")
          }else{
            if(!is.numeric(sim_mat)){
              #coerce to numeric
              sim_mat <- apply(sim_mat,
                               2,
                               as.numeric)
            }
          }

          #check that sim_mat has row/column names that match labels
          tree_labels <- c(base_tree$tip.label,
                           base_tree$node.label)
          if(!(any(tree_labels %in% rownames(sim_mat)))){
            stop("sim_mat needs row and column names that match the tree tip and node labels")
          }

          if(!(any(tree_labels %in% colnames(sim_mat)))){
            stop("sim_mat needs row and column names that match the tree tip and node labels")
          }

          message(
            paste(
              "looking up similarity of taxonomic ancestry from",
              "provided sim_mat for each pair of labels in",
              name_A,
              "and",
              name_B))
          #similarity matrix usually excludes the root node
          #add it back in
          root_label <- base_tree$node.label[1]

          if(root_label %in% c(dat_A_labs,
                               dat_B_labs)){
            if(!(root_label %in% colnames(sim_mat))){
              sim_mat1 <- cbind(sim_mat,
                                rep(NA_real_, nrow(sim_mat))
              )
              sim_mat1 <- rbind(sim_mat1,
                                c(rep(NA_real_, ncol(sim_mat1)-1),
                                  1.0))
            }

            rownames(sim_mat1) <- c(rownames(sim_mat),
                                    root_label)
            colnames(sim_mat1) <- c(colnames(sim_mat),
                                    root_label)
          }

          #use match() here so that labels not in the dimnames of sim_mat return NA,
          #rather than stopping with an error about subscript out of bounds.
          row_inds <- match(dat_A_labs, rownames(sim_mat))
          col_inds <- match(dat_B_labs, colnames(sim_mat))
          sim_mat <- sim_mat[row_inds,
                             col_inds]
        }

        #Each node in the tree is assigned its maximum similarity between the datasets,
        #unless it does not appear in either data set,
        #in which case it is assigned NA.
        #this means similarity will be 1 if a node appears in both datasets.
        cohort_data <- cohort_data %>%
          dplyr::mutate(simil_col =   {
            col_inds <- match(Name, colnames(sim_mat))
            apply(sim_mat[, col_inds],
                  MARGIN = 2,
                  #have to use pmax to return NA without throwing warnings
                  FUN = function(x) do.call(pmax,
                                            args = c(as.list(x),
                                                     list(na.rm = TRUE))))
          }
          ) %>%
          dplyr::mutate(simil_row = {
            row_inds <- match(Name, rownames(sim_mat))
            apply(sim_mat[row_inds, ],
                  MARGIN = 1,
                  #have to use pmax to return NA without throwing warnings
                  FUN = function(x) do.call(pmax,
                                            args = c(as.list(x),
                                                     list(na.rm = TRUE))))
          }
          ) %>%
          dplyr::mutate(simil = pmax(simil_col,
                                     simil_row,
                                     na.rm = TRUE))
      }

      #####################
      # Plot background tree for overlap or similarity
      #######################
      #plot a "background" tree with size a little bigger
      #this will create a "border" around the colored branches
      #useful when color is light and plot has a white background, for example
      bg_opts <- base_opts
      bg_opts$size <- base_opts$size * bg_tree_scale
      tree_plot <- do.call(ggtree,
                           c(list(tr = base_tree,
                                  layout = layout),
                             bg_opts
                           )
      )

      if(show_tippoints %in% TRUE){
        tree_plot <- tree_plot +
          geom_tippoint(size = tippoint_opts$size * bg_tree_scale)
      }

      if(show_nodepoints %in% TRUE){
        tree_plot <- tree_plot +
          geom_nodepoint(size = nodepoint_opts$size * bg_tree_scale)
      }


      #####################
      # Set up color mapping for overlap/similarity highlighting
      ######################
      if(is.null(subtree_mapping)){
        #use a default color scale
        subtree_mapping <- list(
          #this is the result of viridis::viridis(n=20)
          color = c('#440154FF',
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
        )
      }else{
        #Check subtree mapping names
        #It needs to contain "color" or "colour"
        good_subtree_map <- any(c("color", "colour") %in%
                                  names(subtree_mapping))

        if(!is.list(subtree_mapping) |
           is.null(names(subtree_mapping)) |
           !isTRUE(good_subtree_map)){
          stop(paste("subtree_mapping should be a list with at least one element",
                     "named 'color' (or 'colour')."))
        }else{
          #if something reasonable was provided for subtree_map,
          #check if any other aesthetics were provided and will be ignored
          bad_aes <- setdiff(names(subtree_mapping),
                             c("color", "colour"))
          good_aes <- intersect(names(subtree_mapping),
                                c("color", "colour"))
          if(length(bad_aes)>0){
            message(paste("In subtree_mapping, only aesthetics",
                          paste(good_aes, collapse = ", "),
                          "will be used. Aesthetics",
                          paste(bad_aes, collapse = ", "),
                          "were provided, but will be ignored."))
          }

          if("colour" %in% names(subtree_mapping)){
            subtree_mapping$color <- subtree_mapping$colour
            subtree_mapping$colour <- NULL
          }
        } #end if(!is.list(subtree_mapping) |
        # is.null(names(subtree_mapping)) |
        #   !isTRUE(bad_subtree_map))
      } #end if !is.null(subtree_mapping)


      ##########
      # Overlay highlighted tree for overlap or similarity
      ############

      #Set up title for color legend
      color_title <- ifelse(highlight_by %in% "overlap",
                            paste0("Entity overlap % btw\n",
                                   name_A,
                                   "\nand\n",
                                   name_B),
                            paste0("Max similarity btw\n",
                                   name_A,
                                   "\nand\n",
                                   name_B)
      )

      tree_plot <- tree_plot %<+% cohort_data + #add the similarity data
        geom_tree(size = base_opts$size,
                  aes(color = simil)) + #color branches by similarity
        scale_color_gradientn(colors = subtree_mapping$color,
                              limits = c(0,1),
                              guide = "colorbar",
                              name = color_title)

    }else{ #if(is.null(dat_A)), just plot the base tree
      highlight_by <- "none"
      tree_plot <- do.call(ggtree,
                           c(list(tr = base_tree,
                                  layout = layout),
                             base_opts
                           )
      )
    }
  }else if(highlight_by %in% "none"){
    #draw the base tree only
    tree_plot <- do.call(ggtree,
                         c(list(tr = base_tree,
                                layout = layout),
                           base_opts
                         )
    )
  }else{
    message(paste("highlight_by should be one of 'set', 'overlap', 'sim', or 'none'.",
                  paste0("highlight_by = ", highlight_by, "."),
                  "Treating it as 'none'."))
    tree_plot <- do.call(ggtree,
                         c(list(tr = base_tree,
                                layout = layout),
                           base_opts
                         )
    )
  } #end if/else statements for different highlight_by options


  ######
  # Add tip labels (if requested)
  ######
  if (show_tiplabs %in% TRUE){
    if(highlight_by %in% c("set")){
      #then use the specified subtree_aes
      tiplab_args <- c(list(do.call(aes, subtree_aes)),
                       tiplab_opts)
    }else if(highlight_by %in% c("overlap", "sim")){
      tiplab_args <- c(list(aes(color = simil)),
                       tiplab_opts)
    }else{
      #just use base options
      tiplab_args <- tiplab_opts
    }
    tree_plot <- tree_plot +
      do.call(ggtree::geom_tiplab,
              args = tiplab_args)
  }

  ######
  # Add tip points (if requested)
  ######
  if (show_tippoints %in% TRUE){
    if(highlight_by %in% c("set")){
      #then use the specified subtree_aes
      tippoint_args <- c(list(do.call(aes, subtree_aes)),
                         tippoint_opts)
    }else if(highlight_by %in% c("overlap", "sim")){
      tippoint_args <- c(list(aes(color = simil)),
                         tippoint_opts)
    }else{
      #just use base options
      tippoint_args <- tippoint_opts
    }
    tree_plot <- tree_plot +
      do.call(ggtree::geom_tippoint,
              args = tippoint_args)
  }

  ######
  # Add node points (if requested)
  ######
  if (show_nodepoints %in% TRUE){
    if(highlight_by %in% c("set")){
      #then use the specified subtree_aes
      nodepoint_args <- c(list(do.call(aes, subtree_aes)),
                          nodepoint_opts)
    }else if(highlight_by %in% c("overlap", "sim")){
      nodepoint_args <- c(list(aes(color = simil)),
                          nodepoint_opts)
    }else{
      #just use base options
      nodepoint_args <- nodepoint_opts
    }

    tree_plot <- tree_plot +
      do.call(ggtree::geom_nodepoint,
              args = nodepoint_args)
  }

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

  tree_plot <- tree_plot +
    theme(legend.position = "left")

  return(tree_plot)

}


#'@title Display tree with overlap
#'
#'@description Display a tree with branches highlighted according to membership
#'  in two classified data sets, annotated with numbers of entities and
#'  fractional overlap between entities in the two data sets.
#'
#'@details Given a base tree and two classified data sets: Display a tree with
#'  branches highlighted according to membership in neither set, `data_1` only,
#'  `data_2` only, or both sets (as for [display_subtree()] with both `data_1`
#'  and `data_2` provided). Add a heatmap annotation around the tips of the
#'  tree, indicating the number of entities at each tip in each dataset, and the
#'  fractional overlap in entities at each tip between datasets. Fractional
#'  overlap is calculated using [calc_number_overlap()].
#'
#'
#'@param base_tree A `phylo`-class tree object to plot (see [ape::read.tree()]
#'  for details on this class). Default [chemont_tree].
#'@param base_name A name for the base tree (used for the plot title). Default
#'  `NULL` for no title.
#'@param data_1 A `data.frame` containing a classified list of entities
#'@param name_1 A name for `data_1`, used for the plot legend. Default `'Set1'`.
#'@param data_2 Another `data.frame` containing a classified list of entities
#'@param name_2 A name for `data_2`, used for the plot legend. Default `'Set2'`.
#'@param entity_id_col Character: the name of the variable in `data_1` and
#'  `data_2` that identifies entities. Must be the same in both data sets.
#'  Default `NULL` to assume that each row is a unique entity.
#'@param group_level Taxonomy level at which to aggregate entities. Default
#'  `"terminal"`: Calculate number of entities and overlap for each tip label.
#'  Can also be any element of `tax_level_labels` (i.e., a string), or an
#'  integer between 1 and `length(tax_level_labels)`. In this case, entities
#'  will be grouped by unique labels at the specified taxonomic level, rather
#'  than by terminal (tip) labels, for calculation and plotting of number of
#'  chemicals and overlap.
#'@param tax_level_labels Taxonomy level labels. Default [chemont_tax_levels] to
#'  use ChemOnt taxonomy levels.
#'@param annot_angle Angle at which to plot the annotation text (names of
#'  datasets). Default `'auto'` to automatically decide. Otherwise, a numeric.
#'@param ... Additional arguments passed on to [display_subtree()].
#'@return A [ggtree::ggtree()] plot object, branches highlighted as specified in
#'  the `highlight_by` argument to [display_subtree()], with three layers of
#'  heatmap annotation at the tree tips. The innermost layer represents the
#'  number of entities in each tip label (or group of tip labels) in list
#'  `data_1`. The second (middle) layer represents the number of entities in
#'  each tip label (or group of tip labels) for the list `data_2`. The outermost
#'  layer represents the fraction of overlap between the entities at each tip
#'  label (or group of tip labels), calculated as (size of intersection)/(size
#'  of union).
#'@import ggplot2
#'@export
display_overlap <- function(base_tree = chemont_tree,
                            base_name = NULL,
                            data_1,
                            name_1 = 'Set1',
                            data_2,
                            name_2 = 'Set2',
                            entity_id_col = NULL,
                            group_level = "terminal",
                            tax_level_labels = chemont_tax_levels,
                            annot_angle = "auto",
                            ...){
  terminal_label <- NULL

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
    n_1 <- NULL
    n_2 <- NULL
    simil <- NULL
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
                                   entity_id_col = entity_id_col,
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
                            geom = ggplot2::geom_tile,
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


#' @title Add clade labels
#'
#' @description Add clade labels to an existing [ggtree::ggtree()] plot.
#'
#' @details Adds clade labels to an existing [ggtree::ggtree()] plot. Wrapper
#'   for [ggtree::geom_cladelab()] that handles some fiddly data wrangling, and
#'   allows the option to plot clade labels with alternating thick/thin bar
#'   widths, which is useful to distinguish clades when they are separated by
#'   little or no space on the plot.
#'
#'   Currently, [add_cladelab()] cannot be chained using the
#'   `+` operator. This is because [add_cladelab()] is not
#'   currently defined as an S3 class with an associated [ggplot2::ggplot_add()]
#'   method.
#'
#'   #'For example, the following code will *not* work:
#'
#' ```R ggtree(tree) + add_cladelab() ```
#'
#' [add_cladelab()] also cannot be chained using the
#' \code{\link[magrittr]{%>%}} operator if the preceding chain involves
#' \code{`+`}. This is because of operator precedence: R evaluates
#' \code{\link[magrittr]{%>%}} before \code{\link[ggplot2]{\%+\%}}.
#'
#'
#' For example, the following code also will *not* work:
#'
#' ```R ggtree(tree) + layout_circular() %>% add_cladelab() ```
#'
#'   However, the following code *will* work:
#'
#' ```R ggtree(tree) %>% add_cladelab() ```
#'
#' # Clade options
#'
#' The `cladeopt` parameter may include any of the "additional parameters"
#' understood by [ggtree::geom_cladelab()], including
#'
#' * `offset`
#' * `offset.text`
#' * `align`
#' * `extend`
#' * `angle`
#' * `horizontal`
#' * `barsize`
#' * `barcolour`
#' * `fontsize`
#' * `textcolour`
#' * `imagesize`
#' * `imagecolour`
#'
#' Additionally, in [add_cladelab()], the following parameters are
#' understood:
#'
#' * `wrap`: the number of characters to wrap text labels. Default 20.
#'    Use NULL for no wrapping.
#' * `barsize = "alternate"`: Alternates thick and thin bars. This is the
#'    default. It is useful to distinguish tightly-packed clades.
#' * `default_to_tip` Logical: Whether to default to showing tip labels if
#'    no clade label exists at the specified level. For example, if a tip
#'     terminates at level 3, but `clade_level = 4`, should that tip be
#'     labeled with its tip label (`default_to_tip = TRUE`), or should
#'     that tip be unlabeled (`default_to_tip = FALSE`)? The default
#'     value is `TRUE`.
#' * `draw_text` Logical: Whether to draw the text of clade labels,
#'    or only the bars. If `TRUE`, draw text. If `FALSE`, do not. Default `TRUE`.
#'
#'
#' @param tree_plot A [ggtree::ggtree()] plot object, e.g. the output
#'  of [display_subtree()] or of [ggtree::ggtree()].
#' @param tree Optional: the `phylo` tree object plotted in
#'  `tree_plot`. Function will run faster if `tree` is provided;
#'  otherwise `tree` will be reconstructed from `tree_plot$data` using
#'  [generate_taxonomy_tree()].
#' @param clade_level Numeric or `"auto"`: the taxonomy level at which to
#'  label clades, where level 0 is the root. Default value `"auto"` selects
#'  the taxonomy level of the MRCA of the tips plus one, or level 2, whichever
#'  is the greater number. If `clade_level = NULL`, no clade labels will be
#'  added and `tree_plot` will be returned unchanged.
#'@param clade_opts A named list of parameters controlling the appearance of
#'  clade labels. See "Details."
#'@return `tree_plot` with clade labels added.
#'@export
add_cladelab <- function(tree_plot,
                         tree = NULL,
                         clade_level = "auto",
                         clade_opts = list(wrap = 20,
                                           barsize = "alternate",
                                           fontsize = 3,
                                           lineheight = 0.7,
                                           default_to_tip = TRUE,
                                           draw_text = TRUE))
{
  phylo_node <- NULL
  clade_name <- NULL
  clade_name2 <- NULL

  if(!is.null(clade_level)){
    clade_opts_default <- list(wrap = 20,
                               barsize = "alternate",
                               fontsize = 3,
                               lineheight = 0.7,
                               default_to_tip = TRUE,
                               draw_text = TRUE)
    #keep defaults for any clade_opts not specified
    clade_opts <- c(clade_opts,
                    clade_opts_default[setdiff(names(clade_opts_default),
                                               names(clade_opts))])

    #if tree is not passed explicitly
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


    #plot clade bars with alternating widths by default
    #to do this:
    #first need to get order in which clades are plotted.
    #start with order in which *tips* are plotted.
    #ggtree:get_taxa_name() gives us tips in plotting order
    tips_plot <- ggtree::get_taxa_name(tree_view = tree_plot)
    #get clade label corresponding to each of these tips, at the specified level.
    #if there is no clade at the specified level, it will be NA.
    #first get node ID for the clade of each tip at specified level.
    #node ID will be NA if no clade at the specified level.
    clade_plot <- get_clade(node = get_node_from_label(label = tips_plot,
                                                       tree = tree),
                            tree = tree,
                            level = clade_level)
    #now get labels for the clade node IDs. Will be NA if no clade at specified level.
    clade_plot_lab <- get_label_from_node(node = clade_plot,
                                          tree = tree)

    #if we should print tip labels when there is no clade at the specified level:
    if(clade_opts$default_to_tip %in% TRUE){
      #replace NA clade labels with the corresponding tip labels.
      clade_plot[is.na(clade_plot)] <- get_node_from_label(label = tips_plot[is.na(clade_plot)],
                                                           tree = tree)
    }else{ #if we should not print tip labels when no clade at specified level:
      #remove any NA values from clade_plot
      clade_plot <- clade_plot[!is.na(clade_plot)]
      #remove any clade labels that are the *same* as tip labels
      clade_plot <- clade_plot[!(clade_plot_lab %in% tips_plot)]
    }

    #keep only the unique clades, in plotting order corresponding to tips
    clade_plot <- unique(clade_plot)

    #assign alternating bar widths in plotting order
    #create data frame with the node IDs, labels, and bar widths
    clade_dat <- data.frame(phylo_node = clade_plot,
                            clade_name = get_label_from_node(node = clade_plot,
                                                             tree = tree),
                            barsize = rep(1:2,
                                          length.out = length(
                                            clade_plot
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

    #we had it in tip plotting order; now put it in base tree order.
    dat3 <- clade_dat[match(intersect(dat$node,
                                      clade_dat$phylo_node),
                            clade_dat$phylo_node), ]


    if(clade_opts$barsize %in% "alternate"){
      #Re-sort the alternating bar sizes in clade plotting order
      #to correspond to base tree order of clades
      clade_opts$barsize <- dat3$barsize
    }

    #

    #if draw_text is FALSE, set all clade label text to blank
    if(clade_opts$draw_text %in% FALSE){
      dat3$clade_name2 <- ""
    }else{ #if draw_text is TRUE
      if(!is.null(clade_opts$wrap)){ #if wrapping is turned on
        #wrap clade names to width clade_opts$wrap
        dat3$clade_name2 <- stringr::str_wrap(dat3$clade_name,
                                              clade_opts$wrap)
      }else{
        dat3$clade_name2 <- dat3$clade_name2
      }
    }

    tree_plot +
      do.call(geom_cladelab,
              args = c(list(data = dat3,
                            mapping = aes(node = phylo_node,
                                          label = clade_name2,
                                          group = clade_name)),
                       clade_opts))
  }else{
    tree_plot
  }
}


#' Side by side trees
#'
#' This function takes two data sets and constructs both trees and internal
#' layers for data visualization. The left and right out layer of the diagram
#' consists of subtrees corresponding to the left and right data sets. The left
#' and right inner layers display label statistics for the left and right data
#' sets. The center displays label statistics shared by the left and right data
#' sets.
#'
#' @param base_tree The tree to be pruned, as a `phylo`-class
#'   object.
#' @param data_left The left data.frame(or data.table) of chemical
#'   classifications.
#' @param data_right The right data.frame (or data.table) of chemical
#'   classifications.
#' @param name_left Alternate parameter for left data set.
#' @param name_right Alternate parameter for right data set.
#' @param log_trans Alternate parameter for log-transforming data.
#' @param tax_level_labels Vector of the possible taxonomy levels that can
#'   appear as column names in `data_left` and `data_right`` of
#'   classified data.
#' @return An `aplot` consisting of two outer layer `ggtree` objects and three
#'   inner layer `ggplot2` objects.
#' @export
#' @import ggplot2
#' @import ggtree
side_by_side_trees <- function(base_tree = chemont_tree,
                               data_left,
                               data_right,
                               name_left = 'Left tree',
                               name_right = 'Right tree',
                               log_trans = FALSE,
                               tax_level_labels = chemont_tax_levels){
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

    data_left <- data.table::data.table(data_left)
    data_right <- data.table::data.table(data_right)

    terminal_labels_left <- data_left[!is.na(terminal_label), unique(terminal_label)]
    terminal_labels_right <- data_right[!is.na(terminal_label), unique(terminal_label)]

    terminal_labels <- union(terminal_labels_left, terminal_labels_right)

    left_initial_tree <- prune_tree(tree = base_tree, prune_to = data_left,
                                    tax_level_labels = tax_level_labels)
    right_initial_tree <- prune_tree(tree = base_tree, prune_to = data_right,
                                     tax_level_labels = tax_level_labels)

    left_labels <- c(left_initial_tree$tip.label, left_initial_tree$node.label)
    right_labels <- c(right_initial_tree$tip.label, right_initial_tree$node.label)

    all_labels <- union(left_labels, right_labels)

    union_tree <- drop_tips_nodes(tree = chemont_tree,
                                  labels = all_labels,
                                  keep_descendants = FALSE)
    union_tree$edge.length <- adjust_branch_lengths(union_tree)
    terminal_labels <- intersect(terminal_labels, union_tree$tip.label)

    left_ancestors <- lapply(left_labels, function(t) {get_ancestors(chemont_tree, t)})
    right_ancestors <- lapply(right_labels, function(t) {get_ancestors(union_tree, t)})

    all_left_tree <- union(unlist(left_ancestors), left_labels)
    all_right_tree <- union(unlist(right_ancestors), right_labels)

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
    value <- c(tree_data_leftval, tree_data_rightval, tree_data_centerval)
    tree_data <- cbind(tree_data, value)
    names(tree_data)[[3]] <- "value"

    trans <- ifelse(log_trans, 'log1p', 'identity')

    data_plot <- ggplot(tree_data, aes(x = tree, y = tip.label)) +
      geom_tile(aes(fill = value)) +
      ggplot2::scale_fill_viridis_c(trans=trans) +
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

#' @title Circular tree with boxplots
#'
#' @description Plot a circular tree annotated with boxplots.
#'
#' @details This function takes in a `data.frame` of classified entities with an
#'   additional numeric data variable. The `data.frame` will be grouped
#'   according to the terminal classification of the entities, and the specified
#'   numeric data will be plotted as boxplot annotations at the corresponding
#'   tips of the tree.
#'
#' @param data A `data.frame` of classified entities with at least one
#'   additional numeric data variable. Must contain at least the following
#'   variables: variables named for every element of `tax_level_labels`,
#'   containing the (character) labels at each level of classification; and the
#'   variable named in `col`.
#' @param col Character: the name of the numeric data variable to make boxplots
#'   with.
#' @param tax_level_labels Character: the taxonomy levels. Default
#'   [chemont_tax_levels].
#' @param title Character: the title of the plot.
#' @param tree A `phylo`-class tree object specifying the tree to plot.
#' @param layers Which taxonomic layers to display outside of the boxplot layer.
#'   This can be either a string with a single column name, a vector of column
#'   names, or a list of column names. If the input is a vector or a list, it is
#'   fine for it to be length 1.
#' @param adjust_branch_length `TRUE`/`FALSE`: If `TRUE`, adjust branch length
#'   so that all newly-pruned terminal nodes appear at the same length as tips,
#'   even if they were originally internal nodes. Default `FALSE`.
#' @param tippoint_boxplot `TRUE`/`FALSE`: Whether to color the tip points and
#'   the boxplots. If `TRUE`, tip points and boxplots will have matching colors.
#'   If `FALSE` (default), the boxplots will be filled white with black outline
#'   and the tip points will be black.
#' @return A [ggtree::ggtree()] object, plotting the subtree induced by the data, annotated with boxplots
#'   corresponding to the specified numeric column of data.
#' @export
#' @import ggplot2
#' @import ggtree
#' @import ggtreeExtra
#' @import data.table
circ_tree_boxplot <- function(data,
                              col,
                              tax_level_labels = chemont_tax_levels,
                              title = NULL,
                              tree = chemont_tree,
                              layers = NULL,
                              adjust_branch_length = FALSE,
                              tippoint_boxplot = FALSE){
  val <- NULL
  terminal_label <- NULL
  grp <- NULL
  Label <- NULL
  ID <- NULL



  if (!(col %in% names(data)))
    stop(paste('The column', col, 'is not in the input data!'))

  if (!('terminal_label' %in% names(data))){
    data <- add_terminal_label(copy(data))
  }

  new_data <- data[, c(col, 'terminal_label')]
  new_data[[col]] <- as.numeric(new_data[[col]])
  #drop NAs
  new_data <- new_data[!is.na(new_data[[col]]), ]
  #give harmonized names
  new_data[["grp"]] <- new_data$terminal_label
  new_data[["val"]] <- new_data[[col]]
  new_data[["node"]] <- new_data$terminal_label

  #prune tree to supplied data
  new_data_tree <- prune_tree(
    prune_to = data,
    tax_level_labels = tax_level_labels,
    tree = tree,
    adjust_branch_length = adjust_branch_length
  )

  num_nodes_tips <- length(new_data_tree$tip.label) + new_data_tree$Nnode

  tip_node_data <- data.frame('ID' = c(new_data_tree$tip.label,
                                       new_data_tree$node.label),
                              'Label' = c(new_data_tree$tip.label,
                                          new_data_tree$node.label))
  tip_node_data$Label <- factor(tip_node_data$Label)

  circ_plot <- ggtree(new_data_tree,
                      layout = 'circular')
  if (tippoint_boxplot) {
    circ_plot <- circ_plot %<+% tip_node_data +
      geom_tippoint(aes(color = Label), show.legend = FALSE) +
      ggplot2::scale_color_viridis_d(name = 'Terminal label',
                                     option = 'magma')#

    circ_plot <- circ_plot +
      ggtreeExtra::geom_fruit(data = new_data,
                              geom = geom_boxplot,
                              mapping = aes(x = val,
                                            y = terminal_label,
                                            fill = grp),
                              size = 0.2,
                              offset = 0.4,
                              outlier.size = 0.5,
                              outlier.stroke = 0.08,
                              outlier.shape = 21,
                              axis.params = list(axis = 'x',
                                                 text.size = 1.8,
                                                 text.angle = 270,
                                                 hjust = 0),
                              grid.params = list(),
                              show.legend = FALSE) +
      ggplot2::scale_fill_viridis_d(name = 'Group label',
                                    option = 'magma') +
      ggnewscale::new_scale_fill() +
      ggnewscale::new_scale_color()
  } else {
    circ_plot <- circ_plot + geom_tippoint()
    circ_plot <- circ_plot +
      ggtreeExtra::geom_fruit(data = new_data,
                              geom = geom_boxplot,
                              mapping = aes(x = val,
                                            y = terminal_label),
                              size = 0.2,
                              offset = 0.4,
                              outlier.size = 0.5,
                              outlier.stroke = 0.08,
                              outlier.shape = 21,
                              axis.params = list(axis = 'x',
                                                 text.size = 1.8,
                                                 text.angle = 270,
                                                 hjust = 0),
                              grid.params = list(),
                              show.legend = FALSE) +
      ggnewscale::new_scale_fill()
  }

  if (!is.null(layers)){
    label_levels <- get_labels(data, tax_level_labels = tax_level_labels)
    if (is.list(layers) | is.vector(layers)){
      level_names <- names(label_levels)[which(names(label_levels) %in% layers)]
    } else {
      warning('The `layers` parameter must be a list or a vector! Skipping extra layers for now...')
      level_names <- c()
    }

    fruit_data <- data.frame('ID' = c(new_data_tree$tip.label,
                                      new_data_tree$node.label))
    palettes <- c('Blues', 'Oranges', 'BuGn',
                  'OrRd', 'BuPu', 'Reds','GnBu',
                  'RdPu','Greens', 'YlOrBr',
                  'PuBu', 'YlOrRd', 'PuBuGn',
                  'YlGnBu', 'PuRd', 'YlGn', 'Purples', 'Greys')

    for (i in rev(seq_along(level_names))){
      level_index <- which(names(data) %in% level_names[[i]])
      values <- unname(
        as.list(
          data[, unique(.SD),
               .SDcol = level_names[[i]]])
      )[[1]]
      values <- values[!is.na(values)]
      empty_strings <- which(sapply(values, function(t) {t == ''}))
      if (length(empty_strings) > 0){
        values <- values[-which(sapply(values, function(t) {t == ''}))]
      }

      tree_nodes <- lapply(
        c(
          new_data_tree$tip.label,
          new_data_tree$node.label),
        function(x) {x}
      )

      for (j in seq_along(values)){
        total_descendants <- c(
          tree$tip.label,
          tree$node.label)[
            c(
              phangorn::Descendants(
                tree,
                which(
                  c(
                    tree$tip.label,
                    tree$node.label
                  ) %in% values[[j]]
                ),
                type = 'all'
              ),
              which(
                c(
                  tree$tip.label,
                  tree$node.label
                ) %in% values[[j]]
              )
            )
          ]
        name_indices <- which(
          c(
            new_data_tree$tip.label
            , new_data_tree$node.label) %in% total_descendants
        )

        names(tree_nodes)[name_indices] <- values[[j]]

      }
      level_number <- length(unique(names(tree_nodes)))
      level_labels <- levels(factor(names(tree_nodes)))
      names(tree_nodes)[
        which(is.na(names(tree_nodes)))
      ] <- paste0('_',
                  level_names[[i]])

      fruit_data[[level_names[[i]]]] <- factor(names(tree_nodes))

      colors <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(
          n = 9,
          name = palettes[[i]])
      )
      current_palette <- colors(level_number)
      circ_plot <- circ_plot +
        geom_fruit(data = fruit_data,
                   geom = geom_tile,
                   mapping = aes(y = ID,
                                 x = .data[[level_names[[i]]]],
                                 fill = .data[[level_names[[i]]]]),
                   width = 3,
                   pwidth = 0,
                   color = 'white') +
        scale_fill_manual(values = current_palette,
                          labels = level_labels) +
        ggnewscale::new_scale_fill()
    }


  }

  if (is.character(title)){
    circ_plot <- circ_plot +
      ggtitle(title,
              subtitle = paste('Boxplots of the data from', col, 'column')) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
  return(circ_plot)
}

#' Leaf fraction subtree
#'
#' This function takes in two `data.frame`s, plots the subtree induced by the
#' first `data.frame` and colors the tips based on the proportion of chemicals
#' from the second `data.frame` that make up the chemicals from the first, grouped
#' by tip (or terminal) label.
#'
#' @param data_1 A `data.frame` of classified entities (including terminal labels).
#' @param data_2 A `data.frame` of classified entities (including terminal labels).
#' @param name_1 Character: name of first data set (used for plot legend).
#' @param name_2 Character: name of second data set (used for plot legend).
#' @param show_labels `TRUE`/`FALSE`: whether to show tip labels.
#' @param tax_level_labels Character vector of taxonomy levels. Default [chemont_tax_levels].
#' @param tree A `phylo`-class tree object defining the base tree. Default [chemont_tree].
#' @return A [ggtree::ggtree()] object.
#' @export
#' @import ggtree
#' @import ggplot2
#'
#' @seealso \code{\link{add_terminal_label}}
#'
leaf_fraction_subtree <- function(data_1,
                                  data_2,
                                  name_1 = 'Set 1',
                                  name_2 = 'Set 2',
                                  show_labels = FALSE,
                                  tax_level_labels = chemont_tax_levels,
                                  tree = NULL){

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
  terminal_labels <- unique(data_1[!is.na(data_1$terminal_label),
                                   "terminal_label"])

  # For each terminal_label value, determine the chemicals from data_2 that are
  # also in data_1. This checks using the INCHIKEY of each chemical.
  label_percentages <- sapply(terminal_labels,
                              function(t) {
                                data_1_chemicals <- unique(data_1[data_1$terminal_label == t, "INCHIKEY"])
                                data_2_chemicals <- unique(data_2[data_2$terminal_label == t, "INCHIKEY"])
                                shared_chemicals <- intersect(data_1_chemicals, data_2_chemicals)
                                return(length(shared_chemicals)/length(data_1_chemicals))
                              }
  )

  label_data <- data.frame('label' = terminal_labels,
                           'percentages' = unname(label_percentages),
                           'data_1_numbers' <- unname(sapply(terminal_labels,
                                                             function(t) {

                                                               length(
                                                                 unique(data_1[data_1$terminal_label == t, "INCHIKEY"])
                                                               )
                                                             }
                           )
                           ),
                           'data_2_numbers' <- unname(
                             sapply(terminal_labels,
                                    function(t) {
                                      length(
                                        intersect(unique(data_1[data_1$terminal_label == t,
                                                                "INCHIKEY"]),
                                                  unique(data_2[data_2$terminal_label == t,
                                                                "INCHIKEY"])
                                        )
                                      )
                                    }
                             ))
  )

  names(label_data)[3:4] <- c(paste(name_1, 'label numbers'),
                              paste(name_2, 'label numbers in', name_2)
  )

  data_1_tree <- prune_tree(prune_to = data_1,
                            tax_level_labels = tax_level_labels,
                            tree = tree)

  if (length(data_1_tree$tip.label) <= 200){
    tip_size = 3
  } else if (length(data_1_tree$tip.label) <= 500){
    tip_size = 1.5
  } else {
    tip_size = .5
  }

  #add annotation data
  tree_plot <- ggtree::ggtree(data_1_tree,
                              layout = 'circular') %<+%
    label_data

  tree_plot <- tree_plot +
    ggplot2::ggtitle(name_1) +
    theme(plot.title = element_text(hjust = 0.5))

  tree_plot <- tree_plot +
    ggtree::geom_tippoint(aes(color = percentages),
                          size = tip_size)
  tree_plot <- tree_plot +
    ggplot2::scale_color_viridis_c(name = paste0('Overlap %\nwith ',
                                                 name_2),
                                   option = 'plasma')

  if (show_labels){
    tree_plot <- tree_plot +
      ggtree::geom_tiplab(aes(color = percentages), size = 1)
  } else {
    tree_plot <- tree_plot +
      ggtreeExtra::geom_fruit(geom = geom_tile,
                              mapping = aes(color = percentages),
                              width = 10,
                              height = .1,
                              offset = 0.1)
  }

  return(tree_plot)

}
