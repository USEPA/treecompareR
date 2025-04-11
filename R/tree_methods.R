#' Count descendants
#'
#' Count the number of descendants for each tree node
#'
#' Count the number of descendants and direct children for each
#' node of the input tree.
#'
#' @param tree An object of class `phylo` (see [ape::read.tree()] for
#'   description of this class)
#' @return A `data.frame` with variables `node` (the node number), `descendants`
#'   (the number of descendants for that node), `children` (the number of direct
#'   children for that node), and `level`.
#' @export
#'
#' @examples
#'   #a randomly-generated tree
#'   tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'   #count descendants of each node in this tree
#'   count_descendants(tree = tree)
#'
count_descendants <- function(tree){
  if (!inherits(tree, 'phylo')){
    stop("Please input an object of 'phylo' class!")
  }

  node_levels <- get_levels(tree)
  descendant_numbers <- data.frame(node = 1:dim(node_levels)[[1]],
                                   descendants = 0,
                                   children = 0,
                                   level = node_levels[, 'level'])
  for (i in max(node_levels[, 'level']):0){
    # get the tips/internal nodes at level i
    current_indices <- node_levels[node_levels$level == i, 'node']
    # filter out tips
    current_indices <- current_indices[current_indices > length(tree$tip.label)]
    for (j in seq_along(current_indices)){
      # Find the children of the current node
      current_children <- tree$edge[tree$edge[, 1] == current_indices[[j]], 2]
      # Get number of children
      number_of_children <- length(current_children)
      # Get descendant numbers of children
      descendants_of_children <- descendant_numbers[current_children, 'descendants']
      # Calculate total descendants
      total_descendants <- number_of_children + sum(descendants_of_children)
      # Record in the data.frame
      descendant_numbers[current_indices[[j]], 'descendants'] <- total_descendants
      descendant_numbers[current_indices[[j]], 'children'] <- number_of_children

    }
  }
  return(descendant_numbers)
}


#' Get levels
#'
#' Get levels for each node in a tree
#'
#' Levels are counted starting with the root of the tree at level 0. The
#' immediate children of the root node are level 1. Then their children are
#' level 2, etc.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @return A `data.frame` with variables `node` (the node number) and `level`
#'   (the level of each node in the tree). There is one row for each node in
#'   `tree` (including both tips and internal nodes).
#' @export
#'
#' @examples
#'
#'   tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#'   get_levels(tree = tree)
#'
get_levels <- function(tree){

  if (!inherits(tree, 'phylo')){
    stop("Please input an object of 'phylo' class!")
  }

  if(!ape::is.rooted(tree)){
    if(dim(tree$edge)[[1]] > length(tree$tip.label)){# Check if tree is a star
      stop("Please input a rooted tree!")
    }
  }
  root_number = as.integer(length(tree$tip.label) + 1)
  node_levels <- data.frame(node = 1:(dim(tree$edge)[[1]] + 1),
                            level = 0)
  level <- 0
  current_level <- c(root_number)
  while(length(current_level) > 0){
    node_levels[current_level, 'level'] <- level
    children <- tree$edge[tree$edge[, 1] %in% current_level, 2]
    current_level <- children
    level <- level + 1
  }
  return(node_levels)
}

#' Tree as data.frame
#'
#' Get data.frame representation of phylo tree
#'
#' Given a tree, get a data.frame listing node numbers, node names, and
#' level of each node.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @return A data.frame with as many rows as the number of nodes in \code{tree},
#'   and three variables: "node" (the node number), "level" (the hierarchical
#'   level of each node, where the root is level 0), and "Name" (the label
#'   corresponding to each node number).
get_tree_df <- function(tree){

  if (!inherits(tree, 'phylo')){
    stop("Please input an object of 'phylo' class!")
  }

  #get node numbers & levels
  tree_df <- get_levels(tree)
  #get parents of each node
  tree_df$parent <- as.integer(unlist(phangorn::Ancestors(x = tree,
                                        node = tree_df$node,
                                        type = "parent")))
  #add labels
  tree_df$Name <- get_label_from_node(node = tree_df$node,
                                     tree = tree)
  #return
  return(tree_df)
}

#' Generate information content.
#'
#' Generate information content for a tree.
#'
#' Generates a data.frame of information content for the input
#' tree, where information content is as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}. The
#' data.frame also includes the depth of each tip and internal node, the number
#' of descendants of each node, and the number of children for each node.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @return Numeric: a vector of information content for each node of the tree.
#' @export
#'
#' @references \insertRef{seco2004intrinsic}{treecompareR}
#'
#' @examplesIf FALSE
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' generate_information_content(tree = tree)
#'
generate_information_content <- function(tree){
  if (!inherits(tree, 'phylo')){
    stop("Please input an object of 'phylo' class!")
  }

  if (!ape::is.rooted(tree)){
    if(dim(tree$edge)[[1]] > length(tree$tip.label)){# Check if tree is a star
      stop("Please input a rooted tree!")
    }
  }

    descendants <- sapply(phangorn::allDescendants(tree),
                          length)
    #tips will be listed as their own descendants -- remove these
    descendants[seq_along(tree$tip.label)] <- 0
    n_node <- length(descendants)
    log_descendants <- 1 - (log(1 + descendants)/log(n_node))
    return(log_descendants)

}


#' Attach information content
#'
#' This function attaches information content if missing for input tree. This
#' uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param log_descendants Alternate parameter determining type of information
#'   content to use.
#' @return phylo object with information content data.frame attached
#' @export
#'
#' @references
#' \insertRef{seco2004intrinsic}{treecompareR}
#'
#' @examplesIf FALSE
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' tree <- attach_information_content(tree = tree)
#'
attach_information_content <- function(tree, log_descendants = TRUE){
  if (!inherits(tree, 'phylo')){
    stop("Please input an object of 'phylo' class!")
  }

  if (!ape::is.rooted(tree)){
    if(dim(tree$edge)[[1]] > length(tree$tip.label)){# Check if tree is a star
      stop("Please input a rooted tree!")
    }

  }

  if (!is.null(tree$IC)){
    indices <- which(names(tree$IC) %in% c('node', 'descendants', 'children', 'level'))
    if (length(indices) < 4){
      warning('Missing columns in tree$IC!')
      return(tree)
    }
  }

  tree$IC <- generate_information_content(tree = tree)

  return(tree)
}

#' Tree level
#'
#' This function returns the tree level of the given node in a rooted tree.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param node Character (a node name) or integer (a node number).
#' @return The level of the node from the root of the tree.
#' @export
#'
#' @examplesIf FALSE
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' get_tip_level(tree = tree, node = 't2') #specify node by name
#' get_tip_level(tree = tree, node = 2) #specify node by number
get_tip_level <- function(tree,
                          node){

  if(is.character(node)){
    get_node_from_label(label = node,
                        tree = tree)
  }

  return(length(phangorn::Ancestors(x = tree,
                                    node = node)))
}

#' Check similarity inputs
#'
#' Check to make sure nodes are valid.
#'
#' This is a helper function for checking user input values within the
#' \code{\link{jaccard_similarity}}, \code{\link{resnik_similarity}},
#' \code{\link{lin_similarity}}, \code{\link{jiang_conrath_similarity}}, and
#' \code{link{similarity_matrix}} functions.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param node_1 First node of interest. Character (a node name) or integer (a
#'   node number).
#' @param node_2 Second node of interest. Character (a node name) or integer (a
#'   node number).
#'
#' @return Two-element numeric vector: A pair of node numbers corresponding to `node_1` and `node_2`.
check_similarity_inputs <- function(tree = NULL,
                                    node_1 = NULL,
                                    node_2 = NULL){
  if (is.null(tree) | !('phylo' %in% class(tree))){
    stop('Please input a `phylo` object for the tree parameter!')
  }

  tree_labels <- c(tree$tip.label, tree$node.label)

  if(length(node_1)>1){
    stop('Please input a single node for `node_1`!')
  }

  if(length(node_2)>1){
    stop('Please input a single node for `node_2`!')
  }

  if(is.character(node_1)){
    node_1 <- get_node_from_label(label = node_1,
                                  tree = tree)
  }

  if(is.character(node_2)){
    node_2 <- get_node_from_label(label = node_2,
                                  tree = tree)
  }

  if (1 <= min(c(node1, node2))){
    if (max(c(node1, node2)) <= length(tree_labels)){
      return(c(node1, node2))
    } else {
      stop('An input node is out of range!')
    }
  }
  stop('An input node is out of range!')
}

#' Jaccard similarity
#'
#' This function takes in a tree and two input nodes (either labels or node
#' numbers) and returns the Jaccard similarity values of the nodes based on the
#' tree structure. For each node, there is a set of labels along the unique path from the
#' root to the label. These sets are compared using Jaccard similarity. For more
#' information on Jaccard similarity, please consult
#' \href{https://en.wikipedia.org/wiki/Jaccard_index}{Jaccard Index}.
#'
#' For rapid calculation of similarity values of multiple pairs of nodes,
#' consider using \code{\link{similarity_matrix}}.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param node_1 First node of interest. Character (a node name) or integer (a
#'   node number).
#' @param node_2 Second node of interest. Character (a node name) or integer (a
#'   node number).
#' @return Numeric: The Jaccard similarity of the ancestry of the two nodes.
#' @export
#' @examples
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#' jaccard_similarity(tree = tree,
#' node_1 = 1,
#' node_2 = 2)
#'
#'
#' @references
#' \insertRef{pekar2002taxonomy}{treecompareR}
#'
#' \insertRef{pesquita2009semantic}{treecompareR}
#'
#' @seealso \code{\link{resnik_similarity}}, \code{\link{lin_similarity}},
#' \code{\link{jiang_conrath_similarity}}, \code{\link{similarity_matrix}}

jaccard_similarity <- function(tree = NULL,
                              node_1 = NULL,
                              node_2 = NULL){

  #even if node_1 and node_2 were given as labels,
  #they will be converted to node numbers in this step:
  nodes <- check_similarity_inputs(tree = tree,
                                   node_1 = node_1,
                                   node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1

  # Handle the case where the root is both input nodes.
  if (all(c(node1, node2) == root)){
    return(1)
  }

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  jaccard <- get_jaccard(node1 = node1, node2 = node2, tree_nodes = tree_nodes,
                       tree_parents = tree_parents)
  return(jaccard)

}

#' Resnik Similarity
#'
#' This function takes in a tree and two input nodes (either labels or node
#' numbers) and returns the Resnik similarity values of the nodes based on the
#' tree structure. The function uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}. For
#' rapid calculation of similarity values of several pairs of nodes, consider
#' using \code{\link{similarity_matrix}}.
#'
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param node_1 First node of interest. Character (a node name) or integer (a
#'   node number).
#' @param node_2 Second node of interest. Character (a node name) or integer (a
#'   node number).
#' @return The Resnik similarity value between the two input nodes.
#' @seealso \code{\link{jaccard_similarity}}, \code{\link{lin_similarity}},
#'   \code{\link{jiang_conrath_similarity}}, \code{\link{similarity_matrix}}
#' @examples
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 19)
#' resnik_similarity(tree = tree,
#' node_1 = 1,
#' node_2 = 2)
#'
#'
#' @export
#'
#' @references \insertRef{lin1998information}{treecompareR}
#'
#'   \insertRef{resnik1995using}{treecompareR}
#'

resnik_similarity <- function(tree = NULL,
                              node_1 = NULL,
                              node_2 = NULL){

  nodes <- check_similarity_inputs(tree = tree,
                                   node_1 = node_1,
                                   node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- count_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  resnik <- get_resnik(node1 = node1,
                       node2 = node2,
                       tree_nodes = tree_nodes,
                       tree_parents = tree_parents,
                       information_content = information_content)
  return(resnik)

}

#' Lin Similarity
#'
#' This function takes in a tree and two input nodes (either labels or node
#' numbers) and returns the Lin similarity values of the nodes based on the tree
#' structure. The function uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
#' For rapid
#' calculation of similarity values of several pairs of nodes, consider using
#' \code{\link{similarity_matrix}}.
#'
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param node_1 First node of interest. Character (a node name) or integer (a
#'   node number).
#' @param node_2 Second node of interest. Character (a node name) or integer (a
#'   node number).
#' @return The Lin similarity value between the two input nodes.
#' @seealso \code{\link{jaccard_similarity}}, \code{\link{resnik_similarity}},
#'   \code{\link{jiang_conrath_similarity}}, \code{\link{similarity_matrix}}
#'
#' @examples
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 19)
#' lin_similarity(tree = tree,
#' node_1 = 1,
#' node_2 = 2)
#'
#' @export
#'
#' @references \insertRef{lin1998information}{treecompareR}
#'

lin_similarity <- function(tree = NULL,
                              node_1 = NULL,
                           node_2 = NULL){

  nodes <- check_similarity_inputs(tree = tree,
                                   node_1 = node_1,
                                   node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1

  # Handle the case where the root is both input nodes.
  if (all(c(node1, node2) == root)){
    return(1)
  }

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- count_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  lin <- get_lin(node1 = node1,
                 node2 = node2,
                 tree_nodes = tree_nodes,
                       tree_parents = tree_parents,
                 information_content = information_content)
  return(lin)

}

#' Jiang Conrath Similarity
#'
#' This function takes in a tree and two input nodes (either labels or node
#' numbers) and returns the Jiang and Conrath similarity values of the nodes
#' based on the tree structure. The function uses the formulation as described
#' in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
#' Furthermore, this function is a wrapper for a RCPP function. For rapid
#' calculation of similarity values of several pairs of nodes, consider using
#' \code{\link{similarity_matrix}}.
#'
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param node_1 First node of interest. Character (a node name) or integer (a
#'   node number).
#' @param node_2 Second node of interest. Character (a node name) or integer (a
#'   node number).
#' @return The Jiang Conrath similarity value between the two input nodes.
#' @seealso \code{\link{jaccard_similarity}}, \code{\link{resnik_similarity}},
#' \code{\link{lin_similarity}}, \code{\link{similarity_matrix}}
#'
#' @export
#'
#' @references \insertRef{seco2004intrinsic}{treecompareR}
#'
#' \insertRef{jiang1997semantic}{treecompareR}

jiang_conrath_similarity <- function(tree = NULL,
                              node_1 = NULL,
                              node_2 = NULL){

  nodes <- check_similarity_inputs(tree = tree,
                                   node_1 = node_1,
                                   node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1


  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- count_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  jiang_conrath <- get_jiang_conrath(node1 = node1,
                                     node2 = node2,
                                     tree_nodes = tree_nodes,
                       tree_parents = tree_parents,
                       information_content = information_content)
  return(jiang_conrath)

}

#' Similarity matrix
#'
#' This function takes in a list of node labels for rows and for columns (or
#' node indices), a tree, and a similarity measure, and returns a matrix with
#' the calculated similarity values for each pair. Note that the default is to
#' compute the similarity values only for pairs along and above the diagonal, so
#' when the set of nodes for rows and columns differ, one must turn off this
#' feature.
#'
#' To calculate the similarity matrix for all pairs of nodes in a tree, provide
#' `tree` and leave `labels1`, `labels2`, `nodes1`, and `nodes2` all `NULL`.
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class)
#' @param nodes1 A vector of the first set of nodes of interest. Character (node
#'   names) or integer (node numbers). If `NULL` (default), will be taken as all
#'   nodes in `tree`.
#' @param nodes2 A vector of the second set of nodes of interest. Character
#'   (node names) or integer (node numbers). If `NULL` (default), will be taken as all
#'   nodes in `tree`.
#' @param metric A string naming the similarity metric to use. Options are
#'   "jaccard", "resnik", "lin", and "jiang_conrath". Only the first two letters
#'   need be entered (e.g., "ja" for "jaccard", "re" for "resnik", "li" for
#'   "lin", and "ji" for "jiang_conrath"). Case-insensitive.
#' @param upper_tri Logical: If TRUE (default), only the upper triangular part
#'   of the similarity matrix is constructed and returned. If FALSE, the full
#'   matrix is constructed and returned.
#' @return A similarity matrix. The dimension names for the rows and columns
#'   correspond to the associated node labels of the tree.
#' @export
#'
#' @seealso \code{\link{jaccard_similarity}}, \code{\link{resnik_similarity}},
#'   \code{\link{lin_similarity}}, \code{\link{jiang_conrath_similarity}}
similarity_matrix <- function(tree = NULL,
                              nodes1 = NULL,
                              nodes2 = NULL,
                              metric = "jaccard",
                              upper_tri = TRUE){
  if (is.null(tree) | !('phylo' %in% class(tree))){
    stop('Please input a `phylo` object for the tree parameter!')
  }

  tree_labels <- c(tree$tip.label, tree$node.label)
  if(is.null(nodes1)){
    nodes1 <- get_node_from_label(tree_labels, tree)
  }

  if(is.null(nodes2)){
    nodes2 <- get_node_from_label(tree_labels, tree)
  }

  #convert labels to node numbers if necessary
  if(is.character(nodes1)){
    nodes1 <- get_node_from_label(nodes1, tree)
  }

  if(is.character(nodes2)){
    nodes2 <- get_node_from_label(nodes2, tree)
  }

  #convert metric to integer indicator (to pass to C++)
  metric <- substr(tolower(metric), 1, 2)
  metric_list <- c("ja",
                   "re",
                   "li",
                   "ji")
  sim_metric <- match(metric, metric_list)


  if (is.na(sim_metric) | !(sim_metric %in% 1:4)){
    stop('Please input a valid value for sim_metric parameter!')
  }

  root <- length(tree$tip.label) + 1

  nodes1 <- nodes1[nodes1 != root]
  nodes2 <- nodes2[nodes2 != root]

  if(length(nodes1) %in% 0){
    stop('nodes1 length 0 after removing root')
  }

  if(length(nodes2) %in% 0){
    stop('nodes2 length 0 after removing root')
  }

  if(any(!is.finite(nodes1))){
    stop('One or more of nodes1 invalid (non-finite)')
  }

  if(any(!is.finite(nodes2))){
    stop('One or more of nodes2 invalid (non-finite)')
  }

  if(any(nodes1<1)){
    stop('One or more of nodes1 less than 1')
  }

  if(any(nodes2<1)){
    stop('One or more of nodes2 less than 1')
  }

  if(any(nodes1 > length(tree_labels))){
    stop('One or more of nodes1 greater than total nodes + tips in tree')
  }

  if(any(nodes2 > length(tree_labels))){
    stop('One or more of nodes2 greater than total nodes + tips in tree')
  }

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- count_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  similarity <- get_similarity(nodes1 = nodes1,
                               nodes2 = nodes2,
                               tree_nodes = tree_nodes,
                               tree_parents = tree_parents,
                               sim_metric = sim_metric,
                               information_content = information_content,
                               upper_tri = upper_tri)

  dimnames(similarity) <- list(tree_labels[nodes1],
                               tree_labels[nodes2])

  return(similarity)

}

#' Random subtree similarity
#'
#' Random subtree similarity simulation
#'
#' Given a tree, randomly sample (with replacement) `n1` nodes in the tree and
#' treat these as terminal labels for a subtree (tree1). Repeat for `n2`,
#' generating a second randomly sampled subtree tree2. Using a pre-calculated
#' matrix of similarities for all nodes in the full tree, subset its rows to
#' keep only the labels from tree1, and its columns to keep only the labels from
#' tree2.
#'
#' @param n1 Number of nodes to sample for tree1.
#' @param n2 Number of nodes to sample for tree2.
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class). Default `chemont_tree`
#'   to use the full ChemOnt tree.
#' @param sim_tree A pre-computed similarity matrix for `tree`. Default
#'   [chemont_jaccard] to use the full Jaccard similarity matrix for
#'   [chemont_tree]. Other pre-computed options for [chemont_tree] include
#'   [chemont_resnik_IC_SVH], [chemont_lin_IC_SVH],
#'   [chemont_jiangconrath_IC_SVH].
#' @return A data.frame with similarity values for each simulation. The
#'   similarity values reported in each row is the mean similarity value for the
#'   corresponding data set/simulated tree given by the column.
#'
#' @examples
#' set.seed(42)
#' #average similarity between two random subtrees
#' #the same sizes as BIOSOLIDS2021 and USGS_WATER
#' MonteCarlo_similarity(n1 = nrow(biosolids_class), n2 = nrow(usgs_class))
#'
#' @author Paul Kruse, Caroline Ring
#' @export
#'
MonteCarlo_similarity <- function(n1,
                                  n2,
                                  fun = mean,
                                  tree = chemont_tree,
                                  sim_tree = chemont_jaccard){

  Nnode <- length(tree$node.label)
  treelabels <- c(tree$tip.label, tree$node.label[2:Nnode])

  #randomly select the specified number of labels, with replacement
  labs1 <- sample(treelabels, size = n1, replace = TRUE)
  labs2 <- sample(treelabels, size = n2, replace = TRUE)

  #subset the big similarity matrix
  sim_sub <- sim_tree[labs1, labs2]

  #compute average similarity for these two random subtrees
  return(mean(sim_sub))
}

#' Similarity cutoffs
#'
#' This function gives cutoffs for similarity values and percent of
#' representation of a data set. The function determines what percentage of a
#' data set that is represented by the induced subtree of the data set for
#' various values of a fixed similarity measure.
#'
#' @param mat A similarity matrix corresponding to a similarity measure and a
#'   rooted tree.
#' @param data A data.frame of classified entities.
#' @param tax_level_labels Parameter giving classification levels.
#' @param neighbors A parameter giving how many neighbors to use for finding
#'   label average values.
#' @param cutoff Numeric: the cutoff percentage value.
#' @param labels Character: a list of node labels
#'   corresponding to a subtree of a rooted tree.
#' @param counts Integer: the counts of occurrence for each
#'   label.
#' @return Named list of percentage of data represented by similarity values.
#'   The names are the similarity values. The values of the list are percentages
#'   of data represented by allowing similarity values equal to the names.
#' @export
#'
#' @examples
#' get_cutoffs(mat = chemont_jaccard, data = biosolids_class)
#' get_cutoffs(mat = chemont_jaccard, data = biosolids_class, neighbors = 6)
#'
get_cutoffs <- function(mat,
                        data,
                        tax_level_labels = chemont_tax_levels,
                        neighbors = 3,
                        cutoff = NA_real_,
                        labels = NULL,
                        counts = NULL){
    counts_df <- count_entities_per_label(data = data,
                                       tax_level_labels = tax_level_labels) %>%
      dplyr::bind_rows()

    labels <- counts_df[[1]]
    counts <- counts_df[[2]]

  if (!is.numeric(neighbors)){
    warning('Setting `neighbors` to have value 3...')
    neighbors = 3
  }

  if (neighbors - as.integer(neighbors) > 0){
    if (as.integer(neighbors) < 2){
      warning('Neighbors must be greater than 1! Setting value to 3...')
      neighbors = 3
    } else {
      neighbors <- as.integer(neighbors)
    }
  } else if (neighbors < 2){
    warning('Neighbors must be greater than 1! Setting value to 3...')
    neighbors = 3
  }

  indices <- which(dimnames(mat)[[1]] %in% labels)

  temp_mat <- mat[indices, indices]

  average_val <- unname(
    apply(temp_mat,
          MARGIN = 1,
          function(t) {
            sum(
              sort(t, decreasing = TRUE)[1:neighbors]
            )/neighbors
          }
    )
  )

  total = sum(counts)

  temp_counts <- counts[order(average_val, decreasing = TRUE)]

  unique_avgs <- sort(unique(average_val))

  percentages <- sapply(
    rev(unique_avgs),
                        function(t) {
                          sum(temp_counts[sort(average_val,
                                               decreasing = TRUE) >= t]
                              )/total
                          }
                        )

  names(percentages) <- rev(unique_avgs)

  if (!is.na(cutoff)) {
    # returns maximum value that achieves cutoff percentage threshold
    return(names(percentages)[[min(which(percentages >= cutoff))]])
  } else {
    return(percentages)
  }
  }

#' @title Drop tips and nodes
#'
#' @description Drop tips and nodes of a tree as needed to keep only the nodes
#'   of a specified subtree.
#'
#' @details This is a helper function used by [prune_tree()]. Usually it will
#'   not be called by the user directly.
#'
#'   This function takes in a `phylo`-class tree and one of several ways of
#'   specifying a subtree, and drops tips and nodes until the remaining tree
#'   consists solely of the nodes and tips associated with the specified subtree
#'   and their ancestors. The subtree can be specified using argument `data` as
#'   a `data.frame` of classified entities; using argument `label` as a
#'   character vector of subtree labels; using argument `nodes` as an integer
#'   vector of tip and node IDs; or using argument `level` as either a character
#'   string naming one of the taxonomy levels, or an integer giving the taxonomy
#'   level of interest (where root is level 0).
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class).
#' @param data A `data.frame` of classified entities. Default `NULL`.
#' @param labels Character: a set of labels of the subtree. Default `NULL`.
#' @param nodes Integer: a set of nodes of the subtree. Default `NULL`.
#' @param level Integer: the level to which the tree should be pruned. The root
#'   is level zero and each subsequent generation of children nodes is one level
#'   greater. Default `NULL`.
#' @param keep_descendants `TRUE`/`FALSE`: whether to keep all descendants of
#'   the input labels/nodes/level. Default `NULL`, which means behavior depends
#'   on whether `data`, `labels`, `nodes`, or `level` was provided. See
#'   [prune_tree()] for more details.
#' @param tax_level_labels Character vector naming the levels of the taxonomy in
#'   order from most general to most specific (excluding root). Default
#'   [chemont_tax_levels()].
#' @return A `phylo`-class object representing the induced subtree of the data.
#' @author Caroline Ring, Paul Kruse
#'
#' @examples
#' drop_tips_nodes(tree = chemont_tree,
#' data = biosolids_class[1:20,])
#'
#'
#' @references \insertRef{apepackage}{treecompareR}

drop_tips_nodes <- function(tree,
                            data = NULL,
                            labels = NULL,
                            nodes = NULL,
                            level = NULL,
                            keep_descendants = NULL,
                            tax_level_labels = chemont_tax_levels){
  if (!is.null(data)){
    #get terminal labels for each item in this data set
    #these are the labels to keep
    if(is.null(keep_descendants)){
      keep_descendants = FALSE
    }
    if (!"data.frame" %in% class(data)){
      stop("Input parameter `data` must be a data.table or a data.frame!")
    }
    tip_node_labels <- get_terminal_labels(data = as.data.frame(data),
                                  tax_level_labels = tax_level_labels)
    if(isTRUE(keep_descendants)){
      input_nodes <- get_node_from_label(label = tip_node_labels,
                                         tree = tree)
      tip_nodes <- phangorn::Descendants(x = tree,
                                         node = input_nodes)
      tip_node_labels2 <- get_label_from_node(node = unlist(tip_nodes),
                                             tree = tree)
      tip_node_labels <- union(tip_node_labels,
                               tip_node_labels2)
    }
  } else if (!is.null(labels)){
    if(is.null(keep_descendants)){
      keep_descendants <- TRUE
    }
    #get all descendants of input labels
    #in case input labels are internal nodes (i.e. "clades" to keep)
    input_nodes <- get_node_from_label(label = labels,
                                       tree = tree)
    if(isTRUE(keep_descendants)){
    tip_nodes <- phangorn::Descendants(x = tree,
                                             node = input_nodes)
    }else{
    tip_nodes <- input_nodes
    }
    tip_node_labels <- get_label_from_node(node = unlist(tip_nodes),
                                           tree = tree)
    tip_node_labels <- union(labels, tip_node_labels)
  } else if (!is.null(nodes)){ #if user specified one or more nodes to keep
    #keep nodes and their descendants
    if(is.null(keep_descendants)){
      keep_descendants <- TRUE
    }
    if(isTRUE(keep_descendants)){
      tip_nodes <- phangorn::Descendants(x = tree,
                                         node = nodes)
      tip_nodes <- c(nodes, unlist(tip_nodes))
    }else{
      tip_nodes <- nodes
    }

    tip_node_labels <- get_label_from_node(node = tip_nodes,
                                           tree = tree)
  }else if(!is.null(level)){ #is user has specified a level to prune to
    #if it's a numeric level, assume 0 = root, etc
    if(is.numeric(level)){
      #check to make sure it's a valid level
      #first check to make sure it is an integer
      if(!(as.integer(level)==level)){
        stop(paste("'level' was provided as", level,
                   "which is not a valid index for the vector of taxonomy levels",
                   "'tax_level_labels'.",
                   "If numeric, 'level' must be an integer between 1 and the length of",
                   "'tax_level_labels', which is",
                   length(tax_level_labels)))
      }
      #check to make sure it's within the length of tax_level_labels
      if(level > length(tax_level_labels)){
        stop(paste("'level' was provided as", level,
                   ", which is greater than the length of tax_level_labels,",
                   "which is,",
                   length(tax_level_labels),
                   ". tax_level_labels =",
                   paste(tax_level_labels, collapse = ", ")
        ))
      }
      #check to make sure it's not zero or negative
      if(level < 1){
        stop(paste("'level' was provided as", level,
                   "which is not a valid index for the vector of taxonomy levels",
                   "'tax_level_labels'.",
                   "'level' must be an integer between 1 and the length of",
                   "'tax_level_labels', which is",
                   length(tax_level_labels),
                   ". tax_level_labels =",
                   paste(tax_level_labels, collapse = ", ")))
      }
    }else if(is.character(level)){
    #if it's a string, match it with tax_level_labels
      #check to make sure it's a valid level
      if(level %in% tax_level_labels){
      level <- match(level, tax_level_labels)
      }else{
        stop(paste("'level' was provided as", level,
                   ", which is not one of the taxonomy level labels,",
                   "tax_level_labels =",
                   paste(tax_level_labels, collapse = ", ")
        ))
      }
    }
    #find all labels at the given level and do NOT keep descendants by default
    if(is.null(keep_descendants)){
      keep_descendants <- FALSE
    }

    #get tree as data.frame with level numbers for each node
    tree_df <- get_tree_df(tree)
    #get the labels of all nodes at the specified level
    tip_node_labels <- tree_df[tree_df$level %in% level, "Name"]

    #find their descendants only if told to do so
    #(this will result in just keeping the whole tree, which is silly)
    if(isTRUE(keep_descendants)){
      warning(paste("'level' was provided but keep_descendants = TRUE",
                    "which will result in keeping the whole tree,",
                    "and not dropping anything"))
      input_nodes <- get_node_from_label(label = tip_node_labels,
                                         tree = tree)
      tip_nodes <- phangorn::Descendants(x = tree,
                                         node = input_nodes)
      tip_node_labels2 <- get_label_from_node(node = unlist(tip_nodes),
                                              tree = tree)
      tip_node_labels <- union(tip_node_labels,
                               tip_node_labels2)
    }
  }else {
    stop(paste("Please input either a data.frame of chemical classifications",
               "in argument 'data',",
               "one or more tree node labels to keep in argument 'labels',",
               "one or more tree node numbers to keep in argument 'nodes',",
               "or a taxonomic level in argument 'level'!"))
  }

  max_depth <- max(get_levels(tree)$level)
  new_tree <- tree

#do it this way to retain cases where the terminal label was an internal node
#when we drop tips, newly-terminal internal nodes will be promoted to tip
#then we'll need to drop them, too

  for (i in 1:max_depth){
    labels_to_keep <- intersect(tip_node_labels,
                                new_tree$tip.label)
    labels_to_drop <- setdiff(new_tree$tip.label,
                              labels_to_keep)
    new_tree <- ape::drop.tip(new_tree,
                             labels_to_drop,
                              trim.internal = FALSE,
                              collapse.singles = FALSE)
  }


  return(new_tree)
}


#'@title Adjust branch lengths
#'
#'@description Adjust branch lengths for pruned trees.
#'
#'@details This is a helper function used by [prune_tree()]. Usually it will not
#'  be called by the user directly.
#'
#'  This function adjusts the length of branches after pruning. With the
#'  adjustment in this function, internal nodes that have become tip nodes after
#'  pruning (*i.e.*, because their descendants have been pruned away) will be
#'  plotted at the same branch length as other tip nodes.
#'
#'
#'@param tree An object of class `phylo` representing a rooted tree (see
#'  [ape::read.tree()] for description of this class).
#'@return Numeric: A vector of adjusted branch lengths for the input tree.
#'
#' @examples
#' my_tree <- prune_tree(tree = chemont_tree,
#'  prune_to = biosolids_class[1:10,])
#' adjust_branch_lengths(tree = tree)
#'
#'@author Paul Kruse, Caroline Ring

adjust_branch_lengths <- function(tree){
  tree_levels <- get_levels(tree)
  tree_height <- sapply(tree_levels$node, function(t){
    current_level <- tree_levels$level[[t]]
    descendants <- phangorn::Descendants(tree, t, type = 'all')
    max_level <- max(tree_levels$level[descendants])
    return(max_level - current_level + 1)
  })

  plot_height <- numeric(length(tree_height))
  root <- length(tree$tip.label) + 1
  plot_height[[root]] <- 100
  current_level <- c(root)
  while(length(current_level) > 0){
    next_level <- vector('list', length = length(current_level))
    for (i in seq_along(current_level)){
      parent_height <- plot_height[[current_level[[i]]]]
      children <- phangorn::Descendants(tree, current_level[[i]], type = 'children')
      if (length(children) > 0){
        next_level[[i]] <- children
      }

      if(length(children) > 0){
        scale_factor <- (tree_height[children] - 1)/tree_height[children]
        plot_height[children] <- parent_height*scale_factor
      }

    }
    current_level <- unlist(next_level)
  }


  edge.length <- numeric(dim(tree$edge)[1])
  for (i in seq_along(edge.length)){
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    edge.length[[i]] <- plot_height[[parent]] - plot_height[[child]]
  }
  return(edge.length)
}


#' @title Compare similarity measures
#'
#' @description Compare various similarity measures on a variety of tree shapes.
#'
#' @details This function compares the similarity measures of Jaccard, Resnik,
#'   Lin, and Jiang and Conrath on a variety of trees with different shapes. The
#'   tree types include the caterpillar, star, and balanced trees. The input
#'   parameter `n` indicates how many tips for the caterpillar and balanced
#'   tree. The star has 2n tips, and all trees have 2n-1 total nodes and tips.
#'
#' @param n Each tree has 2n-1 total nodes and tips.
#' @return A `data.frame` consisting of the mean self-similarity scores for each
#'   tree and similarity measure.
#' @author Paul Kruse
#' @examples
#' compare_similarity_measures(n=10)
#' compare_similarity_measures(n=20)
#'
compare_similarity_measures <- function(n){
  caterpillar <- generate_caterpillar(n)
  star <- generate_star(2*n)
  balanced <- generate_balanced(n)

  cat_labels <- c(caterpillar$tip.label, caterpillar$node.label)
  star_labels <- c(star$tip.label, star$node.label)
  balanced_labels <- c(balanced$tip.label, balanced$node.label)

  cat_IC <- attach_information_content(caterpillar)
  star_IC <- attach_information_content(star)
  balanced_IC <- attach_information_content(balanced)

  cat_Jaccard <- similarity_matrix(tree = caterpillar,
                                   nodes1 = cat_labels,
                                   nodes2 = cat_labels,
                                   metric = "jaccard")
  star_Jaccard <- similarity_matrix(tree = star,
                                    nodes1 = star_labels,
                                    nodes2 = star_labels,
                                    metric = "jaccard")
  balanced_Jaccard <- similarity_matrix(tree = balanced,
                                        nodes1 = balanced_labels,
                                        nodes2 = balanced_labels,
                                        metric = "jaccard")

  cat_Resnik <- similarity_matrix(tree = caterpillar,
                                  nodes1 = cat_labels,
                                  nodes2 = cat_labels,
                                  metric = "resnik")
  star_Resnik <- similarity_matrix(tree = star,
                                   nodes1 = star_labels,
                                   nodes2 = star_labels,
                                   metric = "resnik")
  balanced_Resnik <- similarity_matrix(tree = balanced,
                                       nodes1 = balanced_labels,
                                       nodes2 = balanced_labels,
                                       metric = "resnik")

  cat_Lin <- similarity_matrix( tree = caterpillar,
                                nodes1 = cat_labels,
                                nodes2 = cat_labels,
                                metric = "lin")
  star_Lin <- similarity_matrix(tree = star,
                                nodes1 = star_labels,
                                nodes2 = star_labels,
                                metric = "lin")
  balanced_Lin <- similarity_matrix(tree = balanced,
                                    nodes1 = balanced_labels,
                                    nodes2 = balanced_labels,
                                    metric = "lin")

  cat_JiangConrath <- similarity_matrix( tree = caterpillar,
                                         nodes1 = cat_labels,
                                         nodes2 = cat_labels,
                                         metric = "jiang")
  star_JiangConrath <- similarity_matrix(  tree = star,
                                           nodes1 = star_labels,
                                           nodes2 = star_labels,
                                           metric = "jiang")
  balanced_JiangConrath <- similarity_matrix( tree = balanced,
                                              nodes1 = balanced_labels,
                                              nodes2 = balanced_labels,
                                              metric = "jiang")

  simulation <- data.frame("Name (number of tips)" = c(paste("Caterpillar", n),
                                      paste("Star", 2*n),
                                      paste("Balanced", n)),
                           "Jaccard" = c(mean(cat_Jaccard[upper.tri(cat_Jaccard, diag = TRUE)]),
                           mean(star_Jaccard[upper.tri(star_Jaccard, diag = TRUE)]),
                           mean(balanced_Jaccard[upper.tri(balanced_Jaccard, diag = TRUE)])),
             "Resnik" = c(mean(cat_Resnik[upper.tri(cat_Resnik, diag = TRUE)]),
                          mean(star_Resnik[upper.tri(star_Resnik, diag = TRUE)]),
                          mean(balanced_Resnik[upper.tri(balanced_Resnik, diag = TRUE)])),
             "Lin" = c(mean(cat_Lin[upper.tri(cat_Lin, diag = TRUE)]),
                       mean(star_Lin[upper.tri(star_Lin, diag = TRUE)]),
                       mean(balanced_Lin[upper.tri(balanced_Lin, diag = TRUE)])),
             "JiangConrath" = c(mean(cat_JiangConrath[upper.tri(cat_JiangConrath, diag = TRUE)]),
                                mean(star_JiangConrath[upper.tri(star_JiangConrath, diag = TRUE)]),
                                mean(balanced_JiangConrath[upper.tri(balanced_JiangConrath, diag = TRUE)])))

  return(simulation)
}

#'@title Get clade
#'
#'@description Get clade (ancestor of a node at a specified level)
#'
#'@details Get the node numbers defining the clades (ancestors at a specified
#'  level) for specified input node numbers.
#'
#'@param node Integer: A vector of one or more node number(s) for which to get
#'  the clade(s)
#'@param tree An object of class `phylo` representing a rooted tree (see
#'  [ape::read.tree()] for description of this class).
#'@param level Integer: The hierarchical taxonomy level at which to get the
#'  clade(s). Root is level 0. Default value is 2 (superclass level, in
#'  ChemOnt).
#'@return Integer vector of node numbers representing the ancestors of the input
#'  nodes at the specified level.
#'@author Caroline Ring, Paul Kruse
#' @examples
#' get_clade(node = 35,
#' tree = chemont_tree,
#' level = 3)
#'
get_clade <- function(node,
                      tree,
                      level = 2){
  #get ancestors back to root for each input node
ancestors <- phangorn::Ancestors(x = tree,
                                 node = node,
                                 type = "all")
if(!is.list(ancestors)){
  ancestors <- list(ancestors)
}
#reverse the order in which ancestors are listed,
#so that root is listed first
ancestors <- lapply(ancestors, rev)
#add the node itself
ancestors <- lapply(seq_along(ancestors),
                    function(i) c(ancestors[[i]], node[i]))
#pull the ancestor at the specified taxonomy level (root = level 0)
clades <- sapply(ancestors, function(x) {
  if(length(x)>=(level+1)){
    x[level+1]
  }else{ #if there is no ancestor at that level, return NA
    NA_real_
  }
}
)

return(clades)
}

#' @title Get all clades
#'
#' @description List all clades in a tree at a specified level
#'
#' @details List all the labels in a specified tree at a specified level of the
#'   taxonomy.
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class).
#' @param level The level at which to display nodes (0 is the root)
#' @return A data.frame with four variables: `node` (the node number in the
#'   tree); `level` (the level of the node in the tree, where root is level 0);
#'   `parent` (the node number of the node's immediate parent); and `Name` (the
#'   text label of the node).
#' @author Caroline Ring
#' @examples
#' #list all ChemOnt superclasses:
#' get_all_clades(tree = chemont_tree, level = 2)
#'
#' @export
get_all_clades <- function(tree, level){
tree_df <- get_tree_df(tree = tree)
clade_df <- tree_df[tree_df$level %in% level, ]
return(clade_df)
}

#' @title Bind entities to a tree
#'
#' @description Bind individual classified entities as new tips to a tree
#'
#' @details Given a `data.frame` of classified entities with identifiers, this
#'   function gets the terminal classifications for each entity, and then binds
#'   the entity ID to each terminal classification as a new tip node of the
#'   tree.
#'
#'   This function may be useful when you plotting classified entities: for
#'   example, if you are using the ChemOnt taxonomy and plotting classified
#'   chemicals, and you want to show all of the individual chemicals with each
#'   classification in the tree, you could use [bind_entities()] first and then
#'   plot the resulting tree using, for example, [display_subtree()].
#'
#' @param tree An object of class `phylo` representing a rooted tree (see
#'   [ape::read.tree()] for description of this class). Tips will be bound to
#'   this tree.
#' @param data Either one `data.frame`, or a list of `data.frames`, containing
#'   classified entities. The data.frames must include the column names
#'   specifiedin `tax_level_labels` and `entity_id_col`.
#' @param entity_id_col Character vector: One or more variable name(s) in
#'   `data` that, together, uniquely identify the entities. Default `NULL`,
#'   which treats each row in `data` as a unique entity named `entity1`,
#'   `entity2`, ... for as many rows as there are in the data. If a vector of
#'   variable names is provided, entities will be named by concatenating rows of
#'   the specified variables.
#' @param tax_level_labels Taxonomy levels used for classification in
#'   `data`. Default is [chemont_tax_levels].
#' @return A \code{phylo}-class object.
#' @author Caroline Ring, Paul Kruse
#' @examples
#'
#' #explicitly specifying entity ID
#' bind_entities(tree = prune_tree(tree = chemont_tree,
#' prune_to = biosolids_class[1:20, ]),
#' data = biosolids_class[1:20, ],
#' entity_id_col = "DTXSID")
#'
#' #not specifying entity ID
#' bind_entities(tree = prune_tree(tree = chemont_tree,
#' prune_to = biosolids_class[1:20, ]),
#' data = biosolids_class[1:20, ])
#'
#' @export
bind_entities <- function(tree,
                          data,
                          entity_id_col = NULL,
                          tax_level_labels = chemont_tax_levels){

  #if a list of data frames is provided, rowbind it all together
  #this will be the "master list" of entities
  if(!is.data.frame(data)){
  if(is.list(data) &
     all(sapply(data, is.data.frame))){
    data <- dplyr::bind_rows(data)
  }
  }

  #if entity_id_col is NULL, then add a row ID
  #if no entity ID column specified,
  #create one with row numbers
  id_null <- FALSE
  if(is.null(entity_id_col)){
    id_null <- TRUE
    #add a new variable to the data
    #ensure it does not conflict with any of the existing variable names
    entity_id_col <- rev(
      make.names(
        names = c(names(data),
                  "id"
        ),
        unique = TRUE
      )
    )[1]
    data[[entity_id_col]] <- paste0("entity",
                                    1:nrow(data))
  }

  #if terminal label not already in data, add it
  if(!"terminal_label" %in% names(data)){
    data <- add_terminal_label(dat = data,
                               entity_id_col = entity_id_col,
                               tax_level_labels = tax_level_labels)
  }

  #get terminal labels for each entity
  #keep only unique entities & terminal labels
  term_labs <- unique(data[c(entity_id_col,
                        "terminal_label")])

  #get node numbers corresponding to terminal labels
  term_labs$terminal_node <- get_node_from_label(label = term_labs$terminal_label,
                                                 tree = tree)

  term_labs <- term_labs[!is.na(term_labs$terminal_node), ]



  #loop over terminal labels
  #create new tree
  #bind to the terminal label node
  #essentially, treat entity as a new level of classification
  tree_df <- get_tree_df(tree)
  new_df_list <- lapply(unique(term_labs$terminal_label),
         function(label){
          tmpdf <- term_labs[term_labs$terminal_label %in% label, ]
          #add these as new nodes whose parent is the terminal label node
          parent_node <- tree_df[tree_df$Name %in% label, "node"]
          parent_level <- tree_df[tree_df$Name %in% label, "level"]
          new_level <- parent_level + 1
          label_df <- data.frame(level = rep(new_level,
                                             nrow(tmpdf)),
                                 parent = rep(parent_node,
                                              nrow(tmpdf)),
                                 Name = do.call(paste,
                                                tmpdf[entity_id_col]
                                                )
          )
  })
#bind all the list of data.frames into one big one
  new_df <- dplyr::bind_rows(new_df_list)

  #find terminal nodes in original tree without any entities
  tips_no_ents <- setdiff(tree$tip.label, term_labs$terminal_label)
  if(length(tips_no_ents)>0){ #if any such entity-less tips
  #create some placeholder entities -- these will be deleted later
  term_fake <- data.frame(terminal_label = tips_no_ents,
                          Name = paste0("fake_entity_",
                                        tips_no_ents))
  fake_df_list <- lapply(term_fake$terminal_label,
                         function(label){
                           tmpdf <- term_fake[term_fake$terminal_label %in% label, ]
                           #add these as new nodes whose parent is the terminal label node
                           parent_node <- tree_df[tree_df$Name %in% label, "node"]
                           parent_level <- tree_df[tree_df$Name %in% label, "level"]
                           new_level <- parent_level + 1
                           label_df <- data.frame(level = rep(new_level,
                                                              nrow(tmpdf)),
                                                  parent = rep(parent_node,
                                                               nrow(tmpdf)),
                                                  Name = tmpdf$Name)
                         })

  fake_df <- dplyr::bind_rows(fake_df_list)

  new_df <- dplyr::bind_rows(new_df, fake_df)
  }

  #new node numbers
  new_df$node <- max(tree_df$node) + 1:nrow(new_df)

  new_df <- dplyr::bind_rows(tree_df,
                             new_df)

  #make into a tree
  new_df <- setNames(new_df,
           c("ID",
             "level",
             "Parent_ID",
             "Name"))
  new_tree <- generate_taxonomy_tree(new_df)

  if(length(tips_no_ents)>0){
  #drop fake entities
  new_tree <- ape::drop.tip(new_tree,
                            term_fake$Name,
                            trim.internal = FALSE,
                            collapse.singles = FALSE)
  }
  return(new_tree)

}

#' @title Prune a tree
#'
#' @description Prune a tree to keep only a specified subtree.
#'
#' @details
#' # How to specify the subtree to keep
#' `prune_to` defines the subtree to *keep* (everything else will be pruned
#' away). It may be specified in several different ways.
#'
#' ## As a data.frame
#'
#' If `prune_to` is a `data.frame` of classified entities, it must have
#' columns corresponding to, and named for, each of the taxonomy levels as
#' defined in the argument `tax_level_labels`, containing the labels at the
#' corresponding level for each entity. It must also have at least one more
#' column, uniquely identifying the entities; the name of the additional column
#' does not matter, as long as it is not the same as one of the taxonomy levels.
#' The result will be to keep only the subtree induced by this classified data
#' set, *i.e.*, only the branches of the tree that occur in this classified data
#' set. By default, any descendants of the node labels in the `data.frame`
#' that do not themselves appear in the `data.frame` will *not* be kept. If
#' you want to keep descendants that do not themselves appear in the
#' `data.frame`, specify \code{keep_descendants = TRUE}.
#'
#' ## As the name of a taxonomy level
#'
#' If `prune_to` is the name of a taxonomy level (one of the levels defined
#' in argument `tax_level_labels`), the result will be to keep only nodes at
#' that taxonomic level or less-specific levels. (For example, for the ChemOnt
#' taxonomy, specifying \code{prune_to = "class"} will keep only nodes at levels
#' "kingdom", "superclass", and "class". Any nodes at level "subclass", "level5",
#' "level6", ... "level11" will be dropped. (If you specify `prune_to` as
#' the name of a taxonomy level, and also specify \code{keep_descendants = TRUE},
#' the result will be to keep the whole tree.)
#'
#' ## As a vector of node/tip labels or numbers
#'
#' If `prune_to` is a vector of node/tip labels (i.e., labels appearing in
#' \code{tree$node.label} and/or \code{tree$tip.label}) or node/tip numbers (i.e.
#' node/tip index numbers between 1 and \code{ape::Ntip(tree) +
#' ape::Nnode(tree)}), the result will be to keep only the nodes/tips that are in
#' that vector, keep their common ancestors, and (by default) also keep their
#' descendants if any. The intention of keeping the descendants by default is to
#' allow the user to prune to specified clades simply by specifying the labels or
#' node numbers of the MRCAs of the clades. For example, using the ChemOnt
#' taxonomy, you could prune to keep all branches in the superclass
#' "Organohalogen compounds" by simply specifying \code{prune_to = "Organohalogen
#' compounds"}.  If you do *not* wish to keep the descendants of the specified
#' node labels/numbers, then specify \code{keep_descendants = FALSE}.
#'
#' @param tree The tree to be pruned, as a `phylo`-class object (see
#'   [ape::read.tree()] for description of this class).
#' @param prune_to What to *keep* from the base tree (everything else will be
#'  pruned away). May be a `data.frame` of classified data; one or more
#'  labels in the tree (tip or internal node labels); one or more node numbers
#'  in the tree (tip or internal nodes); or the name of a taxonomy level (one of
#'  the items in `tax_level_labels`). Default is NULL, which results in no
#'  pruning being done (i.e., the base tree is returned as-is). See Details.
#' @param keep_descendants Whether to keep descendants of what is specified in
#'  `prune_to`. Default NULL will choose the behavior based on the class of
#'  `prune_to`: when `prune_to` is a `data.frame` or one of the
#'  taxonomy level labels, \code{keep_descendants = FALSE} by default. When
#'  `prune_to` is a vector of node/tip labels or numbers in the base
#'  tree,\code{keep_descendants = TRUE} by default. See Details.
#' @param adjust_branch_length Whether to adjust branch length so that all
#'  newly-pruned terminal nodes appear at the same length as tips, even if they
#'  were originally internal nodes. Default FALSE.
#' @param tax_level_labels Vector of the possible taxonomy levels that can appear
#'  as column names in `prune_to` if it is a `data.frame` of
#'  classified data.
#' @param ... Additional arguments, not currently used.
#' @return A `phylo`-class object representing the pruned tree.
#' @author Caroline Ring, Paul Kruse
#' @examples
#'
#' #prune_to as a data.frame of classified entities
#' #prunes to the first 20 chemicals in BIOSOLIDS2021
#' prune_tree(tree = chemont_tree,
#' prune_to = biosolids_class[1:20, ])
#'
#' #prune_to as a character vector of labels
#' #prunes to the first 20 chemicals in BIOSOLIDS2021
#' prune_tree(tree = chemont_tree,
#' prune_to = biosolids_class[1:20, "terminal_label"])
#'
#' #prune_to as an integer vector of node numbers
#' #prunes to the first 20 chemicals in BIOSOLIDS2021
#' my_nodes <- get_node_from_label(label = biosolids_class[1:20, "terminal_label"],
#' tree = chemont_tree)
#' prune_tree(tree = chemont_tree,
#' prune_to = my_nodes)
#'
#' #prune_to as a taxonomy level
#' prune_tree(tree = chemont_tree,
#' prune_to = 2) #prunes to superclasses
#'
#'@export
prune_tree <- function(tree,
                       prune_to = NULL,
                       keep_descendants = NULL,
                       adjust_branch_length = FALSE,
                       tax_level_labels = chemont_tax_levels,
                       ...){
  if(!is.null(prune_to)){ #if user has specified something to prune to
    if(is.data.frame(prune_to)){ #if user has specified a dataset to prune to
      #Prune the tree according to the specified dataset
      if(is.null(keep_descendants)){
        keep_descendants <- FALSE
      }
      pruned_tree <- drop_tips_nodes(tree = tree,
                                     data = prune_to,
                                     keep_descendants = keep_descendants,
                                     tax_level_labels = tax_level_labels)
    }else if(is.character(prune_to)){
      #check if this is one of the tax_level_labels
      #if so, interpret it as a level to prune to
      if(all(prune_to %in% tax_level_labels)){
        #by default do NOT keep descendants
        #since that would just result in keeping the whole tree
        if(is.null(keep_descendants)){
          keep_descendants <- FALSE
        }
        if(isTRUE(keep_descendants)){
        warning(paste("'prune_to' =",
                      paste0('\"', prune_to, '\"'),
        " has been interpreted as a taxonomy level,",
        "because it is in 'tax_level_labels' = ",
        paste(tax_level_labels, collapse = ", "),
        "But 'keep_descendants = TRUE'",
                      "which will result in keeping the whole tree,",
                      "and not pruning anything"))
        }
        pruned_tree <- drop_tips_nodes(tree = tree,
                                       level = prune_to,
                                       keep_descendants = keep_descendants)
      }else{ #if not a tax_level_label,
      #interpret as node/tip labels
      #prune to only the subtree with this label(s)
        #including the descendents of internal label(s) by default
      if(is.null(keep_descendants)){
        keep_descendants <- TRUE
      }
      pruned_tree <- drop_tips_nodes(tree = tree,
                                     labels = prune_to,
                                     keep_descendants = keep_descendants)
    }
      }else if(is.numeric(prune_to)){
      #interpret as node numbers
      #prune to only the subtree with this node(s)
      #including the descendents of internal node(s) by default
      if(is.null(keep_descendants)){
        keep_descendants <- TRUE
      }
      pruned_tree <- drop_tips_nodes(tree = tree,
                                     nodes = prune_to,
                                     keep_descendants = keep_descendants)
    }

    if (adjust_branch_length) {
      pruned_tree$edge.length <- adjust_branch_lengths(pruned_tree)
    }

  }
  return(pruned_tree)
}

#' Convert phylo to classified data.frame
#'
#'Convert a phylo tree into a wide-format "classified" `data.frame`
#'
#'@param tree A `phylo`-class object (see
#'   [ape::read.tree()] for description of this class).
#'@param tax_level_labels Vector of the possible taxonomy levels that can appear
#'  as column names in \code{as_classified.phylo} if it is a `data.frame`
#'  of classified data. Default [chemont_tax_levels].
#'@return A `data.frame` with variables `tip_label` (naming the tips or
#'  entities) and one variable for each item in `tax_level_labels`, giving the
#'  classification for each entity at each taxonomic level.
#' @author Caroline Ring, Paul Kruse
#' @examples
#' #a tree pruned to only the first 10 BIOSOLIDS2021 chemicals
#' my_tree <- prune_tree(tree = chemont_tree, prune_to = biosolids_class[1:10, ])
#' #bind the DTXSIDs as new tips
#' my_tree <- bind_entities(tree = my_tree,
#'                          data = biosolids_class[1:10,],
#'                          entity_id_col = "CASRN")
#' #now convert this tree to a "classified" data.frame
#' as_classified.phylo(tree = my_tree)
#' #compare the result to the original classified data.frame
#' biosolids_class[1:10, c("CASRN", chemont_tax_level)]
#'
#'@export
as_classified.phylo <- function(tree,
                                tax_level_labels = chemont_tax_levels){

  tip_label <- NULL


  foo <- dplyr::bind_rows(
    lapply(1:(ape::Ntip(tree)),
         function(tip_node){
           ancestors <- phangorn::Ancestors(x= tree,
                                            node = tip_node,
                                            type = "all")
           ancestors_rev <- rev(ancestors)[-1] #delete root node
           ancestors_add_tip <- c(ancestors_rev, tip_node)
           ancestors_labels <-  get_label_from_node(node = ancestors_add_tip,
                                                    tree = tree)
           ancestors_levels <- tax_level_labels[seq_along(ancestors_labels)]
           return(data.frame(tip_label = tree$tip.label[tip_node],
                             labels = ancestors_labels,
                             levels = ancestors_levels))
         }
         )
         )

  foo2 <- tidyr::pivot_wider(foo,
                             id_cols = tip_label,
                             names_from = levels,
                             values_from = labels)


return(as.data.frame(foo2))

}

#' Get subtree node numbers
#'
#' This is a helper function that takes a classified data set and a base
#' taxonomy tree and provides all node numbers in the subtree corresponding to
#' the data set.
#'
#' @param data A classified data set.
#' @param base_tree The base tree, as a `phylo`-class object
#'   (see [ape::read.tree()] for description of this class). Default is
#'   \code{\link{chemont_tree}}, the full ChemOnt taxonomy tree.
#' @param tax_level_labels A vector of levels for the taxonomy. Default is
#'   \code{\link{chemont_tax_levels}}, the levels of the ChemOnt taxonomy.
#' @return A vector of node numbers in the base tree that are represented in the
#'   subtree corresponding to the data set.
#' @author Caroline Ring, Paul Kruse
#'   @examples
#'   get_subtree_nodes(data = biosolids_class[1:10,])
#'

get_subtree_nodes <- function(data,
                              base_tree = chemont_tree,
                              tax_level_labels = chemont_tax_levels){
  label <- NULL
  #get all labels represented by the classified dataset
  data_labels <- tidyr::pivot_longer(data,
                                     cols = dplyr::all_of(tax_level_labels),
                                     names_to = "level",
                                     values_to = "label") %>%
    dplyr::pull(label)

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

#' Get label from node
#'
#' Get label for a tip or internal node ID in a phylo tree
#'
#' @param node Integer vector of node ID numbers in phylo tree.
#' @param tree A `phylo`-class object (see
#'   [ape::read.tree()] for description of this class).
#' @return Character vector of tip or internal node labels corresponding to each
#'   node ID. \code{NA_character_} if no node label corresponds to the input
#'   node ID.
#'   @author Caroline Ring, Paul Kruse
#'   @examples
#'   get_label_from_node(node = 35, tree = chemont_tree)
#'   get_label_from_node(node = c(35, 42), tree = chemont_tree)
#'
#' @export
get_label_from_node <-function(node, tree){
  #get total number of nodes in the tree
  N <- dim(tree$edge)[[1]] + 1
  node[node<=0] <- N+100 #this will force return NA label for negative or 0 node IDs
  node[is.na(node)] <- N+100 #same for any NA nodes
  #get all tree labels: tips then nodes
  treelabs <- c(tree$tip.label,
                tree$node.label)
  #get labels corresponding to each node
  label <- treelabs[node]

  return(label)
}

#' Get node from label
#'
#' Get node ID for a label in a phylo tree
#'
#' @param label Character vector of labels for tips or internal nodes in phylo tree
#' @param tree A `phylo`-class object (see
#'   [ape::read.tree()] for description of this class).
#' @return Integer vector of tip or internal node ID numbers.
#' @author Caroline Ring, Paul Kruse
#' @examples
#' get_node_from_label(label = "Benzacridines", tree = chemont_tree)
#' get_node_from_label(label = c("Benzacridines","Dihydrofuranoquinolines"), tree = chemont_tree)
#'
#' @export
get_node_from_label <- function(label, tree){
  #tip labels come first
  tip_nodes <- match(label, tree$tip.label)
  #then internal node labels come
  internal_nodes <- match(label, tree$node.label) + ape::Ntip(tree)
  nodes <- pmin(tip_nodes, internal_nodes, na.rm = TRUE)
  return(nodes)

}
