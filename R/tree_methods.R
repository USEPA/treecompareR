# This generates the number of descendants for each node. The edges are stored
# in a two column data.frame, with the first column indicating the parent and
# the second column indicating the child. Internal numbering of nodes is given
# by tips first and then nodes, with the root given by the first node of the
# internal nodes. The output also includes the number of children for each node
# and the node level. The root defined as level 0 and the level of a child node
# is one greater than its parent node.

#' Generate descendants
#'
#' This function generates a data.frame of the number of descendants for each
#' node of the input tree.
#'
#' @param tree A phylo object
#' @return A data.frame consisting of the node number, descendants, children,
#'   and level.
#' @export
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' generate_descendants(tree = tree)
#'
generate_descendants <- function(tree){
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


#' Generate tree levels
#'
#' This function generates a data.frame of node levels for input rooted tree.
#' The data.frame has a row corresponding to each tip and internal node. It
#' records the depth from the root of each node.
#'
#' @param tree A phylo object representing a rooted tree.
#' @return data.frame consisting of the node number, and the level of each node
#' @export
#' @importFrom ape is.rooted
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' get_levels(tree = tree)
#'
get_levels <- function(tree){
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

#' Get data.frame representation of phylo tree
#'
#' Helper function to get a data.frame listing node numbers, node names, and
#' level of each node
#'
#' @param tree A phylo-class tree object
#' @return A data.frame with as many rows as the number of nodes in \code{tree},
#'   and three variables: "node" (the node number), "level" (the hierarchical
#'   level of each node, where the root is level 0), and "Name" (the label
#'   corresponding to each node number).
get_tree_df <- function(tree){
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
#' This function generates a data.frame of information content for the input
#' tree. This uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}. The
#' data.frame also includes the depth of each tip and internal node, the number
#' of descendants of each node, and the number of children for each node.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param log_descendants Alternate parameter for specifying type of information
#'   content.
#' @return data.frame consisting of node number, children, descendants, level,
#'   and information content for each node
#' @export
#' @importFrom ape is.rooted
#'
#' @references \insertRef{seco2004intrinsic}{treecompareR}
#'
#' @examples
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
#' @param tree A phylo object representing a rooted tree.
#' @param log_descendants Alternate parameter determining type of information
#'   content to use.
#' @return phylo object with information content data.frame attached
#' @export
#' @importFrom ape is.rooted
#'
#' @references
#' \insertRef{seco2004intrinsic}{treecompareR}
#'
#' @examples
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

#' Calculate similarity measures
#'
#' @param tree A phylo object representing a rooted tree.
#' @param labels_A The first (set of) label(s).
#' @param labels_B The second (set of) label(s).
#' @param metric A string naming the similarity metric to use. Options are
#'   "jaccard", "resnik", "lin", and "jiang_conrath". Only the first two letters
#'   need be entered (e.g., "ja" for "jaccard", "re" for "resnik", "li" for
#'   "lin", and "ji" for "jiang_conrath"). Case-insensitive.
#' @return Matrix of similarity values for pairs of labels from `label_A` and
#'   `label_B`.
#' @export
#'
#' @references \insertRef{pekar2002taxonomy}{treecompareR}
#'
#' \insertRef{pesquita2009semantic}{treecompareR}
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' calc_similarity(tree = tree, labels_A = 't2', labels_B = 't4')
#' calc_similarity(tree = tree, labels_A = 'n3', labels_B = 't8')
#'

calc_similarity <- function(tree,
                            labels_A,
                            labels_B,
                            metric = "jaccard"){


  #convert metric to integer indicator (to pass to C++)
  metric <- substr(tolower(metric), 1, 2)
  metric_list <- c("ja",
                   "re",
                   "li",
                   "ji")
  metric_int <- match(metric, metric_list)

  #convert labels A to node ID numbers
  node_A <- get_node_from_label(label = labels_A,
                                tree = tree)
  #convert labels B to node numbers
  node_B <- get_node_from_label(label = labels_B,
                                tree = tree)

root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  #Call C++ function
  outmat <- get_similarity(nodes1 = node_A,
                           nodes2 = node_B,
                           tree_nodes = as.integer(c(tree$edge[,2], root_node)),
                           tree_parents = as.integer(c(tree$edge[,1], 0)),
                           sim_metric = metric_int)
return(outmat)

}

#' Resnik similarity
#'
#' This determines the Resnik similarity for two input nodes in a given tree.
#' This uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
#'
#' @param tree A phylo object representing a rooted tree, with an information
#'   content attribute IC.
#' @param label_A The first node label.
#' @param label_B The second node label.
#' @param node_A Alternate parameter, the first node number.
#' @param node_B Alternate parameter, the second node number.
#' @return The Resnik similarity in the given tree of the pair of nodes.
#' @export
#'
#' @references
#' \insertRef{lin1998information}{treecompareR}
#'
#' \insertRef{resnik1995using}{treecompareR}
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' general_Resnik_similarity(tree = tree, label_A = 't2', label_B = 't4')
#' general_Resnik_similarity(tree = tree, label_A = 'n3', label_B = 't8')
#'
general_Resnik_similarity <- function(tree,
                                      label_A = NULL,
                                      label_B = NULL,
                                      node_A = NULL,
                                      node_B = NULL) {
  if (is.null(label_A) | is.null(label_B)){
    if (is.null(node_A) | is.null(node_B)){
      stop('Please input a pair of labels or a pair of node numbers')
    }else{
      #ensure they are integers
      node_A <- as.integer(node_A)
      node_B <- as.integer(node_B)
    }
  } else {
    #convert labels A to node ID numbers
    node_A <- get_node_from_label(label = label_A,
                                  tree = tree)
    #convert labels B to node numbers
    node_B <- get_node_from_label(label = label_B,
                                  tree = tree)
  }

  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  #Call C++ function
  resSim <- get_resnik(node1 = node_A,
                    node2 = node_B,
                    tree_nodes = as.integer(c(tree$edge[,2], root_node)),
                    tree_parents = as.integer(c(tree$edge[,1], 0)))
}


#' Lin similarity
#'
#' This determines the Lin similarity for two input nodes in a given tree. This
#' uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
#'
#' @param tree A phylo object representing a rooted tree, with an information
#'   content attribute IC.
#' @param label_A The first node label.
#' @param label_B The second node label.
#' @param node_A Alternate parameter, the first node number.
#' @param node_B Alternate parameter, the second node number.
#' @return The Lin similarity in the given tree of the pair of nodes.
#' @export
#'
#' @references
#' \insertRef{lin1998information}{treecompareR}
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' calc_Lin_similarity(tree = tree, label_A = 't2', label_B = 't4')
#' calc_Lin_similarity(tree = tree, label_A = 'n3', label_B = 't8')
#'
calc_Lin_similarity <- function(tree,
                                label_A = NULL,
                                label_B = NULL,
                                node_A = NULL,
                                node_B = NULL){
  if (is.null(label_A) | is.null(label_B)){
    if (is.null(node_A) | is.null(node_B)){
      stop('Please input a pair of labels or a pair of node numbers')
    }else{
      #ensure they are integers
      node_A <- as.integer(node_A)
      node_B <- as.integer(node_B)
    }
  } else {
    #convert labels A to node ID numbers
    node_A <- get_node_from_label(label = label_A,
                                  tree = tree)
    #convert labels B to node numbers
    node_B <- get_node_from_label(label = label_B,
                                  tree = tree)
  }

  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  #Call C++ function
  linSim <- get_lin(node1 = node_A,
                       node2 = node_B,
                    tree_nodes = as.integer(c(tree$edge[,2], root_node)),
                    tree_parents = as.integer(c(tree$edge[,1], 0)))
  return(linSim)
}

#' Jiang and Conrath similarity
#'
#' This determines the Jiang and Conrath similarity for two input nodes in a
#' given tree. This uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}.
#'
#' @param tree A phylo object representing a rooted tree, with an information
#'   content attribute IC.
#' @param label_A The first node label.
#' @param label_B The second node label.
#' @param node_A Alternate parameter, the first node number.
#' @param node_B Alternate parameter, the second node number.
#' @return The Jiang and Conrath similarity in the given tree of the pair of
#'   nodes.
#' @export
#'
#' @references
#' \insertRef{seco2004intrinsic}{treecompareR}
#'
#' \insertRef{jiang1997semantic}{treecompareR}
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' general_JiangConrath_similarity(tree = tree, label_A = 't2', label_B = 't4')
#' general_JiangConrath_similarity(tree = tree, label_A = 'n3', label_B = 't8')
#'
general_JiangConrath_similarity <- function(tree, label_A = NULL, label_B = NULL, node_A = NULL, node_B = NULL){
  if (is.null(label_A) | is.null(label_B)){
    if (is.null(node_A) | is.null(node_B)){
      stop('Please input a pair of labels or a pair of node numbers')
    }else{
      #ensure they are integers
      node_A <- as.integer(node_A)
      node_B <- as.integer(node_B)
    }
  } else {
    #convert labels A to node ID numbers
    node_A <- get_node_from_label(label = label_A,
                                  tree = tree)
    #convert labels B to node numbers
    node_B <- get_node_from_label(label = label_B,
                                  tree = tree)
  }

  root_node <- setdiff(tree$edge[,1], tree$edge[,2])
  #Call C++ function
  jcSim <- get_jiang_conrath(node1 = node_A,
                    node2 = node_B,
                    tree_nodes = as.integer(c(tree$edge[,2], root_node)),
                    tree_parents = as.integer(c(tree$edge[,1], 0)))
}

#' Similarity matrix generator
#'
#' This function takes in a tree and a similarity measure function, and
#' generates a similarity matrix for the given tree and similarity measure. The
#' matrix is symmetric, with rows and columns given by the tip and internal node
#' labels in the ordering given by the tree.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param similarity A similarity measure function that requires as inputs a
#'   tree, label_A, and label_B
#' @return A similarity matrix, which is a symmetric, with values in
#'   \eqn{[0,1]}.
#' @export
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' generate_similarity_matrix(tree = tree, similarity = general_Jaccard_similarity)
#' generate_similarity_matrix(tree = tree, similarity = general_Resnik_similarity)
#'
generate_similarity_matrix <- function(tree, similarity = NULL){
  ifelse(is.null(tree$IC), tree_copy <- attach_information_content(tree), tree_copy <- tree)
  Nnode = length(tree$node.label)
  if (Nnode == 1){
    tree_labels <- tree$tip.label
  } else {
    tree_labels <- c(tree$tip.label, tree$node.label[2:Nnode])
  }


  N <- length(tree_labels)

  sim_matrix <- matrix(nrow = N, ncol = N)
  rownames(sim_matrix) <- tree_labels
  colnames(sim_matrix) <- tree_labels

  #similarity matrix is symmetric so we can save time
  #fill in upper triangular part only
  for (i in 1:N){ #iterate over all rows
    for (j in i:N){ #but only iterate over columns equal to rownum or above
      sim_matrix[i,j] <- similarity(tree = tree_copy,
                                    label_A = tree_labels[i],
                                    label_B = tree_labels[j])
      #assign the symmetric part
      sim_matrix[j,i] <- sim_matrix[i,j]
    }
  }

  return(sim_matrix)

}

check_similarity_inputs <- function(tree = NULL, label_1 = NULL, label_2 = NULL, node_1 = NULL, node_2 = NULL){
  if (is.null(tree) | !('phylo' %in% class(tree))){
    stop('Please input a `phylo` object for the tree parameter!')
  }

  tree_labels <- c(tree$tip.label, tree$node.label)

  if (is.null(label_1)){
    if (is.null(node_1) | !is.numeric(node_1)){
      stop('Please input either a label for `label_1` or a node number for `node_1`!')
    } else {
      node1 <- as.integer(node_1)
    }
  } else {
    node1 <- which(tree_labels %in% label_1)
    if (length(node1) != 1){
      stop('Please input a single node label for `label_1`!')
    }
  }

  if (is.null(label_2)){
    if (is.null(node_2) | !is.numeric(node_2)){
      stop('Please input either a label for `label_2` or a node number for `node_2`!')
    } else {
      node2 <- as.integer(node_2)
    }
  } else {
    node2 <- which(tree_labels %in% label_2)
    if (length(node2) != 1){
      stop('Please input a single node label for `label_2`!')
    }
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

jaccard_similarity <- function(tree = NULL, label_1 = NULL, label_2 = NULL,
                              node_1 = NULL, node_2 = NULL){

  nodes <- check_similarity_inputs(tree = tree, label_1 = label_1, label_2 = label_2,
                                   node_1 = node_1, node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1

  # Handle the case where the root is both input nodes.
  if (all(c(node1, node2) == root)){
    return(1)
  }

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  #information_content <- generate_descendants(tree)
  #information_content$IC <- generate_information_content(tree)
  #information_content <- as.matrix(information_content)

  jaccard <- get_jaccard(node1 = node1, node2 = node2, tree_nodes = tree_nodes,
                       tree_parents = tree_parents)
  return(jaccard)

}


resnik_similarity <- function(tree = NULL, label_1 = NULL, label_2 = NULL,
                              node_1 = NULL, node_2 = NULL){

  nodes <- check_similarity_inputs(tree = tree, label_1 = label_1, label_2 = label_2,
                                   node_1 = node_1, node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- generate_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  resnik <- get_resnik(node1 = node1, node2 = node2, tree_nodes = tree_nodes,
                       tree_parents = tree_parents, information_content = information_content)
  return(resnik)

}

lin_similarity <- function(tree = NULL, label_1 = NULL, label_2 = NULL,
                              node_1 = NULL, node_2 = NULL){

  nodes <- check_similarity_inputs(tree = tree, label_1 = label_1, label_2 = label_2,
                                   node_1 = node_1, node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1

  # Handle the case where the root is both input nodes.
  if (all(c(node1, node2) == root)){
    return(1)
  }

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- generate_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  lin <- get_lin(node1 = node1, node2 = node2, tree_nodes = tree_nodes,
                       tree_parents = tree_parents, information_content = information_content)
  return(lin)

}

jiang_conrath_similarity <- function(tree = NULL, label_1 = NULL, label_2 = NULL,
                              node_1 = NULL, node_2 = NULL){

  nodes <- check_similarity_inputs(tree = tree, label_1 = label_1, label_2 = label_2,
                                   node_1 = node_1, node_2 = node_2)
  node1 <- nodes[[1]]
  node2 <- nodes[[2]]

  root <- length(tree$tip.label) + 1


  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- generate_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  jiang_conrath <- get_jiang_conrath(node1 = node1, node2 = node2, tree_nodes = tree_nodes,
                       tree_parents = tree_parents, information_content = information_content)
  return(jiang_conrath)

}


similarity_matrix <- function(labels_1 = NULL, labels_2 = NULL, nodes_1 = NULL,
                              nodes_2 = NULL, tree = NULL, sim_metric = NA_integer_,
                              upper_tri = TRUE){
  if (is.null(tree) | !('phylo' %in% class(tree))){
    stop('Please input a `phylo` object for the tree parameter!')
  }

  tree_labels <- c(tree$tip.label, tree$node.label)

  if (!is.null(labels_1) & !is.null(labels_2)){
    labels_1 <- labels_1[!is.na(labels_1)]
    labels_2 <- labels_2[!is.na(labels_2)]

    nodes1 <- which(tree_labels %in% labels_1)
    nodes2 <- which(tree_labels %in% labels_2)


  } else if (!is.null(nodes_1) & !is.null(nodes_2)){
    nodes_1 <- nodes_1[!is.na(nodes_1)]
    nodes_2 <- nodes_2[!is.na(nodes_2)]

    nodes1 <- nodes_1[((nodes_1 >= 1) & (nodes_1 <= length(tree_labels)))]
    nodes2 <- nodes_2[((nodes_2 >= 1) & (nodes_2 <= length(tree_labels)))]

    nodes1 <- as.integer(nodes1)
    nodes2 <- as.integer(nodes2)
  } else {
    stop('Please input valid lists for labels_1 and labels_2 or for nodes_1 and nodes_2!')
  }

  if (is.na(sim_metric) | !(sim_metric %in% 1:4)){
    stop('Please input a valid value for sim_metric parameter!')
  }

  root <- length(tree$tip.label) + 1

  node1 <- nodes1[nodes1 != root]
  node2 <- nodes2[nodes2 != root]

  if (any(c(length(nodes1), length(nodes2)) == 0)){
    stop('Please input valid lists for labels_1 and labels_2 or for nodes_1 and nodes_2!')
  }

  tree_nodes <- c(tree$edge[, 2], root)
  tree_parents <- c(tree$edge[, 1], -1)

  information_content <- generate_descendants(tree)
  information_content$IC <- generate_information_content(tree)
  information_content <- as.matrix(information_content)

  similarity <- get_similarity(nodes1 = nodes1, nodes2 = nodes2,
                               tree_nodes = tree_nodes, tree_parents = tree_parents,
                               sim_metric = sim_metric,
                               information_content = information_content,
                               upper_tri = upper_tri)

  return(similarity)

}

#' Random subtree similarity simulation
#'
#' This function takes in a tree, two data sets (or indices representing two
#' subtrees of the tree) and creates random subtrees, calculates their
#' similarity values, and returns the results in a data.frame. This function
#' allows for comparison of the similarity of the input data sets (or subtrees)
#' with similarity of random subtrees of specified sizes.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param data_1 A data.table of chemicals with classifications.
#' @param data_2 A data.table of chemicals with classifications.
#' @param data_1_indices Alternate parameter giving indices of nodes of a
#'   subtree of `tree`.
#' @param data_2_indices Alternate parameter giving indices of nodes of a
#'   subtree of `tree`.
#' @param name_1 An alternate parameter for the name of `data_1`.
#' @param name_2 An alternate parameter for the name of `data_2`.
#' @param label_number Number of labels to use to build the simulated trees.
#' @param repetition Number of simulated trees to build.
#' @param seed Alternate parameter to allow for replication of results.
#' @param only_tips Alternate parameter restricting starting labels to tips or
#'   to tips and internal nodes.
#' @param Jaccard The Jaccard similarity matrix for `tree`.
#' @param Resnik The Resnik similarity matrix for `tree`.
#' @param Lin The Lin similarity matrix for `tree`.
#' @param JiangConrath The JiangConrath similarity matrix for `tree`.
#' @return A data.frame with similarity values for each simulation. The
#'   similarity values reported in each row is the mean similarity value for the
#'   corresponding data set/simulated tree given by the column.
#' @export
#' @importFrom data.table is.data.table
#' @importFrom phangorn Ancestors
#'
#' @examples
#' \donttest{dt1 <- classify_datatable(data.table::data.table(chemical_list_biosolids_2022_05_10)[1:10,])
#' dt1 <- classify_by_smiles(dt1)
#'
#' dt2 <- classify_datatable(data.table::data.table(chemical_list_USGSWATER_2022_05_17)[1:10,])
#' dt2 <- classify_by_smiles(dt2)
#'
#' MonteCarlo_similarity(tree = treecompareR:::chemont_tree, data_1 = dt1, data_2 = dt2,
#'                      name_1 = 'Biosolids', name_2 = 'USGS', seed = 42L, Jaccard = chemont_jaccard,
#'                      Resnik = chemont_resnik_IC_SVH, Lin = chemont_lin_IC_SVH,
#'                      JiangConrath = chemont_jiangconrath_IC_SVH)
#' MonteCarlo_similarity(tree = treecompareR:::chemont_tree, data_1 = dt1, data_2 = dt2,
#'                       name_1 = 'Biosolids', name_2 = 'USGS', label_number = 200, seed = 42L,
#'                       Jaccard = chemont_jaccard,Resnik = chemont_resnik_IC_SVH,
#'                       Lin = chemont_lin_IC_SVH, JiangConrath = chemont_jiangconrath_IC_SVH)}
#'
MonteCarlo_similarity <- function(tree, data_1 = NULL, data_2 = NULL, data_1_indices = NULL, data_2_indices = NULL,
                                  name_1 = 'data_set_1', name_2 =  'data_set_2', label_number = 100,
                                  repetition = 10, seed = NA_real_, only_tips = FALSE, Jaccard = NULL,
                                  Resnik = NULL, Lin = NULL, JiangConrath = NULL){
  if (!is.na(seed) & is.integer(seed)){
    set.seed(seed)
  }
  if(only_tips & (label_number > length(tree$tip.label))){
    stop('Please input a label_number less than the number of tips!')
  } else if (label_number > (length(tree$tip.label) + length(tree$node.label))) {
    stop('Please input a label_number less than the number of nodes and tips!')
  }

  Nnode <- length(tree$node.label)
  dimnames <- c(tree$tip.label, tree$node.label[2:Nnode])

  if (is.null(data_1)){
    if (is.null(data_1_indices)){
      stop('Please input either indices for `data_1_indices` or data.table object for `data_1` parameter!')
    } else {
      temp_labels <- c(tree$tip.label, tree$node.label)[data_1_indices]
      dataset_1_indices <- which(dimnames %in% temp_labels)
    }
  } else if (!data.table::is.data.table(data_1)){
    stop('The `data_1` parameter only takes in a data.table!')
  } else {
    dataset_1_labels <- unlist(get_terminal_labels(data = data_1))
    dataset_1_indices <- which(dimnames %in% dataset_1_labels)
  }

  if (is.null(data_2)){
    if (is.null(data_2_indices)){
      stop('Please input either indices for `data_2_indices` or data.table object for `data_2` parameter!')
    } else {
      temp_labels <- c(tree$tip.label, tree$node.label)[data_2_indices]
      dataset_2_indices <- which(dimnames %in% temp_labels)
    }
  } else if (!data.table::is.data.table(data_2)){
    stop('The `data_2` parameter only takes in a data.table!')
  } else {
    dataset_2_labels <- unlist(get_terminal_labels(data = data_2))
    dataset_2_indices <- which(dimnames %in% dataset_2_labels)
  }


#  if (!(data.table::is.data.table(data_1) & data.table::is.data.table(data_2))){
#    if (is.null(data_1_indices) | is.null(data_2_indices)){
#      stop('Please input either indices for `data_1_indices` and `data_2_indices` or data.table objects for `data_1` and `data_2`')
#    } else {
#      stop('Please input a data.table object for each of the `data_1` and `data_2` parameters!')
#    }
#      }

  if (is.null(Jaccard) & is.null(Resnik) & is.null(Lin) & is.null(JiangConrath)){
    stop('Please input a similarity matrix for at least one of Jaccard, Resnik, Lin, and JiangConrath parameters!')
  }

  simulation_dataframe <- data.frame(Jaccard_all = double(),
                                     Resnik_all = double(),
                                     Lin_all = double(),
                                     JiangConrath_all = double(),
                                     Jaccard_tip_all_data_set_1 = double(),
                                     Resnik_tip_all_data_set_1 = double(),
                                     Lin_tip_all_data_set_1 = double(),
                                     JiangConrath_tip_all_data_set_1 = double(),
                                     Jaccard_tip_all_data_set_2 = double(),
                                     Resnik_tip_all_data_set_2 = double(),
                                     Lin_tip_all_data_set_2 = double(),
                                     JiangConrath_tip_all_data_set_2 = double(),
                                     Jaccard_all_data_set_1 = double(),
                                     Resnik_all_data_set_1 = double(),
                                     Lin_all_data_set_1 = double(),
                                     JiangConrath_all_data_set_1 = double(),
                                     Jaccard_all_data_set_2 = double(),
                                     Resnik_all_data_set_2 = double(),
                                     Lin_all_data_set_2 = double(),
                                     JiangConrath_all_data_set_2 = double(),
                                     all_nodes = integer(),
                                     all_tips = integer())



#  if (is.data.table(data_1) & is.data.table(data_2)) {
#    dataset_1_labels <- unlist(get_terminal_labels(data = data_1))
#    dataset_1_indices <- which(dimnames %in% dataset_1_labels)
#
#    dataset_2_labels <- unlist(get_terminal_labels(data = data_2))
#    dataset_2_indices <- which(dimnames %in% dataset_2_labels)
#  } else {
#    dataset_1_labels <- dimnames[data_1_indices]
#    dataset_1_indices <- data_1_indices
#
#    dataset_2_labels <- dimnames[data_2_indices]
#    dataset_2_labels <- data_2_indices
#  }

  get_indices <- function(indices_1, indices_2) {
    dat <- expand.grid(indices_1, indices_2)
    as.matrix(unique(cbind(pmin(dat[, 1], dat[, 2]),
                           pmax(dat[, 1], dat[,2]))),
              ncol = 2)
  }


  for (i in 1:repetition){
    if (only_tips) {
      label_start <- sample(tree$tip.label, label_number)
    } else {
      label_start <- sample(dimnames, label_number)
    }

    label_start_nodes <- sapply(label_start, function(t) {
      match(t, dimnames)
    })

    label_branches <- sapply(label_start_nodes, function(t) {
      if (!is.na(t)) {
        phangorn::Ancestors(t, x = tree)
      }
    })
    label_ancestors <-unique(unlist(label_branches))

    label_all_nodes <- dimnames[union(label_start_nodes, unlist(label_branches))]
    label_all_tips <- label_start[which(label_start %in% tree$tip.label)]

    all_node_indices <- which(dimnames %in% label_all_nodes)
    all_tip_indices <- which(dimnames %in% label_all_tips)

    #print(all_node_indices)
    #print(all_tip_indices)



    all_nodes <- get_indices(all_node_indices, all_node_indices)
    all_tip_all_dataset_1 <- get_indices(all_tip_indices, dataset_1_indices)
    all_tip_all_dataset_2 <- get_indices(all_tip_indices, dataset_2_indices)
    all_node_all_dataset_1 <- get_indices(all_node_indices, dataset_1_indices)
    all_node_all_dataset_2 <- get_indices(all_node_indices, dataset_2_indices)

    print(i)

    new_row <- double(22L)

    if (!is.null(Jaccard)){
      new_row[[1]] <- mean(Jaccard[all_nodes])#mean(Jaccard[all_node_indices, all_node_indices][upper.tri(Jaccard[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[5]] <- mean(Jaccard[all_tip_all_dataset_1])#mean(Jaccard[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)][upper.tri(Jaccard[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)], diag = TRUE)])
      new_row[[9]] <- mean(Jaccard[all_tip_all_dataset_2])#mean(Jaccard[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)][upper.tri(Jaccard[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)], diag = TRUE)])
      new_row[[13]] <- mean(Jaccard[all_node_all_dataset_1])#mean(Jaccard[union(all_node_indices, dataset_1_indices), union(all_node_indices, dataset_1_indices)][upper.tri(Jaccard[union(all_node_indices, dataset_1_indices), union(all_node_indices, dataset_1_indices)], diag = TRUE)])#mean(Jaccard[all_node_indices, dataset_1_indices])
      new_row[[17]] <- mean(Jaccard[all_node_all_dataset_1])#mean(Jaccard[union(all_node_indices, dataset_2_indices), union(all_node_indices, dataset_2_indices)][upper.tri(Jaccard[union(all_node_indices, dataset_2_indices), union(all_node_indices, dataset_2_indices)], diag = TRUE)])
        #mean(Jaccard[all_node_indices, dataset_2_indices])
    }

    if (!is.null(Resnik)){
      new_row[[2]] <- mean(Resnik[all_nodes])#mean(Resnik[all_node_indices, all_node_indices][upper.tri(Resnik[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[6]] <- mean(Resnik[all_tip_all_dataset_1])#mean(Resnik[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)][upper.tri(Resnik[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)], diag = TRUE)])
      new_row[[10]] <- mean(Resnik[all_tip_all_dataset_2])#mean(Resnik[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)][upper.tri(Resnik[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)], diag = TRUE)])
      new_row[[14]] <- mean(Resnik[all_node_all_dataset_1])#mean(Resnik[union(all_node_indices, dataset_1_indices), union(all_node_indices, dataset_1_indices)][upper.tri(Resnik[union(all_node_indices, dataset_1_indices), union(all_node_indices, dataset_1_indices)], diag = TRUE)])#mean(Resnik[all_node_indices, dataset_1_indices])
      new_row[[18]] <- mean(Resnik[all_node_all_dataset_2])#mean(Resnik[union(all_node_indices, dataset_2_indices), union(all_node_indices, dataset_2_indices)][upper.tri(Resnik[union(all_node_indices, dataset_2_indices), union(all_node_indices, dataset_2_indices)], diag = TRUE)])#mean(Resnik[all_node_indices, dataset_2_indices])
    }

    if (!is.null(Lin)){
      new_row[[3]] <- mean(Lin[all_nodes])#mean(Lin[all_node_indices, all_node_indices][upper.tri(Lin[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[7]] <- mean(Lin[all_tip_all_dataset_1])#mean(Lin[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)][upper.tri(Lin[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)], diag = TRUE)])
      new_row[[11]] <- mean(Lin[all_tip_all_dataset_1])#mean(Lin[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)][upper.tri(Lin[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)], diag = TRUE)])
      new_row[[15]] <- mean(Lin[all_node_all_dataset_1])#mean(Lin[union(all_node_indices, dataset_1_indices), union(all_node_indices, dataset_1_indices)][upper.tri(Lin[union(all_node_indices, dataset_1_indices), union(all_node_indices, dataset_1_indices)], diag = TRUE)])#mean(Lin[all_node_indices, dataset_1_indices])
      new_row[[19]] <- mean(Lin[all_node_all_dataset_2])#mean(Lin[union(all_node_indices, dataset_2_indices), union(all_node_indices, dataset_2_indices)][upper.tri(Lin[union(all_node_indices, dataset_2_indices), union(all_node_indices, dataset_2_indices)], diag = TRUE)])#mean(Lin[all_node_indices, dataset_2_indices])
    }

    if (!is.null(JiangConrath)){
      new_row[[4]] <- mean(JiangConrath[all_nodes])#mean(JiangConrath[all_node_indices, all_node_indices][upper.tri(JiangConrath[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[8]] <- mean(JiangConrath[all_tip_all_dataset_1])#mean(JiangConrath[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)][upper.tri(JiangConrath[union(all_tip_indices, dataset_1_indices), union(all_tip_indices, dataset_1_indices)], diag = TRUE)])
      new_row[[12]] <- mean(JiangConrath[all_tip_all_dataset_1])#mean(JiangConrath[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)][upper.tri(JiangConrath[union(all_tip_indices, dataset_2_indices), union(all_tip_indices, dataset_2_indices)], diag = TRUE)])
      new_row[[16]] <- mean(JiangConrath[all_node_all_dataset_1])#mean(JiangConrath[union(all_node_indices, dataset_1_indices), union(all_node_indices, dataset_1_indices)][upper.tri(JiangConrath[union(all_node_indices, dataset_1_indices),union(all_node_indices, dataset_1_indices)], diag = TRUE)])
      new_row[[20]] <- mean(JiangConrath[all_node_all_dataset_2])#mean(JiangConrath[union(all_node_indices, dataset_2_indices), union(all_node_indices, dataset_2_indices)][upper.tri(JiangConrath[union(all_node_indices, dataset_2_indices),union(all_node_indices, dataset_2_indices)], diag = TRUE)])
        #mean(JiangConrath[union(all_node_indices, dataset_2_indices),union(all_node_indices, dataset_2_indices)])
    }

    new_row[[21]] <- length(all_node_indices)
    new_row[[22]] <- length(all_tip_indices)

    simulation_dataframe[i, ] <- new_row
  }

  names(simulation_dataframe)[5:8] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), 'tip', name_1, 'all')
  names(simulation_dataframe)[9:12] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), 'tip', name_2, 'all')
  names(simulation_dataframe)[13:16] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), 'all', name_1, 'all')
  names(simulation_dataframe)[17:20] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), 'all', name_2, 'all')


  return(simulation_dataframe)
}

#' Similarity cutoffs
#'
#' This function gives cutoffs for similarity values and percent of
#' representation of a data set. The function determines what percentage of a
#' data set that is represented by the induced subtree of the data set for
#' various values of a fixed similarity measure.
#'
#' @param mat A similarity matrix corresponding to a similarity measure and a
#'   rooted tree .
#' @param data A data.table of chemicals with classifications.
#' @param tax_level_labels Parameter giving classification levels.
#' @param neighbors A parameter giving how many neighbors to use for finding
#'   label average values.
#' @param cutoff An alternate parameter giving the cutoff percentage value.
#' @param labels An alternate parameter giving a list of node labels
#'   corresponding to a subtree of a rooted tree..
#' @param counts An alternate parameter giving the counts of occurrence for each
#'   label.
#' @return Named list of percentage of data represented by similarity values.
#'   The names are the similarity values. The values of the list are percentages
#'   of data represented by allowing similarity values equal to the names.
#' @export
#'
#' @examples
#' \donttest{dt <- classify_datatable(data.table::data.table(chemical_list_biosolids_2022_05_10)[1:10,])
#' dt <- classify_by_smiles(dt)
#'
#' get_cutoffs(mat = chemont_jaccard, data = dt)
#' get_cutoffs(mat = chemont_jaccard, data = dt, neighbors = 6)}
#'
get_cutoffs <- function(mat, data, tax_level_labels = NULL, neighbors = 3, cutoff = NA_real_, labels = NULL, counts = NULL){
  if (is.data.table(data)){
    if (is.null(tax_level_labels)){
      tax_level_labels <- c('kingdom', 'superclass', 'class', 'subclass',
                      'level5', 'level6', 'level7', 'level8',
                      'level9', 'level10', 'level11')
    }
    counts <- get_number_of_labels(data = data, tax_level_labels = tax_level_labels)
    labels <- names(counts)

  }

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

  average_val <- unname(apply(temp_mat, MARGIN = 1, function(t) {sum(sort(t, decreasing = TRUE)[1:neighbors])/neighbors}))

  total = sum(counts)

  temp_counts <- counts[order(average_val, decreasing = TRUE)]

  unique_avgs <- sort(unique(average_val))


  #margin <- min(unique_avgs[2:length(unique_avgs)] - unique_avgs[1:(length(unique_avgs)-1)])/3

  percentages <- sapply(rev(unique_avgs), function(t) {sum(temp_counts[sort(average_val, decreasing = TRUE) >= t])/total})

  names(percentages) <- rev(unique_avgs)

  if (!is.na(cutoff)) {
    # returns maximum value that achieves cutoff percentage threshold
    return(names(percentages)[[min(which(percentages >= cutoff))]])
  } else {
    return(percentages)
  }
  }

#' Drop tips and nodes
#'
#' This function is an extension of drop.tip from the APE package. It takes in a
#' data set, collects the labels from classifications of the chemicals in the
#' data set, and drops tips and nodes until the induced subtree that remains
#' consists solely of the labels associated to the data set and their ancestors.
#' Alternatively, one can provide a set of labels instead and the induced
#' subtree is constructed from these labels in the same manner.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param data A data.table of chemicals with classifications.
#' @param labels An alternate parameter for a set of labels of the subtree.
#' @param tax_level_labels An alternate parameter passed to the
#'   \code{\link{get_terminal_labels}} function.
#' @return A phylo object representing the induced subtree of the data.
#'
#' @importFrom ape drop.tip
#'
#' @references
#' \insertRef{apepackage}{treecompareR}

drop_tips_nodes <- function(tree,
                            data = NULL,
                            labels = NULL,
                            nodes = NULL,
                            level = NULL,
                            keep_descendants = NULL,
                            tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                                 'level5', 'level6', 'level7', 'level8',
                                                 'level9', 'level10', 'level11')){
  if (!is.null(data)){
    #get terminal labels for each item in this data set
    #these are the labels to keep
    if(is.null(keep_descendants)){
      keep_descendants = FALSE
    }
    tip_node_labels <- get_terminal_labels(data = data,
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


#' Adjust branch lengths
#'
#' This is a helper function that is used to reset branch lengths for pruned
#' trees.
#'
#' @param tree A phylo object representing a rooted tree.
#' @return A tree with adjust branch lengths.
#'

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
      #next_level <- c(next_level, children)
      #children_height <- tree_height[children]

      #children_tip <- which(children_height == 1)
      #if (length(children_tip) > 0){
      #  plot_height[which(names(plot_height) %in% children[children_tip])] <- 0
      #  children <- children[-children_tip]
      #}

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


#' Compare similarity measures
#'
#' This functions compares the similarity measures of Jaccard, Resnik, Lin, and
#' Jiang and Conrath on a variety of trees. The tree types include the
#' caterpillar, star, and balanced trees. The input parameter `n` indicates how
#' many tips for the caterpillar and balanced tree. The star has 2n tips, and
#' all trees have 2n-1 total nodes and tips.
#'
#' @param n Each tree has 2n-1 total nodes and tips.
#' @return A data.frame consisting of the mean self-similarity scores for each
#'   tree and similarity measure.
compare_similarity_measures <- function(n){
  caterpillar <- generate_caterpillar(n)
  star <- generate_star(2*n)
  balanced <- generate_balanced(n)

  cat_IC <- attach_information_content(caterpillar)
  star_IC <- attach_information_content(star)
  balanced_IC <- attach_information_content(balanced)

  cat_Jaccard <- generate_similarity_matrix(caterpillar, similarity = general_Jaccard_similarity)
  star_Jaccard <- generate_similarity_matrix(star, similarity = general_Jaccard_similarity)
  balanced_Jaccard <- generate_similarity_matrix(balanced, similarity = general_Jaccard_similarity)

  cat_Resnik <- generate_similarity_matrix(cat_IC, similarity = general_Resnik_similarity)
  star_Resnik <- generate_similarity_matrix(star_IC, similarity = general_Resnik_similarity)
  balanced_Resnik <- generate_similarity_matrix(balanced_IC, similarity = general_Resnik_similarity)

  cat_Lin <- generate_similarity_matrix(cat_IC, similarity = general_Lin_similarity)
  star_Lin <- generate_similarity_matrix(star_IC, similarity = general_Lin_similarity)
  balanced_Lin <- generate_similarity_matrix(balanced_IC, similarity = general_Lin_similarity)

  cat_JiangConrath <- generate_similarity_matrix(cat_IC, similarity = general_JiangConrath_similarity)
  star_JiangConrath <- generate_similarity_matrix(star_IC, similarity = general_JiangConrath_similarity)
  balanced_JiangConrath <- generate_similarity_matrix(balanced_IC, similarity = general_JiangConrath_similarity)



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

#'Get node number of clade (ancestor at a specified level)
#'
#'Helper function to get the node number defining a clade for a specified input
#'node number
#'@param node The node number(s) for which to get the clade(s)
#'@param tree The underlying tree (node numbers refer to this tree)
#'@param level The hierarchical taxonomy level at which to get the clade(s). Root
#'  is level 0. Default value is 2 (superclass level, in ChemOnt).
get_clade <- function(node,
                      tree,
                      level){
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

#' List all clades in a tree at a specified level
#'
#' @param tree The \code{\link[ape]{phylo}}-class tree object
#' @param level The level at which to display nodes (0 is the root)
#' @return A data.frame with four variables: \code{node} (the node number in the
#'   tree); \code{level} (the level of the node in the tree, where root is level
#'   0); \code{parent} (the node number of the node's immediate parent); and
#'   \code{Name} (the text label of the node).
#' @export
get_all_clades <- function(tree, level){
tree_df <- get_tree_df(tree = tree)
clade_df <- tree_df[tree_df$level %in% level, ]
return(clade_df)
}

#' Bind individual entities as new tips to a tree
#'
#' @param tree The base tree as a \code{phylo}-class object. Tips will be bound to this tree.
#' @param data Either one data.frame, or a list of data.frames, containing
#'   classified entities. Each row of the data.frame is one entity. The
#'   data.frames must include the column names specifeid in
#'   \code{tax_level_labels} and \code{entity_id_col}.
#' @param entity_id_col The column name in \code{data} containing identifying
#'   labels for the entities.
#' @param tax_level_labels Taxonomy levels used for classification in \code{data}. Default is
#'   the Chemont taxonomy levels: \code{c('kingdom', 'superclass', 'class',
#'   'subclass','level5', 'level6', 'level7', 'level8','level9', 'level10',
#'   'level11')}.
#' @return A \code{phylo}-class object.
#'
bind_entities <- function(tree,
                          data,
                          entity_id_col,
                          tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                               'level5', 'level6', 'level7', 'level8',
                                               'level9', 'level10', 'level11')){

  #if a list of data frames is provided, rowbind it all together
  #this will be the "master list" of entities
  if(!is.data.frame(data)){
  if(is.list(data) &
     all(sapply(data, is.data.frame))){
    data <- dplyr::bind_rows(data)
  }
  }

  #if terminal label not already in data, add it
  if(!"terminal_label" %in% names(data)){
    data <- add_terminal_label(data = data,
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
                                 Name = tmpdf[[entity_id_col]])
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

#'Prune a tree
#'
#'Prune a tree to keep only a specified subtree.
#'
#' # How to specify the subtree to keep
#'\code{prune_to} defines the subtree to *keep* (everything else will be pruned
#'away). It may be specified in several different ways.
#'
#' ## As a data.frame
#'
#'If \code{prune_to} is a \code{data.frame} of classified entities, it must have
#'columns corresponding to, and named for, each of the taxonomy levels as
#'defined in the argument \code{tax_level_labels}, containing the labels at the
#'corresponding level for each entity. It must also have at least one more
#'column, uniquely identifying the entities; the name of the additional column
#'does not matter, as long as it is not the same as one of the taxonomy levels.
#'The result will be to keep only the subtree induced by this classified data
#'set, i.e., only the branches of the tree that occur in this classified data
#'set. By default, any descendants of the node labels in the \code{data.frame}
#'that do not themselves appear in the \code{data.frame} will *not* be kept. If
#'you want to keep descendants that do not themselves appear in the
#'\code{data.frame}, specify \code{keep_descendants = TRUE}.
#'
#' ## As the name of a taxonomy level
#'
#'If \code{prune_to} is the name of a taxonomy level (one of the levels defined
#'in argument \code{tax_level_labels}), the result will be to keep only nodes at
#'that taxonomic level or less-specific levels. (For example, for the ChemOnt
#'taxonomy, specifying \code{prune_to = "class"} will keep only nodes at levels
#'"kingdom", "superclass", and "class". Any nodes at level "subclass", "level5",
#'"level6", ... "level11" will be dropped. (If you specify \code{prune_to} as
#'the name of a taxonomy level, and also specify \code{keep_descendants = TRUE},
#'the result will be to keep the whole tree.)
#'
#' ## As a vector of node/tip labels or numbers
#'
#'If \code{prune_to} is a vector of node/tip labels (i.e., labels appearing in
#'\code{tree$node.label} and/or \code{tree$tip.label}) or node/tip numbers (i.e.
#'node/tip index numbers between 1 and \code{ape::Ntip(tree) +
#'ape::Nnode(tree)}), the result will be to keep only the nodes/tips that are in
#'that vector, keep their common ancestors, and (by default) also keep their
#'descendants if any. The intention of keeping the descendants by default is to
#'allow the user to prune to specified clades simply by specifying the labels or
#'node numbers of the MRCAs of the clades. For example, using the ChemOnt
#'taxonomy, you could prune to keep all branches in the superclass
#'"Organohalogen compounds" by simply specifying \code{prune_to = "Organohalogen
#'compounds"}.  If you do *not* wish to keep the descendants of the specified
#'node labels/numbers, then specify \code{keep_descendants = FALSE}.
#'
#'@param tree The tree to be pruned, as a \code{\link[ape]{phylo}}-class object.
#'@param prune_to What to *keep* from the base tree (everything else will be
#'  pruned away). May be a \code{data.frame} of classified data; one or more
#'  labels in the tree (tip or internal node labels); one or more node numbers
#'  in the tree (tip or internal nodes); or the name of a taxonomy level (one of
#'  the items in \code{tax_level_labels}). Default is NULL, which results in no
#'  pruning being done (i.e., the base tree is returned as-is). See Details.
#'@param keep_descendants Whether to keep descendants of what is specified in
#'  \code{prune_to}. Default NULL will choose the behavior based on the class of
#'  \code{prune_to}: when \code{prune_to} is a \code{data.frame} or one of the
#'  taxonomy level labels, \code{keep_descendants = FALSE} by default. When
#'  \code{prune_to} is a vector of node/tip labels or numbers in the base
#'  tree,\code{keep_descendants = TRUE} by default. See Details.
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
                       keep_descendants = NULL,
                       adjust_branch_length = FALSE,
                       tax_level_labels = chemont_tax_levels){
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

#' Convert a phylo tree into a wide-format "classified" data.frame

as_classified.phylo <- function(tree,
                                tax_level_labels = c('kingdom', 'superclass', 'class', 'subclass',
                                                                      'level5', 'level6', 'level7', 'level8',
                                                                      'level9', 'level10', 'level11')){

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

#' Calculate similarity measures for two datasets
#' @param data_1 A data.frame of classified entities
#' @param data_2 Another data.frame of classified entities
#' @param terminal_label The variable name in the two data.frames that denotes
#'   the terminal label of the classification. Default "terminal_label".
#' @param tree The taxonomy tree to use. Default \code{\link{chemont_tree}}.
#' @param tax_level_labels The set of taxonomy levels to use. Default \code{\link{chemont_tax_levels}}.
#' @param similarity The similarity metric to calculate. Default "jaccard".
calc_similarity_data <- function(data_1,
                            data_2,
                            terminal_label = "terminal_label",
                            tree = chemont_tree,
                            tax_level_labels = chemont_tax_levels,
                            similarity = "jaccard"){

  if(similarity %in% "jaccard"){
    similarity_fun <- "general_Jaccard_similarity"
  }
  #calculate pairwise similarity of ancestry of terminal labels in two data sets

  #check for terminal_label
  if(terminal_label == "terminal_label"){
  if(!(terminal_label %in% names(data_1))){
    data_1 <- add_terminal_label(data = data_1, tax_level_labels = tax_level_labels)
  }

  if(!(terminal_label %in% names(data_2))){
    data_2 <- add_terminal_label(data = data_2, tax_level_labels = tax_level_labels)
  }
  }

  #Keep only data with terminal labels in the tree
  data_1 <- data_1[data_1[[terminal_label]] %in%
                     c(tree$tip.label, tree$node.label), ]
  data_2 <- data_2[data_2[[terminal_label]] %in%
                     c(tree$tip.label, tree$node.label), ]

  #enumerate a matrix of pairs of terminal labels
  #first get all unique labels across both datasets
  #sort them
  mlabs <- sort(union(data_1[[terminal_label]],
                       data_2[[terminal_label]]))
  #keep only the ones that appear in each data set
  mrowlabs <- mlabs[mlabs %in% data_1[[terminal_label]]]
  mcollabs <- mlabs[mlabs %in% data_2[[terminal_label]]]

  m <- matrix(nrow = length(mrowlabs),
              ncol = length(mcollabs))
  rownames(m) <- mrowlabs
  colnames(m) <- mcollabs

  #now calculate matrix elements
  #only need to calc upper triangular part;
  #can then assign

  for(i in 1:(length(mlabs)-1)){
    if(mlabs[i] %in% mrowlabs &
       mlabs[i] %in% mcollabs){
      #entity is 100% similar to itself
      m[mlabs[i], mlabs[i]] <- 1
    }
    for(j in (i+1):(length(mlabs))){
      tmp <- do.call(similarity_fun,
                     list(tree = tree,
                          label_A = mlabs[i],
                          label_B = mlabs[j]))
      if(mlabs[i] %in% mrowlabs &
         mlabs[j] %in% mcollabs){
        m[mlabs[i], mlabs[j]] <- tmp
      }

      if(mlabs[j] %in% mrowlabs &
               mlabs[i] %in% mcollabs){
        m[mlabs[j], mlabs[i]] <- tmp
      }

      if(mlabs[j] %in% mrowlabs &
         mlabs[j] %in% mcollabs){
        #entity is 100% similar to itself
        m[mlabs[j], mlabs[j]] <- 1
      }
    } #end j loop
  } #end i loop

return(m)
}

#' Get subtree node numbers
#'
#' This is a helper function that takes a classified set and a base taxonomy tree and
#' provides all node numbers in the subtree corresponding to the data set.
#'
#' @param data A classified data set.
#' @param base_tree A phylo-class tree object representing the base tree. Default is \code{\link{chemont_tree}},
#'   the full ChemOnt taxonomy tree.
#' @param tax_level_labels A vector of levels for the taxonomy. Default is
#'   \code{\link{chemont_tax_levels}}, the levels of the ChemOnt taxonomy.
#' @return A vector of node numbers in the base tree that are represented in the
#'   subtree corresponding to the data set.

get_subtree_nodes <- function(data,
                              base_tree = chemont_tree,
                              tax_level_labels = chemont_tax_levels){
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
