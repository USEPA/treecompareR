# This generates the number of descendants for each node. The edges are stored
# in a two column data.frame, with the first column indicating the parent and
# the second column indicating the child. Internal numbering of nodes is given
# by tips first and then nodes, with the root given by the first node of the
# internal nodes. The output also includes the number of children for each node
# and the node level.

#' Generates a data.frame of the number of descendants for each node of the input tree
#'
#' @param tree A phylo object
#' @return data.frame consisting of the node number, descendants, children, and level
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

# This function creates a data.frame and records the level of each node,
# starting at the root and exploring each successive level of depth. It takes in
# a rooted tree as an input. It returns a data.frame of length equal to the
# number of nodes of the tree, and depth of each node with values ranging from
# 0 (the root) to the maximum depth of all the nodes, inclusive.

#' Generates a data.frame of node levels for input rooted tree.
#'
#' @param tree A phylo object representing a rooted tree.
#' @return data.frame consisting of the node number, and the level of each node
#' @import ape
get_levels <- function(tree){
  if(!ape::is.rooted(tree))
    stop('Input tree must be rooted!')
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

# This function takes in a rooted-tree 'phylo' object and returns a data.frame
# with information of descendants, children, tree level, and a log-based
# information content for each node. The log-based information content assigns
# to each node the value 1 - log(1+descendants)/log(total nodes), giving tips a
# value of 1, the root a value of 0, and monotonically decreasing in value
# from tip to root.

#' Generates a data.frame of information content
#'
#' @param tree A phylo object representing a rooted tree.
#' @param log_descendants Alternate parameter for specifying type of information content.
#' @return data.frame consisting of node number, children, descendants, level, and information content for each node
#' @export
#' @import ape
generate_information_content <- function(tree, log_descendants = TRUE){
  if (!inherits(tree, 'phylo')){
    stop("Please input an object of 'phylo' class!")
  }

  if (!ape::is.rooted(tree)){
    stop('Please input a rooted tree!')
  }

  if (log_descendants){
    descendants <- generate_descendants(tree)
    n_node <- dim(descendants)[[1]]
    descendants$log_descendants <- 1 - (log(1 + descendants$descendants)/log(n_node))
    return(descendants)
  }
}

# This function checks to see if the rooted tree has an information content attached.
# If not, it then attaches one based on type of information content requested.

#' Attaches information content if missing for input tree
#'
#' @param tree A phylo object representing a rooted tree.
#' @param log_descendants Alternate parameter determining type of information content to use.
#' @return phylo object with information content data.frame attached
#' @export
#' @import ape
attach_information_content <- function(tree, log_descendants = TRUE){
  if (!inherits(tree, 'phylo')){
    stop("Please input an object of 'phylo' class!")
  }

  if (!ape::is.rooted(tree)){
    stop("Please input a rooted tree!")
  }

  if (!is.null(tree$IC)){
    indices <- which(names(tree$IC) %in% c('node', 'descendants', 'children', 'level'))
    if (length(indices) < 4){
      warning('Missing columns in tree$IC!')
      return(tree)
    }
  }

  tree$IC <- generate_information_content(tree = tree, log_descendants = log_descendants)

  return(tree)
}

#' Generates a list of ancestors for a given node in a tree.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label The node label.
#' @param node_number Alternate parameter, the number of the given node.
#' @return A list of nodes back to the root of ancestors for the given node.
get_ancestors <- function(tree, label, node_number = NULL){
  if (!is.null(node_number)){
    ifelse(is.numeric(node_number) & (node_number %in% 1:(1 + length(tree$edge))), index <- node_number, stop('Please input a correct value for node_number'))
  } else {
    if (label %in% c(tree$tip.label, tree$node.label)){
      index <- which(c(tree$tip.label, tree$node.label) == label)
    } else {
      stop(paste0('Label `', label, '` belongs neither to a node nor a tip!'))
    }
  }
  ancestor_nodes <- c()
  temp <- tree$edge[tree$edge[, 2] == index, 1]
  while(length(temp) > 0){
    ancestor_nodes <- c(ancestor_nodes, temp)
    temp <- tree$edge[tree$edge[, 2] == temp, 1]}
  return(sapply(ancestor_nodes, function(t) {tree$node.label[[t-length(tree$tip.label)]]}))
}


#' Returns the tree level of the given node
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label The node label.
#' @param node_number Alternate parameter, the number of the given node.
#' @return The level of the node from the root of the tree.
get_tip_level <- function(tree, label, node_number = NULL){
  return(length(get_ancestors(tree = tree, label = label, node_number = node_number)))
}


# This function takes in a tree and two labels and returns the Jaccard distance
# of the nodes corresponding to the labels within the tree.

#' This determines the Jaccard distance for two input node labels in a given tree.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label_A The first label.
#' @param label_B The second label.
#' @return The Jaccard distance of the label sets for the root to node path.
general_Jaccard_dist <- function(tree, label_A, label_B){
  tree_labels <- c(tree$tip.label, tree$node.label)
  index_A <- which(label_A == tree_labels)
  index_B <- which(label_B == tree_labels)

  if ((length(index_A) + length(index_B)) == 0){
    stop('Please input correct labels for Label_A and Label_B')
  } else if ((length(index_A) + length(index_B)) == 1){
    if (length(index_A) == 0){
      stop('Please input correct label for label_A')
    }
    stop('Please input correct label for label_B')
  }

  # Get path labels for label_A and label_B
  ancestors_A <- c(label_A, get_ancestors(tree, label_A))
  ancestors_B <- c(label_B, get_ancestors(tree, label_B))

  # Get node numbers for label_A and label_B paths
  ancestors_indices_A <- which(tree_labels %in% ancestors_A)
  ancestors_indices_B <- which(tree_labels %in% ancestors_B)
  common_ancestors <- intersect(ancestors_A, ancestors_B)

  # Find the node number of the MRCA, it will be in the order of the ancestors_A
  # which starts at the current label and goes back to the root
  index_C <- which(tree_labels == common_ancestors[[1]])
  label_C <- tree_labels[[index_C]]

  level_A <- length(ancestors_A) - 1
  level_B <- length(ancestors_B) - 1
  level_C <- length(common_ancestors) - 1

  # In the case both label_A and label_B are the root
  if ((level_A + level_B) == 0) {return(0)}

  return((level_A + level_B - 2*level_C)/(level_A + level_B - level_C))
}

# This function takes in a tree and two labels and returns the Jaccard
# similarity of the two labels within the tree

#' This determines the Jaccard similarity for two input node labels in a given tree.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label_A The first label.
#' @param label_B The second label.
#' @return The Jaccard distance of the label sets for the root to node path.
general_Jaccard_similarity <- function(tree, label_A, label_B){
  return(1 - general_Jaccard_dist(tree, label_A, label_B))
}

# This function takes in a tree and two labels and returns the Resnik similarity
# of the two labels within the tree.

#' This determines the Resnik similarity for two input nodes in a given tree.
#'
#' @param tree A phylo object representing a rooted tree, with an information content attribute IC.
#' @param label_A The first node label.
#' @param label_B The second node label.
#' @param node_A Alternate parameter, the first node number.
#' @param node_B Alternate parameter, the second node number.
#' @return The Resnik similarity in the given tree of the pair of nodes.
general_Resnik_similarity <- function(tree, label_A = NULL, label_B = NULL, node_A = NULL, node_B = NULL) {
  labels <- c(tree$tip.label, tree$node.label)
  if (is.null(label_A) | is.null(label_B)){
    if (is.null(node_A) | is.null(node_B)){
      stop('Please input a pair of labels or a pair of node numbers')
    } else {
      index1 <- node_A
      index2 <- node_B
    }
  } else {
    index1 <- which(labels == label_A)
    index2 <- which(labels == label_B)
  }

  # Get the common ancestors. This may include either of the nodes corresponding
  # to the input node labels/numbers.
  commonAncestors <- intersect(c(labels[[index1]],unname(get_ancestors(tree, node_number = index1))),
                               c(labels[[index2]], unname(get_ancestors(tree, node_number = index2))))

  MRCA <- which(labels %in% commonAncestors[[1]])
  value <- tree$IC[MRCA, 'log_descendants']

  return(value)
}

# This function takes in a tree and two labels and returns the Lin similarity
# of the two labels within the tree.

#' This determines the Lin similarity for two input nodes in a given tree.
#'
#' @param tree A phylo object representing a rooted tree, with an information content attribute IC.
#' @param label_A The first node label.
#' @param label_B The second node label.
#' @param node_A Alternate parameter, the first node number.
#' @param node_B Alternate parameter, the second node number.
#' @return The Lin similarity in the given tree of the pair of nodes.
general_Lin_similarity <- function(tree, label_A = NULL, label_B = NULL, node_A = NULL, node_B = NULL){
  labels <- c(tree$tip.label, tree$node.label)
  if (is.null(label_A) | is.null(label_B)){
    if (is.null(node_A) | is.null(node_B)){
      stop('Please input a pair of labels or a pair of node numbers')
    } else {
      index_A <- node_A
      index_B <- node_B
    }
  } else {
    index_A <- which(labels == label_A)
    index_B <- which(labels == label_B)
  }

  resSim <- general_Resnik_similarity(tree, node_A = index_A, node_B = index_B)
  denominator <- sum(tree$IC[c(index_A, index_B), 'log_descendants'])
  return(ifelse(denominator == 0, 1, 2*resSim/denominator))
}

# This function takes in a tree and two labels and returns the Jiang-Conrath
# similarity of the two labels within the tree.

#' This determines the Jiang and Conrath similarity for two input nodes in a given tree.
#'
#' @param tree A phylo object representing a rooted tree, with an information content attribute IC.
#' @param label_A The first node label.
#' @param label_B The second node label.
#' @param node_A Alternate parameter, the first node number.
#' @param node_B Alternate parameter, the second node number.
#' @return The Jiang and Conrath similarity in the given tree of the pair of nodes.
general_JiangConrath_similarity <- function(tree, label_A = NULL, label_B = NULL, node_A = NULL, node_B = NULL){
  labels <- c(tree$tip.label, tree$node.label)
  if (is.null(label_A) | is.null(label_B)){
    if (is.null(node_A) | is.null(node_B)){
      stop('Please input a pair of labels or a pair of node numbers')
    } else {
      index_A <- node_A
      index_B <- node_B
    }
  } else {
    index_A <- which(labels == label_A)
    index_B <- which(labels == label_B)
  }

  resSim <- general_Resnik_similarity(tree, node_A = index_A, node_B = index_B)
  informationSum <- sum(tree$IC[c(index_A, index_B), 'log_descendants'])
  return(1 - ((informationSum - 2*resSim)/2))
}

#' Generates a similarity matrix for a given tree and similarity measure
#' @param tree A phylo object representing a rooted tree.
#' @param similarity A similarity measure function that requires as inputs a tree, label_A, and label_B
#' @return A similarity matrix, which is a symmetric matrix with values in [0,1].
#' @export
generate_similarity_matrix <- function(tree, similarity = NULL){
  ifelse(is.null(tree$IC), tree_copy <- attach_information_content(tree), tree_copy <- tree)
  Nnode = length(tree$node.label)
  tree_labels <- c(tree$tip.label, tree$node.label[2:Nnode])
  N <- length(tree_labels)

  sim_matrix <- matrix(nrow = N, ncol = N)
  rownames(sim_matrix) <- tree_labels
  colnames(sim_matrix) <- tree_labels

  for (i in 1:N){
    for (j in i:N){
      sim_matrix[i,j] <- similarity(tree = tree_copy, label_A = tree_labels[i], label_B = tree_labels[j])
      sim_matrix[j,i] <- sim_matrix[i,j]
    }
  }

  return(sim_matrix)

}
