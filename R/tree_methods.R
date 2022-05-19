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
#' @importFrom ape is.rooted
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
#' @importFrom ape is.rooted
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
#' @importFrom ape is.rooted
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @return A similarity matrix, which is a symmetric matrix with values in \eqn{[0,1]}.
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


#' This function simulates trees for similarity.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param data_1 A data.table of chemicals with classifications.
#' @param data_2 A data.table of chemicals with classifications.
#' @param data_1_indices Alternate parameter giving indices of nodes of a subtree of `tree`.
#' @param data_2_indices Alternate parameter giving indices of nodes of a subtree of `tree`.
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
#' @return A data.frame with similarity values for each simulation.
#' @export
#' @importFrom data.table is.data.table
#' @importFrom phangorn Ancestors
MonteCarlo_similarity <- function(tree, data_1, data_2, data_1_indices = NULL, data_2_indices = NULL, name_1 = 'data_set_1', name_2 =  'data_set_2', label_number = 100, repetition = 10, seed = NA_real_, only_tips = FALSE, Jaccard = NULL, Resnik = NULL, Lin = NULL, JiangConrath = NULL){
  if (!is.na(seed) & is.integer(seed)){
    set.seed(seed)
  }
  if(only_tips & (label_number > length(tree$tip.label))){
    stop('Please input a label_number less than the number of tips!')
  } else if (label_number > (length(tree$tip.label) + length(tree$node.label))) {
    stop('Please input a label_number less than the number of nodes and tips!')
  }
  if (!(data.table::is.data.table(data_1) & data.table::is.data.table(data_2))){
    if (is.null(data_1_indices) | is.null(data_2_indices)){
      stop('Please input either indices for `data_1_indices` and `data_2_indices` or data.table objects for `data_1` and `data_2`')
    } else {
      stop('Please input a data.table object for each of the `data_1` and `data_2` parameters!')
    }
      }
  if (is.null(Jaccard) & is.null(Resnik) & is.null(Lin) & is.null(JiangConrath)){
    stop('Please input a similarity matrix for at least one of Jaccard, Resnik, Lin, and JiangConrath parameters!')
  }

  simulation_dataframe <- data.frame(Jaccard_all = double(),
                                     Resnik_all = double(),
                                     Lin_all = double(),
                                     JiangConrath_all = double(),
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


  dimnames <- c(tree$tip.label, tree$node.label)

  if (is.data.table(data_1) & is.data.table(data_2)) {
    dataset_1_labels <- unlist(get_labels(data = data_1))
    dataset_1_indices <- which(dimnames %in% dataset_1_labels)

    dataset_2_labels <- unlist(get_labels(data = data_2))
    dataset_2_indices <- which(dimnames %in% dataset_2_labels)
  } else {
    dataset_1_labels <- dimnames[data_1_indices]
    dataset_1_indices <- data_1_indices

    dataset_2_labels <- dimnames[data_2_indices]
    dataset_2_labels <- data_2_indices
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

    new_row <- double(14L)

    if (!is.null(Jaccard)){
      new_row[[1]] <- mean(Jaccard[all_node_indices, all_node_indices][upper.tri(Jaccard[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[5]] <- mean(Jaccard[all_tip_indices, dataset_1_indices][upper.tri(Jaccard[all_tip_indices, dataset_1_indices], diag = TRUE)])
      new_row[[9]] <- mean(Jaccard[all_tip_indices, dataset_2_indices][upper.tri(Jaccard[all_tip_indices, dataset_2_indices], diag = TRUE)])
    }

    if (!is.null(Resnik)){
      new_row[[2]] <- mean(Resnik[all_node_indices, all_node_indices][upper.tri(Resnik[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[6]] <- mean(Resnik[all_tip_indices, dataset_1_indices][upper.tri(Resnik[all_tip_indices, dataset_1_indices], diag = TRUE)])
      new_row[[10]] <- mean(Resnik[all_tip_indices, dataset_2_indices][upper.tri(Resnik[all_tip_indices, dataset_2_indices], diag = TRUE)])
    }

    if (!is.null(Lin)){
      new_row[[3]] <- mean(Lin[all_node_indices, all_node_indices][upper.tri(Lin[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[7]] <- mean(Lin[all_tip_indices, dataset_1_indices][upper.tri(lin[all_tip_indices, dataset_1_indices], diag = TRUE)])
      new_row[[11]] <- mean(Lin[all_tip_indices, dataset_2_indices][upper.tri(Lin[all_tip_indices, dataset_2_indices], diag = TRUE)])
    }

    if (!is.null(JiangConrath)){
      new_row[[4]] <- mean(JiangConrath[all_node_indices, all_node_indices][upper.tri(JiangConrath[all_node_indices, all_node_indices], diag = TRUE)])
      new_row[[8]] <- mean(JiangConrath[all_tip_indices, dataset_1_indices][upper.tri(JiangConrath[all_tip_indices, dataset_1_indices], diag = TRUE)])
      new_row[[12]] <- mean(JiangConrath[all_tip_indices, dataset_2_indices][upper.tri(JiangConrath[all_tip_indices, dataset_2_indices], diag = TRUE)])
    }

    new_row[[13]] <- length(all_node_indices)
    new_row[[14]] <- length(all_tip_indices)

    simulation_dataframe[i, ] <- new_row
  }

  names(simulation_datafram)[5:8] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), name_1)
  names(simulation_datafram)[9:12] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), name_2)

  return(simulation_dataframe)
}

#' This function gives cutoffs for similarity level and subtree representation
#'
#' @param tree A phylo object representing a rooted tree.
#' @param matrix A similarity matrix corresponding to `tree`.
#' @param data A data.table of chemicals with classifications.
#' @param neighbors A parameter giving how many neighbors to use for finding label average values.
#' @param cutoff An alternate parameter giving the cutoff percentage value.
#' @param labels An alternate parameter giving a list of node labels corresponding to a subtree of `tree`.
#' @param counts An alternate parameter giving the counts of occurrence for each label.
#' @return Named list of percentage of data represented by similarity values.
#' @export
get_cutoffs <- function(tree, matrix, data, tax_levels = NULL, neighbors = 3, cutoff = NA_real_, labels = NULL, counts = NULL){
  if (is.data.table(data)){
    if (is.null(tax_levels)){
      tax_levels <- c('kingdom', 'superclass', 'class', 'subclass',
                      'level5', 'level6', 'level7', 'level8',
                      'level9', 'level10', 'level11')
    }
    counts <- get_number_of_labels(data = data, tax_levels = tax_levels)
    labels <- names(counts)
    #labels <- get_labels(data = data, tax_levels = tax_levels)
  }


  indices <- which(dimnames(matrix)[[1]] %in% labels)

  temp_mat <- mat[indices, indices]

  average_val <- unname(apply(temp_mat, MARGIN = 1, function(t) {sum(sort(t, decreasing = TRUE)[1:3])/3}))

  total = sum(counts)

  temp_counts <- counts[order(average_val, decreasing = TRUE)]

  unique_avgs <- sort(unique(average_val))

  margin <- min(unique_avgs[2:length(unique_avgs)] - unique_avgs[1:(length(unique_avgs)-1)])/3

  percentages <- sapply(rev(unique_avgs), function(t) {sum(temp_counts[sort(average_three, decreasing = TRUE) > (t - margin)])/total})

  names(percentages) <- rev(unique_avgs)

  if (!is.na(cutoff)) {
    # returns maximum value that achieves cutoff percentage threshold
    return(names(percentages)[[min(which(percentages >= cutoff))]])
  } else {
    return(percentages)
  }
  }
