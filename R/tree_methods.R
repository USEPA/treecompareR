# This generates the number of descendants for each node. The edges are stored
# in a two column data.frame, with the first column indicating the parent and
# the second column indicating the child. Internal numbering of nodes is given
# by tips first and then nodes, with the root given by the first node of the
# internal nodes. The output also includes the number of children for each node
# and the node level.

#' Generate descendants
#'
#' This function generates a data.frame of the number of descendants for each
#' node of the input tree.
#'
#' @param tree A phylo object
#' @return data.frame consisting of the node number, descendants, children, and
#'   level
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

#' Generate information content.
#'
#' This function generates a data.frame of information content for the input
#' tree. This uses the formulation as described in
#' \href{https://www.researchgate.net/publication/220837848_An_Intrinsic_Information_Content_Metric_for_Semantic_Similarity_in_WordNet/stats}{An
#' Intrinsic Information Content Metric for Semantic Similarity in WordNet}. The
#' data.frame also includes the depth of each tip and internal node.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param log_descendants Alternate parameter for specifying type of information
#'   content.
#' @return data.frame consisting of node number, children, descendants, level,
#'   and information content for each node
#' @export
#' @importFrom ape is.rooted
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' generate_information_content(tree = tree)
#'
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

#' Ancestors
#'
#' Generates a list of ancestors for a given node in a tree.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label The node label.
#' @param node_number Alternate parameter, the number of the given node.
#' @return A list of nodes back to the root of ancestors for the given node.
#' @export
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' get_ancestors(tree = tree, label = 't2')
#' get_ancestors(tree = tree, label = 'n1')
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


#' Tree level
#'
#' This function returns the tree level of the given node in a rooted tree.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label The node label.
#' @param node_number Alternate parameter, the number of the given node.
#' @return The level of the node from the root of the tree.
#' @export
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' get_tip_level(tree = tree, label = 't2')
#' get_tip_level(tree = tree, label = 'n1')
get_tip_level <- function(tree, label, node_number = NULL){
  return(length(get_ancestors(tree = tree, label = label, node_number = node_number)))
}


#' Jaccard distance
#'
#' This determines the Jaccard distance for two input node labels in a given
#' tree. For each node, there is a set of labels along the unique path from the
#' root to the label. These sets are compared using Jaccard distance. For more
#' information on Jaccard distance, please consult
#' \href{https://en.wikipedia.org/wiki/Jaccard_index}{Jaccard Index}.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label_A The first label.
#' @param label_B The second label.
#' @return The Jaccard distance of the label sets for the root to node path.
#' @export
#'
#' @seealso \code{\link{general_Jaccard_similarity}}
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' general_Jaccard_dist(tree = tree, label_A = 't2', label_B = 't4')
#' general_Jaccard_dist(tree = tree, label_A = 'n3', label_B = 't8')
#'
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

#' Jaccard similarity
#'
#' This determines the Jaccard similarity for two input node labels in a given
#' tree. For each node, there is a set of labels along the unique path from the
#' root to the label. These sets are compared using Jaccard similarity. For more
#' information on Jaccard similarity, please consult
#' \href{https://en.wikipedia.org/wiki/Jaccard_index}{Jaccard Index}.
#'
#' @param tree A phylo object representing a rooted tree.
#' @param label_A The first label.
#' @param label_B The second label.
#' @return The Jaccard distance of the label sets for the root to node path.
#' @export
#'
#' @seealso \code{\link{general_Jaccard_dist}}
#'
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' general_Jaccard_similarity(tree = tree, label_A = 't2', label_B = 't4')
#' general_Jaccard_similarity(tree = tree, label_A = 'n3', label_B = 't8')
#'
general_Jaccard_similarity <- function(tree, label_A, label_B){
  return(1 - general_Jaccard_dist(tree, label_A, label_B))
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
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' general_Resnik_similarity(tree = tree, label_A = 't2', label_B = 't4')
#' general_Resnik_similarity(tree = tree, label_A = 'n3', label_B = 't8')
#'
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
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' general_Lin_similarity(tree = tree, label_A = 't2', label_B = 't4')
#' general_Lin_similarity(tree = tree, label_A = 'n3', label_B = 't8')
#'
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
#' @examples
#'
#' tree <- generate_topology(n = 8, rooted = TRUE, seed = 42)
#'
#' general_JiangConrath_similarity(tree = tree, label_A = 't2', label_B = 't4')
#' general_JiangConrath_similarity(tree = tree, label_A = 'n3', label_B = 't8')
#'
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
#' MonteCarlo_similarity(tree = treecompareR:::chemont_tree, data_1 = dt1, data_2 = dt2, name_1 = 'Biosolids',
#'                      name_2 = 'USGS', label_number = 200, seed = 42L, Jaccard = chemont_jaccard,
#'                      Resnik = chemont_resnik_IC_SVH, Lin = chemont_lin_IC_SVH,
#'                      JiangConrath = chemont_jiangconrath_IC_SVH)}
#'
MonteCarlo_similarity <- function(tree, data_1, data_2, data_1_indices = NULL, data_2_indices = NULL,
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


  Nnode <- length(tree$node.label)
  dimnames <- c(tree$tip.label, tree$node.label[2:Nnode])

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

    #print(all_node_indices)
    #print(all_tip_indices)

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
      new_row[[7]] <- mean(Lin[all_tip_indices, dataset_1_indices][upper.tri(Lin[all_tip_indices, dataset_1_indices], diag = TRUE)])
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

  names(simulation_dataframe)[5:8] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), name_1)
  names(simulation_dataframe)[9:12] <- paste0(c('Jaccard', 'Resnik', 'Lin', 'JiangConrath'), name_2)

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


  margin <- min(unique_avgs[2:length(unique_avgs)] - unique_avgs[1:(length(unique_avgs)-1)])/3

  percentages <- sapply(rev(unique_avgs), function(t) {sum(temp_counts[sort(average_val, decreasing = TRUE) > (t - margin)])/total})

  names(percentages) <- rev(unique_avgs)

  if (!is.na(cutoff)) {
    # returns maximum value that achieves cutoff percentage threshold
    return(names(percentages)[[min(which(percentages >= cutoff))]])
  } else {
    return(percentages)
  }
  }
