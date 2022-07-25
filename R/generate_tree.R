
#' Generate topology
#'
#' This function generates a rooted tree with internal nodes of degree at least
#' two. The sampling of potential topologies and labels assigned to each
#' topology are not completely random. The root is degree two, and all other
#' internal nodes have degrees (potentially) governed by the alternate
#' parameters `max_deg` and `min_deg`.
#'
#' @param n The number of tips, a positive integer.
#' @param rooted Whether the tree is rooted or not.
#' @param max_deg The maximum degree any node can have.
#' @param min_deg The minimum degree any node can have (aside from the root, if
#'   rooted).
#' @param seed A seed to allow for replication of results.
#' @return A 'phylo' object representing the generated tree.
#' @export
#' @importFrom ape rtree
#' @importFrom ape root.phylo
#'
generate_topology <- function(n, rooted = FALSE, max_deg = NULL, min_deg = NULL, seed = NA){
  if (!is.na(seed)){
    set.seed(seed = seed)
  }
  n <- as.integer(n)
  if (!is.null(max_deg)){
    max_deg <- as.integer(max_deg)
    if (is.na(max_deg) | (max_deg < 2))
      stop('The value for max_deg must be an integer at least 2!')
  }
  if (!is.null(min_deg)){
    min_deg <- as.integer(min_deg)
    if (is.na(min_deg) | (min_deg < 2))
      stop('The value for min_deg must be an integer at least 2!')
  }
  if (is.integer(min_deg) & is.integer(max_deg)){
    if (min_deg > max_deg)
      stop('The value for min_deg cannot exceed the value for max_deg!')
  }
  if (n < 1)
    stop("a tree must have at least 1 tip")
  if (n < 3 && !rooted)
    stop("an unrooted tree must have at least 3 tips")
  if (n < 4)
    return(ape::rtree(n, rooted = rooted))
  # Start with two leaves and build on after
  nb <- n - 2L
  if (is.null(max_deg)){
    if (is.null(min_deg)){
      partition <- generate_partition_2(nb, seed = seed)
    } else {
      partition <- generate_partition_2(nb, min_deg = (min_deg-1), seed = seed)
    }
  } else {
    if (is.null(min_deg)){
      partition <- generate_partition_2(nb, max_deg = (max_deg-1), seed = seed)
    } else {
      partition <- generate_partition_3(n = nb, max_deg = (max_deg-1), min_deg = (min_deg - 1), seed = seed)
    }
  }
  # Select the number of surviving choices for splitting and adding new nodes
  splits <- cumsum(c(2, partition)[1:length(partition)])
  x <- as.integer(runif(length(partition)) * splits) + 1L
  breaks <- c(2, cumsum(partition + 2) + 2)
  add_tips <- c(2, cumsum(partition) + 2)

  tip.label <- paste0('t', seq_len(n))
  Nnode <- length(partition) + 1L
  node.label <- paste0('n', seq_len(Nnode))

  TIPS <- sample.int(n)
  N <- 2 + sum(partition + 2)

  edge <- matrix(NA_integer_, N, 2L)
  alive <- logical(N)
  alive[1:2] <- TRUE
  Nalive <- 2L
  e <- 1:2
  ROOT <- n + 1L
  edge[1:2] <- ROOT
  nextnode <- ROOT + 1L
  edge[1:2 + N] <- TIPS[1:2]
  i <- 1L

  #print(partition)
  #print(x)
  #print(breaks)
  #print(add_tips)
  #print(edge)
  #print(alive)
  while (i <= length(partition)) {
    ## draw a branch among the live ones
    k <- which(alive)[x[i]]
    #print(k)
    alive[k] <- FALSE
    e <- (breaks[[i]]+1):breaks[[i+1]]
    alive[e] <- TRUE
    edge[e[1]] <- edge[k]
    edge[e[1] + N] <- nextnode
    edge[e[2:length(e)]] <- nextnode
    edge[e[2] + N] <- edge[k + N]
    edge[e[3:length(e)] + N] <- TIPS[(add_tips[i]+1):add_tips[i+1]]
    nextnode <- nextnode + 1L
    Nalive <- Nalive + partition[i]
    i <- i + 1L
    #print(edge)
    #print(alive)
  }



  edge <- edge[alive, ]
  phy <- list(edge = edge,
              tip.label = tip.label,
              node.label = node.label,
              Nnode = Nnode)
  class(phy) <- "phylo"
  phy <- reorder(phy)

  if (rooted)
    phy <- ape::root.phylo(phy, node = ROOT)

  phy
}






# In this function we generate a sequence of integers, each is at least 2, with
# sum equal to the input positive integer. This does not return the input
# integer unless it is 2, 3, and possibly 4.

#' Partition generator
#'
#' This function generates a partition of an input positive integer. Each
#' constituent of the partition has value at least 2.
#' 2.
#'
#' @param n positive integer at least 2.
#' @param seed A seed to allow for replication of results.
#' @return A vector of integers each at least 2 and with sum equal to the input.
#' @export
#'
#' @examples
#'
#' generate_partition(n = 4, seed = 42)
#' generate_partition(n = 4, seed = 24)
#'
generate_partition <- function(n, seed = NA){
  if(!is.numeric(n) | as.integer(n) < 2)
    stop('Please input an integer at least 2!')
  if (n - as.integer(n) > 0){
    warning(paste('Setting n =', as.integer(n)))
    n <- as.integer(n)
  }


  # if n = 2 or 3, return itself
  if (n < 4) return(n)
  if(!is.na(seed)){
    set.seed(seed)
  }
  partition <- c()
  total <- n
  # Test whether to
  while(total > 4){
    # sample an integer while ensuring what remains is at least 2
    current <- sample(2:(total-2), 1L)
    partition <- c(partition, current)
    total <- total - current
  }
  # If remaining total is 2 or 3, append to the end. If it is 4, randomly assign
  # either 4 or c(2, 2) to the end.
  if (total > 0){
    # randomly assign 4 or c(2,2) to the end if total is 4
    if (total == 4){
      ifelse(runif(1) > .5, partition <- c(partition, 4), partition <- c(partition, 2, 2))
    } else {
      partition <- c(partition, total)
    }
  }
  return(partition)
}

# In this function we generate a sequence of integers, each is at least 1, with
# sum equal to the input positive integer. An optional value, max_deg, can be input to limit the
# maximum value for a constituent member of the partition. An optional value
# min_deg, can be input to limit the minimum value for a constituent member of
# the partition. If min_deg > max_deg or min_deg = max_deg and max_deg does not
# divide n, an error is thrown. If min_deg > n, min_deg is set to 0.

#' Partition generator 2
#'
#' This function generates a partition with optional constraints placed on it.
#' If no such partition can be generated, the function stops and reports this.
#'
#' @param n Positive integer.
#' @param max_deg Maximal value an element of the partition can take.
#' @param min_deg Minimal value an element of the partition can take.
#' @param seed A seed to allow for replication of results.
#' @return Vector of generated partition.
#' @export
#'
#' @examples
#'
#' generate_partition_2(n = 7, min_deg = 3, max_deg = 5, seed = 4)
#' generate_partition_2(n = 7, min_deg = 3, max_deg = 5, seed = 2)
#'
generate_partition_2 <- function(n, max_deg = NULL, min_deg = 0, seed = NA){
  if(!is.numeric(n) | as.integer(n) < 1)
    stop('Please input an integer at least 1!')
  if (n - as.integer(n) > 0){
    warning(paste('Setting n =', as.integer(n)))
    n <- as.integer(n)
  }
  if(!is.na(seed)){
    set.seed(seed)
  }
  if(!is.null(max_deg)){
    max_deg <- as.integer(max_deg)
    if(is.na(max_deg) | max_deg < 1)
      stop('When using max_deg, please input an integer value at least 1 for it!')
  }
  if(min_deg != 0){
    min_deg <- as.integer(min_deg)
    if(is.na(min_deg) | min_deg < 1)
      stop('When using min_deg, please input an integer value at least 1 for it!')
  }
  if(is.integer(min_deg) & is.integer(max_deg)){
    if (min_deg > max_deg){
      stop('When using both min_deg and max_deg, please ensure min_deg <= max_deg!')
    } else if (min_deg == max_deg){
      if (n %% min_deg != 0){
        stop('The min_deg and max_deg are equal, but n is not a multiple!')
      } else {
        return(rep(min_deg, n %/% min_deg))
      }
    } else if (min_deg > 0){
      lower <- min(ceiling(n/max_deg), floor(n/min_deg))
      upper <- max(ceiling(n/max_deg), floor(n/min_deg))

      # For lower <= j <= upper, if there is no such j for which
      # min_deg*j <= n <= max_deg*j, then it is not possible to create a partition
      # with the input parameters.
      if(!any((((lower:upper)*min_deg) <= n) & (((lower:upper)*max_deg) >= n))){
        stop(paste('Such a partition of', n, 'with values between', min_deg,
                   'and', max_deg, 'is not possible to construct!'))
      }
    }
    #else if (max_deg == (2*min_deg)) {
    #quotient <- n %/% min_deg
    #remainder <- n %% min_deg
    #if (remainder > 0) {
    #  if (quotient > 1){
    #    return(c(rep(min_deg, (quotient - 1)), (min_deg + remainder)))
    #  } else if (quotient == 1){
    #    return(n)
    #    } else {
    #  min_deg <- 0
    #  warning('Setting min_deg = 0 since min_deg > n!')
    #}
    #} else {
    #  return(rep(min_deg, quotient))
    #}


  }
  if(is.integer(min_deg) & (min_deg > n)){
    min_deg <- 0
    warning('Setting min_deg = 0 since min_deg > n!')
  }

  # if n = 1, return itself
  if (n < 2) return(n)
  partition <- c()
  total <- n
  while(total > max(1, 2*min_deg)){
    # sample an integer
    current <- sample(max(1, min_deg):(min(total - min_deg, max_deg)), 1L)
    partition <- c(partition, current)
    total <- total - current
  }
  if (total == 1){
    partition <- c(partition, 1)
  } else if (total == (2*min_deg) & (min_deg !=0)){
    ifelse(runif(1) > .5, partition <- c(partition, as.integer(2*min_deg)), partition <- c(partition, min_deg, min_deg))
  } else if (total > 0){
    partition <- c(partition, total)
  }
  return(partition)
}

#' Partition generator 3
#'
#' This function generates a partition of an input positive integer subject to
#' both minimum degree and maximum degree constraints.
#'
#' @param n positive integer.
#' @param max_deg Maximum value an element of the partition can take.
#' @param min_deg Minimum value an element of the partition can take.
#' @param seed A seed to allow for replication of results.
#' @return Vector with partition or error in case specified partition is
#'   impossible.
#' @export
#'
generate_partition_3 <- function(n, max_deg, min_deg, seed = NA){
  # Creating a partition x_1 + \cdots + x_k = n, min_deg <= x_j <= max_deg
  # bounds ceiling(n/max_deg) <= k <= floor(n/min_deg). For now, we sample k
  # uniformly and later determine the actual distribution of such partitions.
  lower <- min(floor(n/min_deg), ceiling(n/max_deg))
  upper <- max(floor(n/min_deg), ceiling(n/max_deg))

  # For lower <= j <= upper, if there is no such j for which
  # min_deg*j <= n <= max_deg*j, then it is not possible to create a partition
  # with the input parameters.
  if(!any((((lower:upper)*min_deg) <= n) & (((lower:upper)*max_deg) >= n))){
    stop(paste('Such a partition of', n, 'with values between', min_deg,
               'and', max_deg, 'is not possible to construct!'))
  }
  if(!is.na(seed)){
    set.seed(seed)
  }
  if (lower == upper){
    k <- lower
  } else {
    k <- sample(lower:upper, 1)
  }

  # If min_deg*upper >

  #print(lower)
  #print(upper)
  #print(k)

  # We transform the problem as follows: y_j = x_j - min_deg, and
  # y_1 + \cdots + y_k = n - k(min_deg). We place max_deg - min_deg copies of
  # each index 1, 2, \cdots , k into a bag and choose n - k(min_deg) from these.
  # we then determine how many copies of each index j are chose and assign that
  # value to y_j. Then we add min_deg to get the value of x_j.
  bag <- rep(1:k, (max_deg - min_deg))
  #print(length(bag))
  values <- sample(x = bag, size = (n - k*min_deg))
  #print(values)
  y_values <- integer(k)
  for (i in 1:k){
    y_values[i] <- length(which(values == i))
  }
  x_values <- y_values + min_deg

  return(x_values)
}

#' Tree simulator
#'
#' This function generates several trees subject to input constraints.
#'
#' @param n_trees The number of trees to generate.
#' @param simulation_seed Optional parameter to allow for replication of
#'   simulations.
#' @param n The number of tips each tree will have.
#' @param ... Parameters passed to the generate_topology() function.
#' @return A list of 'phylo' objects, each representing a simulated tree.
#' @export
#'
simulate_trees <- function(n_trees = 1, simulation_seed = NA, n, ...) {
  if(!is.na(simulation_seed)){
    set.seed(simulation_seed)
  }
  runs <- list(n_trees)
  for (i in 1:n_trees) {
    new_tree <- tryCatch(
      {
        generate_topology(n = n, ...)
        },
      error = function(e) {
        message(paste(e, '\n'))
        return(NA)
      }#,
      #warning = function(w){
      #  message(w)
      #}
    )
    #print(str(new_tree))
    #print(length(new_tree))
    if (length(new_tree) > 1){
      runs[[i]] <- new_tree
    } else {
      stop('There was an error... see above!')
    }

  }
  return(runs)
}

#' Generate caterpillar tree
#'
#' This function generates a caterpillar tree, a rooted binary tree of maximal
#' depth with n tips.
#'
#' @param n The number of tips.
#' @return A `phylo` object representing the generated tree.
#' @export
generate_caterpillar <- function(n){
  n <- as.integer(n)
  if (n < 2){
    stop('Please input an integer at least 2!')
  }

  edge <- matrix(NA_integer_, nrow = (2*n - 2), ncol = 2)

  for (i in 1:(n-1)){
    edge[i, 1] <- n + i
    edge[i, 2] <- i
  }

  for (i in 1:(n-2)){
    edge[i + n, 1] <- i + n
    edge[i + n, 2] <- i + n + 1
  }

  edge[n, 1] <- 2*n-1
  edge[n, 2] <- n


  tip.label <- paste0("t", 1:n)
  node.label <- paste0('n', 1:(n-1))
  Nnode <- n - 1

  phy <- list(edge = edge,
              tip.label = tip.label,
              node.label = node.label,
              Nnode = Nnode)
  class(phy) <- "phylo"
  phy <- reorder(phy)

  phy <- ape::root.phylo(phy, node = n+1)

  return(phy)
}

#' Generate Star tree
#'
#' This function generates a star tree, a rooted tree of minimal
#' depth with n tips.
#'
#' @param n The number of tips.
#' @return A `phylo` object representing the generated tree.
#' @export
generate_star <- function(n){
  n <- as.integer(n)

  if (n < 2){
    stop('Please input an integer at least 2!')
  }

  edge <- matrix(NA_integer_, nrow = n, ncol = 2)
  edge[, 1] <- n+1
  edge[, 2] <- 1:n

  tip.label <- paste0('t', 1:n)
  node.label <- c('n1')
  Nnode <- 1

  phy <- list(edge = edge,
              tip.label = tip.label,
              node.label = node.label,
              Nnode = Nnode)
  class(phy) <- "phylo"
  phy <- reorder(phy)

  phy <- ape::root.phylo(phy, node = n+1)

  return(phy)

}

#'Generate balanced tree
#'
#'This function will generate a rooted binary tree that is balanced (or as close
#'to balanced) with n tips.
#'
#'@param n The number of tips
#'@return A `phylo` object representing the generated tree.
generate_balanced <- function(n){
  n <- floor(as.integer(n))

  if (n < 2){
    stop('Please input an integer at least 2!')
  }


  Nnode <- n - 1

  tip.label <- paste0('t', 1:n)
  node.label <- paste0('n', 1:Nnode)



  edge <- matrix(NA_integer_, nrow = 2*n-2, ncol = 2)
  edge[, 1] <- rep((n+1):(2*n-1), each = 2)



  if (identical(n, 2^floor(log2(n)))){
    edge[, 2] <- c((n+2):(2*n-1), 1:n)
  } else {
    next_layer <- 1:(2*(n - 2^floor(log2(n))))
    previous_layer <- (2*(n-2^floor(log2(n))) + 1):n
    edge[, 2] <- c((n+2):(2*n - 1), previous_layer, next_layer)
  }

  phy <- list(edge = edge,
              tip.label = tip.label,
              node.label = node.label,
              Nnode = Nnode)

  class(phy) <- "phylo"
  phy <- reorder(phy)

  phy <- ape::root.phylo(phy, node = n+1)

  #for (i in 2:(n-1)){
  #  phy <- ape::rotate(phy = phy, node = (n+i))
  #}

  #phy <- ape::rotate(phy = phy, node = (n+1))

  return(phy)

}
