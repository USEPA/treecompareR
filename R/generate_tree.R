
# In this function, we generate a random tree with a specified number of tips.
# The internal nodes have degree at least 2 while the root is given degree 2.

#' Generates a random tree with internal nodes of degree at least two.
#' @param n The number of tips, a positive integer.
#' @param rooted Whether the tree is rooted or not.
#' @param max_deg The maximum degree any node can have.
#' @param min_deg The minimum degree any node can have (aside from the root, if rooted)
#' @param seed A seed to allow for replication of results.
#' @return A 'phylo' object representing the generated tree.
#' @export
#' @importFrom ape rtree
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
    phy <- root.phylo(phy, node = ROOT)

  phy
}






# In this function we generate a sequence of integers, each is at least 2, with
# sum equal to the input positive integer. This does not return the input
# integer unless it is 2, 3, and possibly 4.

#' Generates a partition of a positive integer with each element at least 2.
#' @param n positive integer at least 2.
#' @param seed A seed to allow for replication of results.
#' @return A vector of integers each at least 2 and with sum equal to the input.
#' @export
generate_partition <- function(n, seed = NA){
  if(!is.integer(n) | n < 2)
    stop('Please input an integer at least 2!')
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

#' Generates a partition with optional constraints placed on it.
#' @param n Positive integer.
#' @param max_deg Maximal value an element of the partition can take.
#' @param min_deg Minimal value an element of the partition can take.
#' @param seed A seed to allow for replication of results.
#' @return Vector of generated partition.
generate_partition_2 <- function(n, max_deg = NULL, min_deg = 0, seed = NA){
  if(!is.integer(n) | n < 1)
    stop('Please input an integer at least 1!')
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

# This function handles the case where both a max degree and a min degree are
# setting constraints on the partition of n. We assume that n, min_deg, max_deg
# are all positive integers.

#' Generates partition of positive integer, with optional constraints to partition.
#' @param n positive integer.
#' @param max_deg Maximum value an element of the partition can take.
#' @param min_deg Minimum value an element of the partition can take.
#' @param seed A seed to allow for replication of results.
#' @return Vector with partition or error in case specified partition is impossible.
#' @export
generate_partition_3 <- function(n, max_deg, min_deg, seed = NA){
  # Creating a partition x_1 + \cdots + x_k = n, min_deg <= x_j <= max_deg
  # bounds ceiling(n/max_deg) <= k <= floor(n/min_deg). For now, we sample k
  # uniformly and later determine the actual distribution of such partitions.
  lower <- ceiling(n/max_deg)
  upper <- floor(n/min_deg)

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

#' Function to generate several trees
#'
#' @param n_trees The number of trees to generate.
#' @param simulation_seed Optional parameter to allow for replciation of simulations.
#' @param ... Parameters passed to the generate_topology() function.
#' @return A list of 'phylo' objects, each representing a simulated tree.
#' @export
simulate_trees <- function(n_trees = 1, simulation_seed = NA, ...) {
  if(!is.na(simulation_seed)){
    set.seed(simulation_seed)
  }
  runs <- c()
  for (i in 1:n_trees) {
    runs <- c(runs, generate_topology(...))
  }
  return(runs)
}
