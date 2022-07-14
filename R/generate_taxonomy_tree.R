
# In this function, we read in the JSON file that contains the parent child
# relationships in a taxonomy. The JSON file should have three pieces of data
# for each entry: 'Name', 'ID', 'Parent_ID'

#' Parent-child relationship
#'
#' This is a helper function that reads in a JSON file with the parent-child
#' relationships describing a taxonomy and returns a data.frame with this data.
#'
#' @param file_name Path to JSON file
#' @return A data.frame with parent-object relationships.
#' @importFrom jsonlite fromJSON
#'
get_parent_child <- function(file_name){
  temp <- jsonlite::fromJSON(txt = file_name)
  names(temp) <- c('Name', 'ID', 'Parent_ID')
  temp
}

#' Path to root
#'
#' This helper function generates a path from the root to a given node from a
#' data.frame of parent-child relationships. This is used in creating a tree
#' based off of the data.frame of parent-child relationships.
#'
#' @param df data.frame with parent-child relationships
#' @param x Numeric identification number for node
#' @return Path from node to root, as a string delimited by '/'.
#'
get_path_to_root <- function(df, x){
  if (x < 0){
    return('root_')
  }

  current_index <- which(df[, 'id_number'] %in% x)
  #print(current_index)
  current_parent <- df[current_index, 'parent_id_number']
  #print(current_parent)
  path_to_root <- paste(get_path_to_root(df, current_parent), x, sep = '/')

  return(path_to_root)
}

# Slow conversion workaround suggested on github

#' Plot height
#'
#' This helper function allows for fast tree construction. See
#' \href{https://github.com/gluc/data.tree/issues/92#issuecomment-299571142}{Slow
#' conversion between data.tree and Newick format} for background on this helper
#' function.
#'
#' @param node A data.tree object.
#' @param rootHeight A positive number placeholder for root height, default 100.
#' @return A data.tree object with plotHeight attribute.
#' @import data.tree
#'
SetPlotHeight <- function(node, rootHeight = 100) {

  #traverse from leaves towards root to calculate the height and store it in height2
  #Note: cannot call it height, as that attribute already exists on `Node`
  data.tree::Set(node$leaves, height2 = 1)
  node$Do(function(x) x$height2 <- data.tree::Aggregate(x, "height2", max) + 1, traversal = "post-order", filterFun = data.tree::isNotLeaf)

  #traverse from root to leaves to calculate the plotHeight
  #(where plotHeight uses the semantics of dendrogram/phylo, where leaves should have height 0 and root height = rootHeight. This meaning is not the same as in the data.tree package)
  node$plotHeight <- rootHeight
  node$Do(function(x) x$plotHeight <- x$parent$plotHeight * (1 - 1 / x$height2), filterFun = data.tree::isNotRoot)
}


#' Generate taxonomy tree
#'
#' This function will convert taxonomy relationships from a data.frame to a
#' 'phylo' object. In addition to the 'phylo' object, there is a data.frame
#' returned that includes data relevant to the creation of the 'phylo' object.
#'
#' This function takes a slightly circuitous route. It converts the
#' \code{data.frame} input argument \code{tax_nodes} to an object of class
#' \code{\link[data.tree]{Node}}. Then, it converts that
#' \code{\link[data.tree]{Node}} object into a Newick-format tree
#' representation. Finally, it converts that Newick-format tree representation
#' into an object of class \code{\link[phytools]{phylo}}.
#'
#' This is necessary because \code{\link[ggtree]{ggtree}} requires
#' \code{\link[phytools]{phylo}}-class objects as input, but
#' \code{\link[phytools]{phytools}} only has methods to create
#' \code{\link[phytools]{phylo}}-class objects from Newick, SIMMAP, or
#' Nexus-formatted trees. However, the ClassyFire ontology is not available in
#' any of those three formats, but only as files defining the name and ID number
#' of each node in the ontology, and giving the ID number of each node's parent
#' (if any). For example, the ClassyFire ontology can be downloaded in JSON
#' format at \url{http://classyfire.wishartlab.com/tax_nodes.json} or in OBO
#' format at \url{http://classyfire.wishartlab.com/downloads}. (Both JSON and
#' OBO files contain the same information.) The most-accessible method for
#' creating a Newick-format tree from this information is found in
#' \code{\link[data.tree]{ToNewick}}, which requires a
#' \code{\link[data.tree]{Node}} or \code{\link[ape]{phylo}}-class input.
#'
#' The argument \code{tax_nodes} is a \code{data.frame} that defines the
#' taxonomy. For each node in the taxonomy, it contains a name, an ID, and the
#' ID of the node's parent. This \code{data.frame} may be read from an OBO file
#' using e.g. \code{\link[ontologyIndex]{get_ontology}} (and
#' \code{\link[ontologyIndex]{as.data.frame.ontology_index}}), or from a JSON
#' file using for example \code{\link[jsonlite]{fromJSON}}.
#'
#' @param tax_nodes A data.frame containing parent-child relationships for the
#'   taxonomy. Must contain columns "Name", "ID", and "Parent_ID", which
#'   respectively provide a name for each node, an ID for each node, and the ID
#'   of the parent of each node. Each row represents one node of the taxonomy.
#' @return A list containing two elements. \code{tree_object} is a
#'   \code{phylo}-class tree object.\code{tax_nodes} is a data.frame with as
#'   many rows as there are nodes in the taxonomy, and nine variables.
#'   \describe{ \item{Name}{Text label for the node, as in the input data.frame}
#'   \item{ID}{Original ID for the node, as in the input data.frame}
#'   \item{Parent_ID}{Parent ID for the node, as in the input data.frame}
#'   \item{id_number}{Internally generated node ID number}
#'   \item{parent_id_number}{Internally generated parent ID number}
#'   \item{pathString}{Path from root to the node as a slash-separated string of
#'   internally-generated ID numbers} \item{plotHeight}{Internally-generated
#'   plot height for each node; currently does not represent any real data, only
#'   used internally to speed plotting} \item{level}{Taxonomic level of the
#'   node, where level 0 is the root, level 1 are the root's immediate children,
#'   etc.} \item{phylo_node}{Node number for this label in the phylo tree
#'   object, which is different from the internally-generated node number. This
#'   is useful for plotting.} }
#' @export
#'
generate_taxonomy_tree <- function(tax_nodes = NULL){
  # Get the root of the tree, which will have no parent and thus have null value
  # for the 'Parent_ID'
  root <- tax_nodes[which(!(tax_nodes[, 'Parent_ID'] %in% tax_nodes[, 'ID'])), 'ID']

  # Generate columns that enumerate the 'ID' and 'Parent_ID' values. Change the
  # id_number value of the root to -1.
  tax_nodes$id_number <- 1:(dim(tax_nodes)[[1]])
  tax_nodes[tax_nodes$ID == root, 'id_number'] <- -1
  tax_nodes$parent_id_number <- tax_nodes$id_number[match(tax_nodes$Parent_ID, tax_nodes$ID)]

  # Generate the pathString column that is used to build the tree
  tax_nodes$pathString <- sapply(tax_nodes$id_number, function(s) {
    get_path_to_root(tax_nodes, s)
  })

  # Convert the pathString information into an object of 'Node' class
  tax_tree <- data.tree::as.Node(tax_nodes, mode = 'table')


  # Set plot heights for nodes and convert into Newick string
  SetPlotHeight(tax_tree)
  plot_heights <- data.frame(name = tax_tree$Get('name'),
                             plot_height = tax_tree$Get('plotHeight'))
  plot_heights[1, 'name'] <- 	'-1'
  tax_nodes['plotHeight'] <- plot_heights[match(tax_nodes[,'id_number'], plot_heights[,'name']) ,'plot_height']

  Newick_tax_tree <- data.tree::ToNewick(tax_tree,
                                         heightAttribute = 'plotHeight')
  tree_object <- phytools::read.newick(text = Newick_tax_tree)

  # Rename tip and node labels to original names from the JSON file
  #tip labels:
  temp_index <- match(tree_object$tip.label,
                      tax_nodes$id_number)
  temp_labels <- tax_nodes$Name[temp_index]
  tree_object$tip.label <- temp_labels

  #node labels:
  temp_index <- match(tree_object$node.label,
                      tax_nodes$id_number)
  temp_labels <- tax_nodes$Name[temp_index]
  #handle root node label
  temp_labels[tree_object$node.label %in% "root_"] <- tax_nodes[tax_nodes$id_number == -1, 'Name']
  tree_object$node.label <- temp_labels

  #Get tree object node numbers (different from internally-generated)
  #at each level of the taxonomy
  foo <- phangorn::Ancestors(tree_object,
                             node = 1:Ntip(tree_object))

  #this is a list of ancestors for each tip node
  #read them backwards
  foo_rev <- lapply(foo, rev)
  #add the tip node number too
  foo_rev <- sapply(1:Ntip(tree_object),
                function(i) c(foo_rev[[i]], i))
  max_level <- max(sapply(foo, length))
  level_nodes <- vector(mode = "list", length = max_level)
  level_names <- vector(mode = "list", length = max_level)
  for (i in 1:max_level){
    #get the node numbers at this level
    level_nodes[[i]] <- unlist(
      unique(
        sapply(foo_rev,
               function(x){
                 if(length(x)>=i){
                   x[i]
                 }else{
                   NULL
                 }
               }
        )
      )
    )
#get the corresponding node labels at this level
    level_names[[i]] <- select.tip.or.node(element= level_nodes[[i]],
                                           tree = tree_object)
    #assign level number to each label in tax_nodes
    tax_nodes[tax_nodes$Name %in% level_names[[i]], "level"] <- i
  }

  nodenames <- data.frame(phylo_node = unlist(level_nodes),
                          Name = unlist(level_names))

  tax_nodes <- merge(tax_nodes, nodenames, by = "Name")
  tax_nodes$level <- tax_nodes$level - 1 #so that root is level 0

  return(list("tree_object" = tree_object,
              "tax_nodes" = tax_nodes))
}

#' Get label for a tip or internal node of a phylo tree
#'
#' @param element Vector of node numbers in phylo tree
#' @param tree phylo tree object
#' @return Character vector of tip or internal node labels
select.tip.or.node <- function(element, tree) {
  label <- vector(mode = "character", length = length(element))
  label[element < (Ntip(tree)+1)] <- tree$tip.label[element[element <( Ntip(tree)+1)]]
  label[element >  Ntip(tree)] <- tree$node.label[element[element >  Ntip(tree)]-Ntip(tree)]
  return(label)
}
