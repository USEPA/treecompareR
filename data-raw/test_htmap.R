my_htmap <- generate_heatmap(
  tree_object = chemont_tree,
  matrix = chemont_jaccard,
  row_data = biosolids_class,
  column_data = usgs_class,
  terminal_only = TRUE,
  row_split = "superclass",
  column_split = "superclass",
  row_title_rot = 0,
  column_title_rot = 90,
  row_title_gp = grid::gpar(fontsize = 6),
  column_title_gp = grid::gpar(fontsize = 6),
  row_gap = unit(0.1, "mm"),
  column_gap = unit(0.1, "mm"))

#To get labels by cluster:
row_inds <- row_order(my_htmap)
names(row_inds)
row_inds[["Homogeneous non-metal compounds"]]

#Convert row indexes into labels:
#get matrix
my_mat <- my_htmap@ht_list[[1]]@matrix
#pull rownames by cluster
row_labs <- lapply(row_inds,
                   function(this_cluster){
                     rownames(my_mat)[this_cluster]
                   })
row_labs[["Homogeneous non-metal compounds"]]

#similarly for columns
column_inds <- column_order(my_htmap)
names(column_inds)
column_inds[["Homogeneous non-metal compounds"]]

#convert column indexes into labels
column_labs <- lapply(column_inds,
                   function(this_cluster){
                     colnames(my_mat)[this_cluster]
                   })
column_labs[["Homogeneous non-metal compounds"]]

#Now, choose a cluster, draw the tree, and highlight branches according to their row/column set membership.

display_subtree(base_tree = chemont_tree,
                prune_to = list(biosolids_class, usgs_class),
                data_1 = row_labs[["Homogeneous non-metal compounds"]],
                name_1 = "Biosolids cluster 1/1",
                data_2 = column_labs[["Homogeneous non-metal compounds"]],
                name_2 = "USGS Water cluster 1/1",
                highlight_by = "set",
                show_tiplabs = FALSE,
                clade_level = 2)


# Try a different heatmap & clustering
my_htmap <- generate_heatmap(
  tree_object = chemont_tree,
  matrix = chemont_jaccard,
  row_data = biosolids_class,
  column_data = usgs_class,
  terminal_only = TRUE,
  row_split = 6L,
  column_split = 6L,
  row_title = "Biosolids",
  column_title = "USGS Water",
  row_title_rot = 0,
  column_title_rot = 90,
  row_title_gp = grid::gpar(fontsize = 6),
  column_title_gp = grid::gpar(fontsize = 6),
  row_gap = unit(0.1, "mm"),
  column_gap = unit(0.1, "mm"))

#To get labels by cluster:
row_inds <- row_order(my_htmap)

#Convert row indexes into labels:
#get matrix
my_mat <- my_htmap@ht_list[[1]]@matrix
#pull rownames by cluster
row_labs <- lapply(row_inds,
                   function(this_cluster){
                     rownames(my_mat)[this_cluster]
                   })

#similarly for columns
column_inds <- column_order(my_htmap)

#convert column indexes into labels
column_labs <- lapply(column_inds,
                      function(this_cluster){
                        colnames(my_mat)[this_cluster]
                      })

#Now, choose a cluster, draw the tree, and highlight branches according to their row/column set membership.
#Examine cluster 4/4
display_subtree(base_tree = chemont_tree,
                prune_to = list(biosolids_class, usgs_class),
                prune_args = list(adjust_branch_length = TRUE),
                data_1 = row_labs[[4]],
                name_1 = "Biosolids cluster 4/4",
                data_2 = column_labs[[4]],
                name_2 = "USGS Water cluster 4/4",
                highlight_by = "set",
                show_tiplabs = FALSE,
                clade_level = 2,
                bg_tree_scale = 0,
                base_opts = list(size = 1))

#Examine cluster 2/5 -- some low-similarity terminal nodes
display_subtree(base_tree = chemont_tree,
                prune_to = list(biosolids_class, usgs_class),
                prune_args = list(adjust_branch_length = TRUE),
                data_1 = row_labs[[2]],
                name_1 = "Biosolids cluster 2/5",
                data_2 = column_labs[[5]],
                name_2 = "USGS Water cluster 2/5",
                highlight_by = "set",
                show_tiplabs = FALSE,
                clade_level = 2,
                bg_tree_scale = 0,
                base_opts = list(size = 1))

#Examine cluster 1/1 -- a high-similarity cluster
display_subtree(base_tree = chemont_tree,
                prune_to = list(biosolids_class, usgs_class),
                prune_args = list(adjust_branch_length = TRUE),
                data_1 = row_labs[[1]],
                name_1 = "Biosolids cluster 1/1",
                data_2 = column_labs[[1]],
                name_2 = "USGS Water cluster 1/1",
                highlight_by = "set",
                show_tiplabs = FALSE,
                clade_level = 2,
                bg_tree_scale = 0,
                base_opts = list(size = 1))

#Examine cluster 2/3 -- both high and low similarity
display_subtree(base_tree = chemont_tree,
                prune_to = list(biosolids_class, usgs_class),
                prune_args = list(adjust_branch_length = TRUE),
                data_1 = row_labs[[2]],
                name_1 = "Biosolids cluster 2/3",
                data_2 = column_labs[[3]],
                name_2 = "USGS Water cluster 2/3",
                highlight_by = "set",
                show_tiplabs = FALSE,
                clade_level = 2,
                bg_tree_scale = 0,
                base_opts = list(size = 1))

cluster_tree <- function(htmap,
                         cluster_row,
                         cluster_column,
                         tree,
                         prune_to = NA,
                         ...){
  htmap <- draw(htmap)
  #To get labels by cluster:
  row_inds <- row_order(htmap)

  #Convert row indexes into labels:
  #get matrix
  my_mat <- htmap@ht_list[[1]]@matrix
  #pull rownames by cluster
  row_labs <- lapply(row_inds,
                     function(this_cluster){
                       rownames(my_mat)[this_cluster]
                     })

  #similarly for columns
  column_inds <- column_order(htmap)

  #convert column indexes into labels
  column_labs <- lapply(column_inds,
                        function(this_cluster){
                          colnames(my_mat)[this_cluster]
                        })

  name_1 <- htmap@row_title
  name_2 <- htmap@column_title

  if(is.na(prune_to)){
    #by default, prune to labels in the heatmap
      prune_to <- do.call(union, dimnames(my_mat))
  }

  display_subtree(base_tree = tree,
                  prune_to = prune_to,
                  prune_args = list(adjust_branch_length = TRUE),
                  data_1 = row_labs[[cluster_row]],
                  name_1 = "Biosolids cluster 4/4",
                  data_2 = column_labs[[cluster_column]],
                  name_2 = "USGS Water cluster 4/4",
                  highlight_by = "set",
                  ...)

}


