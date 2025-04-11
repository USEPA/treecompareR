my_tree <- prune_tree(chemont_tree,
                      prune_to = dplyr::bind_rows(biosolids_class,
                                                  usgs_class))
my_tree <- bind_entities(my_tree,
                         data = dplyr::bind_rows(
                           biosolids_class, usgs_class),
                         entity_id_col = "CASRN")
my_tree_mat <- similarity_matrix(tree = my_tree,
                                 upper_tri = FALSE)
generate_heatmap(tree_object = my_tree,
                 matrix = my_tree_mat,
                 row_data = biosolids_class,
                 column_data = usgs_class,
                 entity_id_col = "CASRN",
                 row_split = 5L,
                 column_splot = 5L)
