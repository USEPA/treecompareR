my_tree <- prune_tree(chemont_tree,
                      prune_to = dplyr::bind_rows(biosolids_class,
                                                  usgs_class))
my_tree2 <- bind_entities(my_tree,
                         data = dplyr::bind_rows(
                           biosolids_class, usgs_class),
                         entity_id_col = "PREFERRED.NAME")

#add entities as new layer of classification
my_biosolids_class <- biosolids_class %>%
  #add a potential 12th level
  dplyr::mutate(level12 = NA_character_) %>%
  #pivot longer
  tidyr::pivot_longer(cols = tidyselect::all_of(
    c(chemont_tax_levels,
      "level12")),
    names_to = "level",
    values_to = "label") %>%
  dplyr::mutate(level_number = match(level, c(chemont_tax_levels, "level12"))) %>%
  #find the first NA label for each entity and substitute with CASRN
  dplyr::mutate(label = dplyr::if_else(level_number == (terminal_level + 1),
                                       PREFERRED.NAME,
                                       label)) %>%
  dplyr::select(-level_number) %>%
  tidyr::pivot_wider(names_from = level,
                     values_from = label) %>%
  add_terminal_label(entity_id_col = "DTXSID", tax_level_labels = c(chemont_tax_levels,
                                                                    "level12"))

my_usgs_class <- usgs_class %>%
  #add a potential 12th level
  dplyr::mutate(level12 = NA_character_) %>%
  #pivot longer
  tidyr::pivot_longer(cols = tidyselect::all_of(
    c(chemont_tax_levels,
      "level12")),
    names_to = "level",
    values_to = "label") %>%
  dplyr::mutate(level_number = match(level, c(chemont_tax_levels, "level12"))) %>%
  #find the first NA label for each entity and substitute with CASRN
  dplyr::mutate(label = dplyr::if_else(level_number == (terminal_level + 1),
                                       PREFERRED.NAME,
                                       label)) %>%
  dplyr::select(-level_number) %>%
  tidyr::pivot_wider(names_from = level,
                     values_from = label) %>%
  add_terminal_label(entity_id_col = "DTXSID", tax_level_labels = c(chemont_tax_levels,
                                                                    "level12"))

system.time(
  my_tree_mat2 <- similarity_matrix(tree = my_tree2,
                                 nodes1 = NULL,
                                 nodes2 = NULL,
                                 upper_tri = FALSE)
)

#plot tree highlighted by set membership
treeplot_set <- display_subtree(base_tree = my_tree2,
                prune_to = "Organohalogen compounds",
tax_level_labels = c(chemont_tax_levels,
"level12"),
entity_id_col = "PREFERRED.NAME",
data_1 = my_biosolids_class,
data_2 = my_usgs_class,
sim_mat = my_tree_mat2,
highlight_by = "set",
base_opts = list(size = 1))

#and the same tree highlighted by similarity
treeplot_sim <- display_subtree(base_tree = my_tree2,
                prune_to = "Organohalogen compounds",
tax_level_labels = c(chemont_tax_levels,
"level12"),
entity_id_col = "PREFERRED.NAME",
data_1 = my_biosolids_class,
data_2 = my_usgs_class,
sim_mat = my_tree_mat2,
highlight_by = "sim",
base_opts = list(size = 1))

#now... create heatmap of similarity of entities only.
#

my_htmap <- generate_heatmap(row_data = my_biosolids_class,
                 column_data = my_usgs_class,
                 terminal_only = TRUE,
                 tree_object = my_tree2,
                 matrix = my_tree_mat2,
                 entity_id_col = "CASRN",
                 row_split = 9L,
                 column_split = 9L)


#tree highlighted by membership in cluster 3/2
row_names <- dimnames(my_htmap@ht_list[[1]]@matrix)[[1]][
  stats::order.dendrogram(
    ComplexHeatmap::row_dend(my_htmap)[[3L]]
  )
]
col_names <- dimnames(my_htmap@ht_list[[1]]@matrix)[[2]][
  stats::order.dendrogram(
    ComplexHeatmap::column_dend(my_htmap)[[2L]]
  )
]

display_subtree(base_tree = my_tree2,
                tax_level_labels = c(chemont_tax_levels,
                                     "level12"),
                entity_id_col = "PREFERRED.NAME",
                data_1 = row_names,
                data_2 = col_names,
                name_1 = "Cluster 3/2 biosolids",
                name_2 = "Cluster 3/2 USGS",
                highlight_by = "set")
