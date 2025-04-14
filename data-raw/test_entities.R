my_tree <- prune_tree(chemont_tree,
                      prune_to = dplyr::bind_rows(biosolids_class,
                                                  usgs_class))
my_tree <- bind_entities(my_tree,
                         data = dplyr::bind_rows(
                           biosolids_class, usgs_class),
                         entity_id_col = "CASRN")

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
                                       CASRN,
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
                                       CASRN,
                                       label)) %>%
  dplyr::select(-level_number) %>%
  tidyr::pivot_wider(names_from = level,
                     values_from = label) %>%
  add_terminal_label(entity_id_col = "DTXSID", tax_level_labels = c(chemont_tax_levels,
                                                                    "level12"))

#we really only need the similarity of the terminal classifications.
#which should cut down substantially on the size of the problem.
system.time(
  my_tree_mat2 <- similarity_matrix(tree = my_tree,
                                 nodes1 = my_biosolids_class$terminal_label,
                                 nodes2 = my_usgs_class$terminal_label,
                                 upper_tri = FALSE)
)

#now... create heatmap of similarity of entities only.
my_htmap <- generate_heatmap(row_data = my_biosolids_class,
                 column_data = my_usgs_class,
                 terminal_only = TRUE,
                 tree_object = my_tree,
                 matrix = my_tree_mat2,
                 entity_id_col = "CASRN",
                 row_split = 5L,
                 column_split = 5L)
