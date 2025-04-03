knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(treecompareR)
library(kableExtra)

options(width = 80)
chemont_taxonomy <- generate_taxonomy_tree(tax_nodes = chemont_df)
str(chemont_taxonomy[[1]], width = 60)
str(chemont_taxonomy[[2]], width = 60)

## biosolids <- data.table(chemical_list_biosolids_2022_05_10)
## biosolids[INCHIKEY == '' & CASRN != '', INCHIKEY := {
##   temp = ''
##   attempt <- get_chemical_identifiers(unique(CASRN))
##   if (!is.null(attempt)){
##     temp = attempt@meta$inchikey
##   }
##   ifelse(is.null(temp), '', temp)
## }, by = CASRN]
## biosolids[SMILES == '' & CASRN != '', SMILES := {
##   temp = ''
##   attempt <- get_chemical_identifiers(unique(CASRN))
##   if (!is.null(attempt)){
##     temp = attempt@meta$smiles
##   }
##   ifelse(is.null(temp), '', temp)
## }, by = CASRN]
## 
## biosolids_classified <- classify_datatable(biosolids)
## biosolids_classified <- classify_by_smiles(biosolids_classified)

biosolids_classified <- biosolids_class

## usgswater <- data.table(chemical_list_USGSWATER_2022_05_17)
## usgswater[is.na(INCHIKEY) & !is.na(CASRN), INCHIKEY := {
##   temp = ''
##   attempt <- get_chemical_identifiers(unique(CASRN))
##   if (!is.null(attempt)){
##     temp = attempt@meta$inchikey
##   }
##   ifelse(is.null(temp), '', temp)
## }, by = CASRN]
## usgswater[is.na(SMILES) & !is.na(CASRN), SMILES := {
##   temp = ''
##   attempt <- get_chemical_identifiers(unique(CASRN))
##   if (!is.null(attempt)){
##     temp = attempt@meta$smiles
##   }
##   ifelse(is.null(temp), '', temp)
## }, by = CASRN]
## 
## usgswater_classified <- classify_datatable(usgswater)
## usgswater_classified <- classify_by_smiles(usgswater_classified)

usgswater_classified <- usgs_class

kableExtra::kbl(biosolids_classified[1:5, c('PREFERRED.NAME', 'INCHIKEY', 
                                         'AVERAGE.MASS', 'kingdom',
                                         'superclass', 'class', 
                                         'subclass', 'level5')]) %>%
  kableExtra::kable_styling(full_width = FALSE, bootstrap_options = 'striped', font_size = 9)
#  kbl() %>%
#  kable_styling(bootstrap_options = c("striped", "hover")) %>%
#  kable_classic()

kableExtra::kbl(usgswater_classified[1:5, c('PREFERRED.NAME', 'INCHIKEY',
                                         'AVERAGE.MASS', 'kingdom',
                                         'superclass', 'class',
                                         'subclass', 'level5')])  %>%
  kableExtra::kable_styling(full_width = FALSE, bootstrap_options = 'striped', font_size = 9)
#  kbl() %>%
#  kable_styling(bootstrap_options = c("striped", "hover")) %>%
#  kable_classic()

data_list <- list(biosolids_classified, 
                  usgswater_classified)
names(data_list) <- c('Biosolids', 'USGS Water')
label_bars(data_list)
## fig_1_2 <- label_bars(data_list)
## 
## pdf(file = 'label-bars_fig1.pdf',
##     width = 16,
##     height = 8.27)
## fig_1_2[[1]]
## dev.off()
## 
## pdf(file = 'label-bars_fig2.pdf',
##     width = 16,
##     height = 8.27)
## fig_1_2[[2]]
## dev.off()

ggtree(chemont_tree) + layout_circular()
## fig_3 <- ggtree(chemont_tree) + layout_circular()
## 
## pdf(file = 'chemont-tree_fig3.pdf',
##     width = 16,
##     height = 8.27)
## fig_3
## dev.off()

display_subtree(data_1 = biosolids_classified, 
                name_1 = 'Biosolids')
display_subtree(data_1 = usgswater_classified, 
                name_1 = 'USGS Water')
## fig_4 <- display_subtree(data_1 = biosolids_classified,
##                 name_1 = 'Biosolids')
## fig_5 <- display_subtree(data_1 = usgswater_classified,
##                 name_1 = 'USGS Water')
## 
## pdf(file = 'biosolids-subtree_fig4.pdf',
##     width = 16,
##     height = 8.27)
## fig_4
## dev.off()
## 
## pdf(file = 'usgs-water-subtree_fig5.pdf',
##     width = 16,
##     height = 8.27)
## fig_5
## dev.off()

display_subtree(data_1 = biosolids_classified, 
                data_2 = usgswater_classified, 
                name_1 = 'Biosolids', 
                name_2 = 'USGS Water')

## fig_6 <- display_subtree(data_1 = biosolids_classified,
##                 data_2 = usgswater_classified,
##                 name_1 = 'Biosolids',
##                 name_2 = 'USGS Water')
## 
## pdf(file = 'subtree-overlap_fig6.pdf',
##     width = 16,
##     height = 8.27)
## fig_6
## dev.off()

prune_and_display_subtree(prune_to = biosolids_classified)
prune_and_display_subtree(prune_to = usgswater_classified)
## fig_7 <- prune_and_display_subtree(prune_to = biosolids_classified)
## fig_8 <- prune_and_display_subtree(prune_to = usgswater_classified)
## 
## pdf(file = 'biosolids-prune-subtree_fig7.pdf',
##     width = 16,
##     height = 8.27)
## fig_7
## dev.off()
## 
## pdf(file = 'usgswater-prune-subtree_fig8.pdf',
##     width = 16,
##     height = 8.27)
## fig_8
## dev.off()

data_set_subtrees(data_1 = biosolids_classified, 
                  data_2 = usgswater_classified, 
                  name_1 = 'Biosolids', 
                  name_2 = 'USGS water')
## fig_9_10 <- data_set_subtrees(data_1 = biosolids_classified,
##                   data_2 = usgswater_classified,
##                   name_1 = 'Biosolids',
##                   name_2 = 'USGS water')
## pdf(file = 'biosolids-usgs-shaded-subtree_fig9.pdf',
##     width = 16,
##     height = 8.27)
## fig_9_10[[1]]
## dev.off()
## 
## pdf(file = 'usgs-biosolids-shaded-subtree_fig10.pdf',
##     width = 16,
##     height = 8.27)
## fig_9_10[[2]]
## dev.off()

biosolids_leaf_fraction <- leaf_fraction_subtree(data_1 = biosolids_classified, 
                                                 data_2 = usgswater_classified, 
                                                 name_1 = 'Biosolids', 
                                                 name_2 = 'USGS water')
biosolids_leaf_fraction[[1]]
kableExtra::kbl(head(biosolids_leaf_fraction[[2]]))%>%
  kableExtra::kable_styling(bootstrap_options = 'striped', font_size = 9)
## fig_11 <- leaf_fraction_subtree(data_1 = biosolids_classified,
##                                                  data_2 = usgswater_classified,
##                                                  name_1 = 'Biosolids',
##                                                  name_2 = 'USGS water')
## pdf(file = 'biosolids-leaf-fraction_fig11.pdf',
##     width = 16,
##     height = 8.27)
## fig_11[[1]]
## dev.off()

usgswater_leaf_fraction <- leaf_fraction_subtree(data_1  = usgswater_classified, 
                                                 data_2 = biosolids_classified, 
                                                 name_1 = 'USGS water', 
                                                 name_2 = 'Biosolids')
usgswater_leaf_fraction[[1]]
kableExtra::kbl(usgswater_leaf_fraction[[2]][1:5, ])%>%
  kableExtra::kable_styling(bootstrap_options = 'striped', font_size = 9)
#  kbl() %>%
#  kable_styling(bootstrap_options = c("striped", "hover")) %>%
#  kable_classic()

circ_tree_boxplot(biosolids_classified, 
                  col = 'AVERAGE.MASS', 
                  title = 'Biosolids', 
                  tippoint_boxplot = TRUE, 
                  layers = c('kingdom', 'superclass'))
## fig_12 <- circ_tree_boxplot(biosolids_classified,
##                   col = 'AVERAGE.MASS',
##                   title = 'Biosolids',
##                   tippoint_boxplot = TRUE,
##                   layers = c('kingdom', 'superclass'))
## pdf(file = 'biosolids-circ-boxplot_fig12.pdf',
##     width = 16,
##     height = 8.27)
## fig_12
## dev.off()

circ_tree_boxplot(biosolids_classified, 
                  col = 'AVERAGE.MASS', 
                  title = 'Biosolids', 
                  layers = c('kingdom', 'superclass'))

biosolids_usgs_ht <- generate_heatmap(tree_object = chemont_tree, 
                                      matrix = chemont_jaccard, 
                                      row_data = usgswater_classified, 
                                      column_data = biosolids_classified, 
                                      row_split = 9L, column_split = 9L, 
                                      row_title = 'USGS Water', 
                                      column_title = 'Biosolids', 
                                      name = 'Jaccard Similarity')
biosolids_usgs_ht
## fig_13 <- generate_heatmap(tree_object = chemont_tree,
##                                       matrix = chemont_jaccard,
##                                       row_data = USGSWATER_class,
##                                       column_data = BIOSOLIDS2021_class,
##                                       row_split = 9L, column_split = 9L,
##                                       row_title = 'USGS Water',
##                                       column_title = 'Biosolids',
##                                       name = 'Jaccard Similarity')
## pdf(file = 'heatmap-biosolids-usgs_fig13.pdf',
##     width = 16,
##     height = 8.27)
## fig_13
## dev.off()

generate_tree_cluster(tree = chemont_tree, tree_object = chemont_tree, 
                      htmap = biosolids_usgs_ht, row_cluster = 2, 
                      column_cluster = 2, row_name = 'USGS Water', 
                      column_name = 'Biosolids', isolate_subtree = TRUE)
## fig_14_15 <- generate_tree_cluster(tree = chemont_tree, tree_object = chemont_tree,
##                       htmap = biosolids_usgs_ht, row_cluster = 2,
##                       column_cluster = 2, row_name = 'USGS Water',
##                       column_name = 'Biosolids', isolate_subtree = TRUE)
## pdf(file = 'row2-col2-full_fig14.pdf',
##     width = 16,
##     height = 8.27)
## fig_14_15[[1]]
## dev.off()
## 
## pdf(file = 'row2-col2-prune_fig15.pdf',
##     width = 16,
##     height = 8.27)
## fig_14_15[[2]]
## dev.off()

generate_tree_cluster(tree = chemont_tree, tree_object = chemont_tree, 
                      htmap = biosolids_usgs_ht, row_cluster = 9, 
                      column_cluster = 8, row_name = 'USGS Water', 
                      column_name = 'Biosolids', isolate_subtree = TRUE)

## fig_16_17 <- generate_tree_cluster(tree = chemont_tree, tree_object = chemont_tree,
##                       htmap = biosolids_usgs_ht, row_cluster = 9,
##                       column_cluster = 8, row_name = 'USGS Water',
##                       column_name = 'Biosolids', isolate_subtree = TRUE)
## pdf(file = 'row9-col8-full_fig16.pdf',
##     width = 16,
##     height = 8.27)
## fig_16_17[[1]]
## dev.off()
## 
## pdf(file = 'row9-col8-prune_fig17.pdf',
##     width = 16,
##     height = 8.27)
## fig_16_17[[2]]
## dev.off()
