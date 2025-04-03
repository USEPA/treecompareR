## ----include = FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------------------------------------
library(treecompareR)
library(kableExtra)


## ----chemont-tree, fig.cap="chemont-tree"-------------------------------------------------------------------
options(width = 80)
chemont_taxonomy <- generate_taxonomy_tree(tax_nodes = chemont_df)
str(chemont_taxonomy[[1]], width = 60)
str(chemont_taxonomy[[2]], width = 60)


## ----display_classifications--------------------------------------------------------------------------------
kableExtra::kbl(biosolids_class[1:5, c('PREFERRED.NAME', 'INCHIKEY', 
                                         'AVERAGE.MASS', 'kingdom',
                                         'superclass', 'class', 
                                         'subclass', 'level5')]) %>%
  kableExtra::kable_styling(full_width = FALSE, bootstrap_options = 'striped', font_size = 9)
#  kbl() %>%
#  kable_styling(bootstrap_options = c("striped", "hover")) %>%
#  kable_classic()

kableExtra::kbl(usgs_class[1:5, c('PREFERRED.NAME', 'INCHIKEY',
                                         'AVERAGE.MASS', 'kingdom',
                                         'superclass', 'class',
                                         'subclass', 'level5')])  %>%
  kableExtra::kable_styling(full_width = FALSE, bootstrap_options = 'striped', font_size = 9)
#  kbl() %>%
#  kable_styling(bootstrap_options = c("striped", "hover")) %>%
#  kable_classic()


## ----label-bars, fig.cap="label-bars",  fig.align='center', fig.dim=c(6,4)----------------------------------
data_list <- list(biosolids_class, 
                  usgs_class)
names(data_list) <- c('Biosolids', 'USGS Water')
label_bars(data_list)

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
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


## ----chemont-tree-plot, fig.cap="chemont-tree-plot"---------------------------------------------------------
ggtree(chemont_tree) + layout_circular()

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
## fig_3 <- ggtree(chemont_tree) + layout_circular()
## 
## pdf(file = 'chemont-tree_fig3.pdf',
##     width = 16,
##     height = 8.27)
## fig_3
## dev.off()


## ----subtree-plots, fig.cap="subtree-plots",  fig.align='center', fig.dim=c(6,4)----------------------------
display_subtree(data_1 = biosolids_class, 
                name_1 = 'Biosolids')
display_subtree(data_1 = usgs_class, 
                name_1 = 'USGS Water')

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
## fig_4 <- display_subtree(data_1 = biosolids_class,
##                 name_1 = 'Biosolids')
## fig_5 <- display_subtree(data_1 = usgs_class,
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


## ----subtree-plots-overlap, fig.cap="subtree-plots-overlap", fig.align='center', fig.dim=c(6,4)-------------
display_subtree(data_1 = biosolids_class, 
                data_2 = usgs_class, 
                name_1 = 'Biosolids', 
                name_2 = 'USGS Water')


## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
## fig_6 <- display_subtree(data_1 = biosolids_class,
##                 data_2 = usgs_class,
##                 name_1 = 'Biosolids',
##                 name_2 = 'USGS Water')
## 
## pdf(file = 'subtree-overlap_fig6.pdf',
##     width = 16,
##     height = 8.27)
## fig_6
## dev.off()


## ----pruned-subtree,  fig.align='center', fig.dim=c(12,12), out.height=600, out.width=600-------------------
prune_and_display_subtree(prune_to = biosolids_class)
prune_and_display_subtree(prune_to = usgs_class)

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
## fig_7 <- prune_and_display_subtree(prune_to = biosolids_class)
## fig_8 <- prune_and_display_subtree(prune_to = usgs_class)
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


## ----pruned-subtree-overlap,  fig.align='center', fig.dim=c(12,12), out.height=600, out.width=600-----------
data_set_subtrees(data_1 = biosolids_class, 
                  data_2 = usgs_class, 
                  name_1 = 'Biosolids', 
                  name_2 = 'USGS water')

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
## fig_9_10 <- data_set_subtrees(data_1 = biosolids_class,
##                   data_2 = usgs_class,
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


## ----biosolids-leaf-fraction,  fig.align='center', fig.dim=c(12,8), out.height=400, out.width=600-----------

  biosolids_leaf_fraction <- leaf_fraction_subtree(data_1 = biosolids_class, 
                                                 data_2 = usgs_class, 
                                                 name_1 = 'Biosolids', 
                                                 name_2 = 'USGS water')

print(biosolids_leaf_fraction[[1]])

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
## fig_11 <- leaf_fraction_subtree(data_1 = biosolids_class,
##                                                  data_2 = usgs_class,
##                                                  name_1 = 'Biosolids',
##                                                  name_2 = 'USGS water')
## pdf(file = 'biosolids-leaf-fraction_fig11.pdf',
##     width = 16,
##     height = 8.27)
## fig_11[[1]]
## dev.off()


## ----usgs-leaf-fraction,  fig.align='center', fig.dim=c(12,8), out.height=400, out.width=600----------------

  usgswater_leaf_fraction <- leaf_fraction_subtree(data_1  = usgs_class, 
                                                 data_2 = biosolids_class, 
                                                 name_1 = 'USGS water', 
                                                 name_2 = 'Biosolids')

usgswater_leaf_fraction[[1]]
kableExtra::kbl(usgswater_leaf_fraction[[2]][1:5, ])%>%
  kableExtra::kable_styling(bootstrap_options = 'striped', font_size = 9)
#  kbl() %>%
#  kable_styling(bootstrap_options = c("striped", "hover")) %>%
#  kable_classic()


## ----circ-tree-biosolids-mass-tip-color, fig.align='center', fig.dim=c(12,8), out.height=400, out.width=600----
circ_tree_boxplot(biosolids_class, 
                  col = 'AVERAGE.MASS', 
                  title = 'Biosolids', 
                  tippoint_boxplot = TRUE, 
                  layers = c('kingdom', 'superclass'))

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
## fig_12 <- circ_tree_boxplot(biosolids_class,
##                   col = 'AVERAGE.MASS',
##                   title = 'Biosolids',
##                   tippoint_boxplot = TRUE,
##                   layers = c('kingdom', 'superclass'))
## pdf(file = 'biosolids-circ-boxplot_fig12.pdf',
##     width = 16,
##     height = 8.27)
## fig_12
## dev.off()


## ----circ-tree-biosolids-mass-no-tip-color, fig.align='center', fig.dim=c(12,8), out.height=400, out.width=600----
circ_tree_boxplot(biosolids_class, 
                  col = 'AVERAGE.MASS', 
                  title = 'Biosolids', 
                  layers = c('kingdom', 'superclass'))


## ----generate-similarity-chemont, eval = FALSE--------------------------------------------------------------
## #the following will take a long time to run!
## #however, this is how the similarity matrixes can be generated using treecompareR functions.
## 
## chemont_jaccard <- similarity_matrix(tree = chemont_tree,
##                                               metric = "jaccard")
## 
## chemont_resnik_IC_SVH <- similarity_matrix(tree = chemont_tree,
##                                                    metric = "resnik")
## 
## chemont_lin_IC_SVH <- similarity_matrix(tree = chemont_tree,
##                                                     metric = "lin")
## 
## chemont_jiangconrath_IC_SVH <- similarity_matrix(tree = chemont_tree,
##                                                     metric = "jiang_conrath")


## ----heatmap-biosolids-usgs, fig.align='center', fig.dim=c(6,4)---------------------------------------------
biosolids_usgs_ht <- generate_heatmap(tree_object = chemont_tree, 
                                      matrix = chemont_jaccard, 
                                      row_data = usgs_class, 
                                      column_data = biosolids_class, 
                                      row_split = 9L, column_split = 9L, 
                                      row_title = 'USGS Water', 
                                      column_title = 'Biosolids', 
                                      name = 'Jaccard Similarity')
biosolids_usgs_ht

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
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


## ----rc-cluster-2-2, fig.align='center', fig.dim=c(18,12), out.height=480, out.width=720--------------------
generate_tree_cluster(tree = chemont_tree, tree_object = chemont_tree, 
                      htmap = biosolids_usgs_ht, row_cluster = 2, 
                      column_cluster = 2, row_name = 'USGS Water', 
                      column_name = 'Biosolids', isolate_subtree = TRUE)

## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
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


## ----rc-cluster-9-8, fig.align='center', fig.dim=c(18,10), out.height=350, out.width=630--------------------
generate_tree_cluster(tree = chemont_tree, tree_object = chemont_tree, 
                      htmap = biosolids_usgs_ht, row_cluster = 9, 
                      column_cluster = 8, row_name = 'USGS Water', 
                      column_name = 'Biosolids', isolate_subtree = TRUE)


## ----include=FALSE, eval=FALSE------------------------------------------------------------------------------
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

