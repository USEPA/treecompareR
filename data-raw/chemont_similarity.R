#chemont_similarity.R

#run this after running "data-raw/make_chemont_taxonomy.R"
#calculate similarity metrics for all pairs of nodes in the ChemOnt tree
#note that this will take a long time to run!

devtools::load_all("treecompareR") #to get chemont_tree

chemont_jaccard <- similarity_matrix(tree = chemont_tree,
                                     metric = "jaccard")
usethis::use_data(chemont_jaccard, overwrite = TRUE)

chemont_resnik_IC_SVH <- similarity_matrix(tree = chemont_tree,
                                           metric = "resnik")
usethis::use_data(chemont_resnik_IC_SVH, overwrite = TRUE)

chemont_lin_IC_SVH <- similarity_matrix(tree = chemont_tree,
                                        metric = "lin")
usethis::use_data(chemont_lin_IC_SVH, overwrite = TRUE)

chemont_jiangconrath_IC_SVH <- similarity_matrix(tree = chemont_tree,
                                                 metric = "jiang_conrath")
usethis::use_data(chemont_jiangconrath_IC_SVH, overwrite = TRUE)
