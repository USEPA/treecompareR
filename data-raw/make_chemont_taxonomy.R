library(jsonlite)

chemont_df <- jsonlite::fromJSON("http://classyfire.wishartlab.com/tax_nodes.json")
names(chemont_df) <- c("Name", "ID", "Parent_ID")
usethis::use_data(chemont_df, overwrite = TRUE)

chemont_tree <- generate_taxonomy_tree(tax_nodes = chemont_df)
usethis::use_data(chemont_tree, overwrite = TRUE)

chemont_tax_levels <-  c('kingdom', 'superclass',
                                           'class', 'subclass',
                                           'level5', 'level6',
                                           'level7', 'level8',
                                           'level9', 'level10',
                                           'level11')
usethis::use_data(chemont_tax_levels, overwrite = TRUE)
