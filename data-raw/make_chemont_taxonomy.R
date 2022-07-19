library(jsonlite)

chemont_df <- jsonlite::fromJSON("http://classyfire.wishartlab.com/tax_nodes.json")
names(chemont_df) <- c("Name", "ID", "Parent_ID")
usethis::use_data(chemont_df, overwrite = TRUE)

chemont_tree <- generate_taxonomy_tree(tax_nodes = chemont_df)
usethis::use_data(chemont_tree, overwrite = TRUE)

