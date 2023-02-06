library(jsonlite)

chemont_df <- jsonlite::fromJSON("http://classyfire.wishartlab.com/tax_nodes.json")
names(chemont_df) <- c("Name", "ID", "Parent_ID")
#Correct one error in parentage that we found through trial and error:
#"Ureides" (CHEMONTID:0002022) is listed with parent ID CHEMONTID:0000364,
#but it should be listed with parent ID CHEMONTID:0000517 ("Ureas")
chemont_df[chemont_df$ID %in% "CHEMONTID:0002022",
           "Parent_ID"] <- "CHEMONTID:0000517"

usethis::use_data(chemont_df, overwrite = TRUE)

chemont_tree <- generate_taxonomy_tree(tax_nodes = chemont_df)
usethis::use_data(chemont_tree, overwrite = TRUE)

#Define the list of ChemOnt taxonomy levels
chemont_tax_levels <-  c('kingdom',
                         'superclass',
                         'class',
                         'subclass',
                         'level5',
                         'level6',
                         'level7',
                         'level8',
                         'level9',
                         'level10',
                         'level11')
usethis::use_data(chemont_tax_levels, overwrite = TRUE)
