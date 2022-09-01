

# USGSWATER dataset
# Go to https://comptox.epa.gov/dashboard/chemical-lists/USGSWATER
# Click the "Export" button
# Select the "CSV" option
# Save the result in data-raw/ as USGSWATER.csv
# This was done on 2022-07-13 at 11:27 AM EDT

#Save in rda format in data/
usgs <- read.csv("data-raw/USGSWATER.csv",
                 check.names = FALSE)
usgs <- setNames(usgs,
                 gsub(pattern = "#",
                      replacement = "Num",
                      x = names(usgs),
                      fixed = TRUE))
usgs <- setNames(usgs,
                 gsub(pattern = "%",
                      replacement = "Pct",
                      x = names(usgs),
                      fixed = TRUE))
usgs <- setNames(usgs,
                 make.names(names(usgs)))
usethis::use_data(usgs, overwrite = TRUE)

# BIOSOLIDS2021 dataset
# Go to https://comptox.epa.gov/dashboard/chemical-lists/BIOSOLIDS2021
# Click the "Export" button
# Select the "CSV" option
# Save the result in data-raw/ as BIOSOLIDS2021.csv
# This was done on 2022-07-13 at 11:30 AM EDT
biosolids <- read.csv("data-raw/BIOSOLIDS2021.csv",
                      check.names = FALSE)
biosolids <- setNames(biosolids,
                 gsub(pattern = "#",
                      replacement = "Num",
                      x = names(biosolids),
                      fixed = TRUE))
biosolids <- setNames(biosolids,
                 gsub(pattern = "%",
                      replacement = "Pct",
                      x = names(biosolids),
                      fixed = TRUE))
biosolids <- setNames(biosolids,
                 make.names(names(biosolids)))
usethis::use_data(biosolids, overwrite = TRUE)

#ClassyFire classified versions of each list
#load package functions for classification
devtools::load_all()

#Biosolids classifications
#First lookup inchikeys
biosolids_inchi <- classify_inchikeys(inchikeys = unique(biosolids$INCHIKEY))
biosolids_class <- merge(biosolids,
                         biosolids_inchi,
                         by = "INCHIKEY")
#for anything that could not be looked up by inchikeys,
#query by smiles
biosolids_smiles <- biosolids_class[is.na(biosolids_class$kingdom) &
                                      nzchar(biosolids_class$SMILES) &
                  !is.na(biosolids_class$SMILES), ]
valid_smiles <- sapply(biosolids_smiles$SMILES,
                         webchem::is.smiles)
biosolids_smiles_class <- classify_structures(input = biosolids_smiles[valid_smiles, "SMILES"],
                                              queued_wait = 10,
                                              processing_wait_per_input = 3)

biosolids_smiles_class2 <- dplyr::left_join(x = biosolids_smiles,
y = biosolids_smiles_class[, c("structure",
                               chemont_tax_levels)],
by = c("SMILES"="structure",
chemont_tax_levels)
)

#rowbind with inchikey classified datta
biosolids_inchi_class <- dplyr::anti_join(x = biosolids_class,
                                          y = biosolids_smiles)

biosolids_smiles_class2 <- biosolids_smiles_class2[names(biosolids_inchi_class)]

#save classified biosolids data
biosolids_class <- rbind(biosolids_inchi_class, biosolids_smiles_class2)
biosolids_class <- add_terminal_label(biosolids_class)

usethis::use_data(biosolids_class, overwrite = TRUE)

#USGS classifications
#First lookup inchikeys
usgs_inchi <- classify_inchikeys(inchikeys = unique(usgs$INCHIKEY))
usgs_class <- merge(usgs,
                         usgs_inchi,
                         by = "INCHIKEY")

#for anything that could not be looked up by inchikeys,
#query by smiles
usgs_smiles <- usgs_class[is.na(usgs_class$kingdom) &
                                      nzchar(usgs_class$SMILES) &
                                      !is.na(usgs_class$SMILES), ]
valid_smiles <- sapply(usgs_smiles$SMILES,
                       webchem::is.smiles)
usgs_smiles_class <- classify_structures(input = usgs_smiles[valid_smiles, "SMILES"],
                                         queued_wait = 10,
                                         processing_wait_per_input = 3)

usgs_smiles_class2 <- dplyr::left_join(x = usgs_smiles,
                                       y = usgs_smiles_class[, c("structure",
                                            chemont_tax_levels)],
                                       by = c("SMILES"="structure",
                                              chemont_tax_levels)
)

#rowbind with inchikey classified datta
usgs_inchi_class <- dplyr::anti_join(x = usgs_class,
                                          y = usgs_smiles)

usgs_smiles_class2 <- usgs_smiles_class2[names(usgs_inchi_class)]

#save classified usgs data
usgs_class <- rbind(usgs_inchi_class, usgs_smiles_class2)
usgs_class <- add_terminal_label(usgs_class)

#save classified USGS water data
usethis::use_data(usgs_class, overwrite = TRUE)
