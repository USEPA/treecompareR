

# USGSWATER dataset
# Go to https://comptox.epa.gov/dashboard/chemical-lists/USGSWATER
# Click the "Export" button
# Select the "CSV" option
# Save the result in data-raw/ as USGSWATER.csv
# This was done on 2022-07-13 at 11:27 AM EDT

#Save in rda format in data/
usgs <- read.csv("data-raw/USGSWATER.csv")
#keep only relevant columns
usgs <- usgs[, c("DTXSID",
                 "PREFERRED.NAME",
                 "CASRN",
                 "INCHIKEY",
                 "SMILES")]
USGSWATER <- usgs
usethis::use_data(USGSWATER, overwrite = TRUE)

# BIOSOLIDS2021 dataset
# Go to https://comptox.epa.gov/dashboard/chemical-lists/BIOSOLIDS2021
# Click the "Export" button
# Select the "CSV" option
# Save the result in data-raw/ as BIOSOLIDS2021.csv
# This was done on 2022-07-13 at 11:30 AM EDT
biosolids <- read.csv("data-raw/BIOSOLIDS2021.csv")
#keep only relevant columns
biosolids <- biosolids[, c("DTXSID",
                 "PREFERRED.NAME",
                 "CASRN",
                 "INCHIKEY",
                 "SMILES")]
BIOSOLIDS2021 <- biosolids
usethis::use_data(BIOSOLIDS2021, overwrite = TRUE)

#ClassyFire classified versions of each list
#load package functions for classification
devtools::load_all()

#Biosolids classifications
#First lookup inchikeys
biosolids_inchi <- classify_inchikeys(inchikeys = unique(BIOSOLIDS2021$INCHIKEY))
biosolids_class <- merge(BIOSOLIDS2021,
                         biosolids_inchi,
                         by = "INCHIKEY")
#for anything that could not be looked up by inchikeys,
#query by smiles
biosolids_smiles <- biosolids_class[is.na(biosolids_class$kingdom) &
                                      nzchar(biosolids_class$SMILES) &
                  !is.na(biosolids_class$SMILES), ]
valid_smiles <- sapply(biosolids_smiles$SMILES,
                         webchem::is.smiles)
biosolids_smiles_class <- classify_structures(input = biosolids_smiles[valid_smiles, "SMILES"])

biosolids_smiles_class2 <- dplyr::left_join(x = biosolids_smiles[, 1:5],
y = biosolids_smiles_class,
by = c("SMILES"="structure"))

#rowbind with inchikey classified datta
biosolids_inchi_class <- dplyr::anti_join(x = biosolids_class,
                                          y = biosolids_smiles)

biosolids_smiles_class2 <- biosolids_smiles_class2[names(biosolids_inchi_class)]

#save classified biosolids data
BIOSOLIDS2021_class <- rbind(biosolids_inchi_class, biosolids_smiles_class2)

usethis::use_data(BIOSOLIDS2021_class, overwrite = TRUE)

#USGS classifications
#First lookup inchikeys
usgs_inchi <- classify_inchikeys(inchikeys = unique(USGSWATER$INCHIKEY))
usgs_class <- merge(USGSWATER,
                         usgs_inchi,
                         by = "INCHIKEY")

#for anything that could not be looked up by inchikeys,
#query by smiles
usgs_smiles <- usgs_class[is.na(usgs_class$kingdom) &
                                      nzchar(usgs_class$SMILES) &
                                      !is.na(usgs_class$SMILES), ]
valid_smiles <- sapply(usgs_smiles$SMILES,
                       webchem::is.smiles)
usgs_smiles_class <- classify_structures(input = usgs_smiles[valid_smiles, "SMILES"])

usgs_smiles_class2 <- dplyr::left_join(x = usgs_smiles[, 1:5],
                                            y = usgs_smiles_class,
                                            by = c("SMILES"="structure"))

#rowbind with inchikey classified datta
usgs_inchi_class <- dplyr::anti_join(x = usgs_class,
                                          y = usgs_smiles)

usgs_smiles_class2 <- usgs_smiles_class2[names(usgs_inchi_class)]

#save classified usgs data
USGSWATER_class <- rbind(usgs_inchi_class, usgs_smiles_class2)

#save classified USGS water data
usethis::use_data(USGSWATER_class, overwrite = TRUE)
