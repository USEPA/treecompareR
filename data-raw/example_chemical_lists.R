
setwd("data-raw/")
# USGSWATER dataset
# Go to https://comptox.epa.gov/dashboard/chemical-lists/USGSWATER
# Click the "Export" button
# Select the "CSV" option
# Save the result in data-raw/ as USGSWATER.csv
# This was done on 2022-07-13 at 11:27 AM EDT

#Save in rda format in data/
usgs <- read.csv("USGSWATER.csv")
#keep only relevant columns
usgs <- usgs[, c("DTXSID",
                 "PREFERRED.NAME",
                 "CASRN",
                 "INCHIKEY",
                 "SMILES")]
save(usgs, file = "../data/USGSWATER.RData")

# BIOSOLIDS2021 dataset
# Go to https://comptox.epa.gov/dashboard/chemical-lists/BIOSOLIDS2021
# Click the "Export" button
# Select the "CSV" option
# Save the result in data-raw/ as BIOSOLIDS2021.csv
# This was done on 2022-07-13 at 11:30 AM EDT
biosolids <- read.csv("BIOSOLIDS2021.csv")
#keep only relevant columns
biosolids <- biosolids[, c("DTXSID",
                 "PREFERRED.NAME",
                 "CASRN",
                 "INCHIKEY",
                 "SMILES")]
save(biosolids, file = "../data/BIOSOLIDS2021.RData")

#ClassyFire classified versions of each list
#load package function classify_inchikeys()
devtools::load_all("../")

biosolids_classified <- classify_inchikeys(inchikeys = biosolids$INCHIKEY)

usgs_classified <- classify_inchikeys(inchikeys = usgs$INCHIKEY)


