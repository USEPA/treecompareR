library(ctxR)
devtools::load_all()

biosolids <- ctxR::get_chemicals_in_list(list_name = "BIOSOLIDS2022")

bio_class <- biosolids |>
  dplyr::filter(isMarkush %in% 0) |>
  dplyr::pull(smiles) |>
 treecompareR::classify_structures(wait_sec_get = 5)

biosolids <- biosolids |>
  dplyr::left_join(bio_class$classified,
                   by = c(smiles = "identifier_orig"))

usethis::use_data(biosolids, overwrite = TRUE)

biosolids_alternative_parents <- bio_class$alternative_parents
usethis::use_data(biosolids_alternative_parents, overwrite = TRUE)

biosolids_molecular_framework <- bio_class$molecular_framework
usethis::use_data(biosolids_molecular_framework, overwrite = TRUE)

biosolids_substituents <- bio_class$substituents
usethis::use_data(biosolids_substituents, overwrite = TRUE)

biosolids_external_descriptors <- bio_class$external_descriptors
usethis::use_data(biosolids_external_descriptors, overwrite = TRUE)

biosolids_predicted_chebi <- bio_class$predicted_chebi_terms
usethis::use_data(biosolids_predicted_chebi, overwrite = TRUE)

biosolids_predicted_lipidmaps <- bio_class$predicted_lipidmaps_terms
usethis::use_data(biosolids_predicted_lipidmaps, overwrite = TRUE)

#ToxCast invitroDB
toxcast <- ctxR::get_chemicals_in_list(list_name = "ToxCast_invitroDB_v3_0")

toxcast_class <- toxcast |>
  dplyr::filter(isMarkush %in% 0) |>
  dplyr::pull(smiles) |>
  treecompareR::classify_structures(wait_sec_get = 5)

saveRDS(toxcast_class, 'toxcast_class.Rds')

toxcast <- toxcast |>
  dplyr::left_join(toxcast_class$classified,
                   by = c(smiles = "identifier_orig"))

#There was a problem with page 3 of batch 1 (items 200-300).
#returned error 500.
#re-trying these in batches of 10 to isolate the problem.
tc2 <- toxcast |>
       dplyr::filter(isMarkush %in% 0) |>
       dplyr::slice(200:300) |>
      dplyr::pull(smiles) |>
      treecompareR::classify_structures(sizelim = 10, wait_sec_get = 5)

#looks like the problem is in batch 2 here, or items 20-30, or really items 220-230.

#query these one inchikey at a time

tc3 <- toxcast |>
       dplyr::filter(isMarkush %in% 0) |>
       dplyr::slice(220:230) |>
       dplyr::pull(inchikey) |>
       treecompareR::classify_inchikeys(wait_sec_get = 5)

#okay, well... all of those work just fine?

#what can we do to make it fail more gracefully in these cases?
#we can only infer which structures are on which page
#(right now it's 100 per page but the pagination is still as though it were 10)
#so what if getting a single page fails?
#can we keep track and try again? how do we know if we should try again?
#can we infer entities that *should* be on each page and then mark failed pages,
#and let the user re-try?
#what if they change the pagination again?

bad_id <- toxcast_class$classified |>
  dplyr::filter(entity_type %in% "not queried" & exclude %in% FALSE) |>
  dplyr::distinct(identifier_orig) |>
  dplyr::pull(identifier_orig)

#retry these in batches of 10
bad_class <- treecompareR::classify_structures(input = bad_id,
  sizelim = 10,
  wait_sec_get = 5)

#batch 2 is the one that fails with error 500 (items 11-20)
#batch 18 also fails
#batch 25 also fails

bad_id2 <- bad_class$classified |>
  dplyr::filter(entity_type %in% "not queried" & exclude %in% FALSE) |>
  dplyr::distinct(identifier_orig) |>
  dplyr::pull(identifier_orig)

#try these again one at a time?
bad_class2 <- treecompareR::classify_structures(input = bad_id2,
                                               sizelim = 1,
                                               wait_sec_get = 5)

#7, 26, 39 fails
#yet there are 23 of these 49 marked "not queried"
#because a lot of the queries presumably had 0 pages??


usethis::use_data(toxcast, overwrite = TRUE)
