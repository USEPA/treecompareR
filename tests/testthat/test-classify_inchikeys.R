test_that("classify_inchikeys fails gracefully if network is down",{
  expect_no_error(
    classify_inchikeys(inchikeys = "",
                       network_check_url = "http://httpstat.us/503"),
    regexp = "No\\sinternet\\sconnection\\sor\\ssomething\\swrong\\swith\\sthe\\snetwork"
  )

  expect_no_warning(
    classify_inchikeys(inchikeys = "",
                       network_check_url = "http://httpstat.us/503"),
    regexp = "No\\sinternet\\sconnection\\sor\\ssomething\\swrong\\swith\\sthe\\snetwork"
  )

  expect_message(
    classify_inchikeys(inchikeys = "",
                       network_check_url = "http://httpstat.us/503"),
    regexp = "No\\sinternet\\sconnection\\sor\\ssomething\\swrong\\swith\\sthe\\snetwork"
  )
}

)

test_that("classify_inchikeys fails gracefully with no error if ClassyFire is down",{
  expect_no_error(
    classify_inchikeys(inchikeys = "",
                       classyfire_check_url = "http://httpstat.us/503"),
    regexp = "ClassyFire\\sappears\\sto\\sbe\\sdown"
  )
}
)

test_that("classify_inchikeys fails gracefully with no warning if ClassyFire is down",{
  expect_no_warning(
    classify_inchikeys(inchikeys = "",
                       classyfire_check_url = "http://httpstat.us/503"),
    regexp = "ClassyFire\\sappears\\sto\\sbe\\sdown"
  )
}
)

test_that("classify_inchikeys fails gracefully with a message if ClassyFire is down",{
  expect_message(
    classify_inchikeys(inchikeys = "",
                       classyfire_check_url = "http://httpstat.us/503"),
    regexp = "ClassyFire\\sappears\\sto\\sbe\\sdown"
  )
}
)

test_that("classify_inchikeys works with a chemical that has intermediate nodes", {
  expect_no_error(
    classify_inchikeys(
      inchikeys = "SEMRCUIXRUXGJX-UHFFFAOYSA-N"
    )
  )
}
)

test_that("classify_inchikeys works with a chemical that has no intermediate nodes", {
  expect_no_error(
    classify_inchikeys(
      inchikeys = "PLDWAJLZAAHOGG-UHFFFAOYSA-N"
    )
  )
}
)

test_that("classify_inchikeys handles a bad inchikey", {
  expect_no_error(
    classify_inchikeys(
      inchikeys = "BAD-INCHIKEY-X"
    )
  )
}
)

test_that("classify_inchikeys handles a combination of good and bad inchikeys", {
  expect_no_error(
    classify_inchikeys(
      inchikeys = c("SEMRCUIXRUXGJX-UHFFFAOYSA-N",
                    "PLDWAJLZAAHOGG-UHFFFAOYSA-N",
                    "BAD-INCHIKEY-X",
                    NA_character_,
                    "   ",
                    "\t\r")
    )
  )
}
)

test_that("classify_inchikeys handles all missing inchikeys", {
  expect_no_error(
    classify_inchikeys(
      inchikeys = c(NA_character_,
                    "   ",
                    "\t")
    )
  )
}
)
