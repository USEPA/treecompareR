test_that("classify_structures fails gracefully if network is down",{
  expect_no_error(
    classify_structures(input = "CCC",
                        check_avail = TRUE,
                       network_check_url = "http://httpstat.us/503"),
    regexp = "No\\sinternet\\sconnection\\sor\\ssomething\\swrong\\swith\\sthe\\snetwork"
  )

  expect_no_warning(
    classify_structures(input = "CCC",
                        check_avail = TRUE,
                        network_check_url = "http://httpstat.us/503"),
    regexp = "No\\sinternet\\sconnection\\sor\\ssomething\\swrong\\swith\\sthe\\snetwork"
  )

  expect_message(
    classify_structures(input = "CCC",
                        check_avail = TRUE,
                        network_check_url = "http://httpstat.us/503"),
    regexp = "No\\sinternet\\sconnection\\sor\\ssomething\\swrong\\swith\\sthe\\snetwork"
  )
}
)

test_that("classify_structures fails gracefully if ClassyFire is down",{
  expect_no_error(
    classify_structures(input = "CCC",
                        check_avail = TRUE,
                       classyfire_check_url = "http://httpstat.us/503"),
    regexp = "ClassyFire\\sappears\\sto\\sbe\\sdown"
  )

  expect_no_warning(
    classify_structures(input = "CCC",
                        check_avail = TRUE,
                        classyfire_check_url = "http://httpstat.us/503"),
    regexp = "ClassyFire\\sappears\\sto\\sbe\\sdown"
  )

  expect_message(
    classify_structures(input = "CCC",
                        check_avail = TRUE,
                        classyfire_check_url = "http://httpstat.us/503"),
    regexp = "ClassyFire\\sappears\\sto\\sbe\\sdown"
  )
}
)

test_that("classify_structures works with one good SMILES", {
  expect_no_error(
    classify_structures(input = "COC1=CC(Br)=CC=C1")
  )
})

test_that("classify_structures works with multiple good SMILES", {
  expect_no_error(
    classify_structures(
      input = c(
        "COC1=CC(Br)=CC=C1",
        "SCCSCCS",
        "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O"
      )
    )
  )
})

test_that("classify_structures handles a bad SMILES", {
  expect_no_error(
    classify_structures(
      input = "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|"
    )
  )
}
)

test_that("classify_structures handles multiple bad SMILES", {
  expect_no_error(
    classify_structures(
      input = c("C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
  )
}
)

test_that("classify_structures handles a mixture of good and bad SMILES", {
  expect_no_error(
    classify_structures(
      input = c("COC1=CC(Br)=CC=C1",
                "SCCSCCS",
                "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O",
                "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
  )
}
)

test_that("classify_structures handles a mixture of good, bad, and NA/blank SMILES", {
  expect_no_error(
    classify_structures(
      input = c("COC1=CC(Br)=CC=C1",
                "SCCSCCS",
                NA_character_,
                "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O",
                "  ",
                "CCCCCCCCCCCCCCCC(O)=O",
                "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
  )
}
)

test_that("classify_structures works when lipidmap and external descriptors present", {
  expect_no_error(
    classify_structures(
      input = c("CCCCCCCCCCCCCCCC(O)=O",
                "CCCC(O)=O")
    )
  )
}
)

test_that("classify_structures handles all NA/blank SMILES", {
  expect_no_error(
    classify_structures(
      input = c(NA_character_,
                "  ",
                "\t",
                "\n",
                "\r")
    )
  )
}
)

test_that("classify_structures handles batch mode", {
  expect_no_error(
    classify_structures(
      input = c("COC1=CC(Br)=CC=C1",
                "SCCSCCS",
                NA_character_,
                "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O",
                "  ",
                "CCCCCCCCCCCCCCCC(O)=O",
                "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|"),
      sizelim = 3
    )
  )
}
)
