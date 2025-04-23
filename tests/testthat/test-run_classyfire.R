test_that("query_structure works with one good SMILES", {
  expect_no_error(
    query_structure(input = "COC1=CC(Br)=CC=C1")
  )
})


test_that("query_structure works with multiple good SMILES", {
  expect_no_error(
    query_structure(
      input = c(
        "COC1=CC(Br)=CC=C1",
        "SCCSCCS",
        "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O"
      )
    )
  )
})

test_that("query_structure handles a bad SMILES", {
  expect_no_error(
    query_structure(
      input = "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|"
    )
  )
}
)

test_that("query_structure handles multiple bad SMILES", {
  expect_no_error(
    query_structure(
      input = c("C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
  )
}
)

test_that("query_structure handles a mixture of good and bad SMILES", {
  expect_no_error(
    query_structure(
      input = c("COC1=CC(Br)=CC=C1",
                "SCCSCCS",
                "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O",
                "CCCCNP(N)(N)=S",
                "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
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
                "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
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
        "BAD-INCHIKEY-X")
    )
  )
}
)

#somethign with a lipidmap term
test_that("classify_structures works when lipidmap and external descriptors present", {
  expect_no_error(
classify_structures(
  input = c("CCCCCCCCCCCCCCCC(O)=O",
            "CCCC(O)=O")
)
)
}
)

#one with a subclass and one without
test_that("query_structure works when one item has a subclass and one does not",{
expect_no_error(
  query_structure(input = c("CC1=NC2=C(C3=CC=CC=C3C=C2)C1(C)C",
                           "CCCCCCCCCCCCCCCC(O)=O"))
)
})

#two smiles without subclasses
test_that("query_structure works when neither item has a subclass",{
  expect_no_error(
    query_structure(input = c("CC1=NC2=C(C3=CC=CC=C3C=C2)C1(C)C",
                               "CCCCNP(N)(N)=S"))
  )
})

