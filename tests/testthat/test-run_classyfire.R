test_that("query_classyfire works with one good SMILES", {
  expect_no_error(
    query_classyfire(input = "COC1=CC(Br)=CC=C1")
  )
})


test_that("query_classyfire works with multiple good SMILES", {
  expect_no_error(
    query_classyfire(
      input = c(
        "COC1=CC(Br)=CC=C1",
        "SCCSCCS",
        "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O"
      )
    )
  )
})

test_that("query_classyfire handles a bad SMILES", {
  expect_no_error(
    query_classyfire(
      input = "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|"
    )
  )
}
)

test_that("query_classyfire handles multiple bad SMILES", {
  expect_no_error(
    query_classyfire(
      input = c("C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
  )
}
)

test_that("query_classyfire handles a mixture of good and bad SMILES", {
  expect_no_error(
    query_classyfire(
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

test_that("classify_structure works with one good SMILES", {
  expect_no_error(
    classify_structure(input = "COC1=CC(Br)=CC=C1")
  )
})


test_that("classify_structure works with multiple good SMILES", {
  expect_no_error(
    classify_structure(
      input = c(
        "COC1=CC(Br)=CC=C1",
        "SCCSCCS",
        "[K+].[K+].[O-]S(=O)(=O)OOS([O-])(=O)=O"
      )
    )
  )
})

test_that("classify_structure handles a bad SMILES", {
  expect_no_error(
    classify_structure(
      input = "C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|"
    )
  )
}
)

test_that("classify_structure handles multiple bad SMILES", {
  expect_no_error(
    classify_structure(
      input = c("C*.CSC1=NC=CN=C1 |c:5,7,t:3,lp:3:2,5:1,8:1,m:1:9.7.6|",
                "[Cl-].C*.CC(C)(C)CC(C)(C)C1=CC=C(OCCOCC[N+](C)(C)CC2=CC=CC=C2)C=C1 |c:25,27,30,t:9,11,23,lp:0:4,15:2,18:2,m:2:13.12|",
                "[*]OC(=O)C=C |$_R1;;;;;$,lp:1:2,3:2,RG:_R1={CC(O)C* |$;;;;_AP1$,lp:2:2|},{CC(*)CO |$;;_AP1;;$,lp:4:2|}|")
    )
  )
}
)

test_that("classify_structure handles a mixture of good and bad SMILES", {
  expect_no_error(
    classify_structure(
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
