test_that(".hpo_procesar filtra genes y quita duplicados", {
  h <- .hpo_procesar(fixture("hpo_mini.tsv"), c("BTK", "RAG1"))
  expect_named(h, c("gen", "hpo_id", "hpo_nombre"))
  expect_equal(nrow(h), 2)
  expect_false("NAT2" %in% h$gen)
})

test_that(".hpo_procesar avisa si cambia el formato", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("a\tb\n1\t2", tmp)
  expect_error(.hpo_procesar(tmp, "BTK"), "columnas")
})

test_that("la URL de HPO apunta al release", {
  expect_match(.hpo_url("v2026-09-01"), "download/v2026-09-01/genes_to_phenotype.txt$")
})
