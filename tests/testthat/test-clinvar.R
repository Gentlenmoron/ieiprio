test_that(".spdi_a_vcf convierte coordenadas 0-based a VCF", {
  r <- .spdi_a_vcf(c("NC_000023.11:101354147:C:T", "NC_000001.11:999:G:A",
                     "NC_012920.1:10:A:G", "NC_000001.11:5:CA:C", "algo"))
  expect_equal(r[[1]], data.frame(chr = "X", pos = 101354148L, ref = "C", alt = "T"))
  expect_equal(r[[2]]$pos, 1000L)
  expect_equal(r[[3]]$chr, "MT")
  expect_null(r[[4]])
  expect_null(r[[5]])
})

test_that(".clinvar_procesar se queda con SNV de un solo gen", {
  j <- jsonlite::read_json(fixture("clinvar_esummary.json"))$result
  t <- .clinvar_procesar(j[setdiff(names(j), "uids")])
  expect_equal(t$clinvar_id, c(11L, 44L))
  expect_equal(t$gen, c("BTK", "RAG1"))
  expect_equal(t$clasificacion, c("Pathogenic", "Uncertain significance"))
  expect_equal(t$pos[1], 101354148L)
  expect_equal(nrow(.clinvar_procesar(list())), 0)
})

test_that(".clinvar_elegir usa la busqueda y elige el menor ID valido", {
  j <- jsonlite::read_json(fixture("clinvar_esummary.json"))$result
  terminos <- character(0)
  local_mocked_bindings(
    .clinvar_buscar = function(termino, retmax = 300) { terminos <<- c(terminos, termino); c(44L, 11L, 22L) },
    .clinvar_resumen = function(ids) j[as.character(ids)]
  )
  v <- .clinvar_elegir("BTK", "Pathogenic")
  expect_equal(v$clinvar_id, 11L)
  expect_match(terminos[1], "BTK\\[gene\\] AND clinsig_pathogenic\\[Properties\\]")
  expect_error(.clinvar_elegir("BTK", "Pathogenic", excluir = 11L), "No encontre")
})

test_that("ClinVar en vivo", {
  skip_if_not(Sys.getenv("IEIPRIO_EN_VIVO") == "true", "test en vivo desactivado")
  v <- .clinvar_elegir("BTK", "Pathogenic")
  expect_equal(v$gen, "BTK")
  expect_equal(v$chr, "X")
})

test_that("iei_ejemplo_vcf avisa si no esta generado o lo devuelve", {
  ruta <- system.file("extdata", "ejemplo.vcf.gz", package = "ieiprio")
  if (nzchar(ruta)) {
    v <- iei_leer_vcf(iei_ejemplo_vcf(), sexo = "M")
    expect_equal(unique(v$muestra), "EJEMPLO01")
    expect_gt(nrow(v), 3000)
    s <- iei_verificar_sexo(v)
    expect_equal(s$sexo_inferido, "M")
    expect_true(s$concordante)
  } else {
    expect_error(iei_ejemplo_vcf(), "no esta instalado")
  }
})
