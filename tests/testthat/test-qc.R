vcf_x <- function(n_het, n_hom, sexo = NULL) {
  n <- n_het + n_hom
  x <- tibble::tibble(muestra = "S1", chr = "X", pos = 5e6 + seq_len(n), ref = "A", alt = "G",
                      genotipo = c(rep("het", n_het), rep("hom_alt", n_hom)),
                      dp = 30, gq = 99, ab = c(rep(0.5, n_het), rep(1, n_hom)),
                      filtro = "PASS", build = "GRCh38")
  attr(x, "build") <- "GRCh38"
  attr(x, "ieiprio_sexo") <- c(S1 = if (is.null(sexo)) NA_character_ else sexo)
  x
}

test_that("estima el sexo desde el X", {
  expect_equal(iei_verificar_sexo(vcf_x(1, 40))$sexo_inferido, "M")
  expect_equal(iei_verificar_sexo(vcf_x(20, 20))$sexo_inferido, "F")
  expect_equal(iei_verificar_sexo(vcf_x(7, 33))$sexo_inferido, "indeterminado")
  expect_true(is.na(iei_verificar_sexo(vcf_x(2, 5))$sexo_inferido))
})

test_that("avisa cuando el sexo declarado no coincide", {
  expect_warning(r <- iei_verificar_sexo(vcf_x(20, 20, "M")), "no coincide")
  expect_false(r$concordante)
  expect_true(iei_verificar_sexo(vcf_x(1, 40, "M"))$concordante)
})

test_that("ignora la region pseudoautosomica", {
  x <- vcf_x(30, 0)
  x$pos <- 50000 + seq_len(nrow(x))
  expect_true(is.na(iei_verificar_sexo(x)$sexo_inferido))
})

test_that("iei_leer_vcf guarda sexo y huella del VCF", {
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  expect_equal(attr(v, "ieiprio_sexo"), c(P01 = "M", P02 = "F"))
  expect_equal(attr(v, "ieiprio_vcf")$archivo, "mini.vcf")
  expect_match(attr(v, "ieiprio_vcf")$md5, "^[0-9a-f]{32}$")
})

test_that("iei_filtrar calcula el control de calidad por muestra", {
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  f <- iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba())
  qc <- attr(f, "ieiprio_qc")
  expect_named(qc, c("P01", "P02"))
  expect_equal(qc$P01$n_variantes, 10)
  expect_equal(qc$P02$sin_llamar, 1)
  expect_equal(sum(qc$P01$hist_dp$n), sum(!is.na(v$dp[v$muestra == "P01"])))
  expect_equal(qc$P01$sexo_declarado, "M")
  expect_true(is.na(qc$P01$sexo_inferido))
  p <- attr(f, "ieiprio_parametros")
  expect_equal(p$min_dp, 10)
  expect_match(p$genes, "lista propia")
  expect_equal(attr(f, "ieiprio_vcf")$archivo, "mini.vcf")
})

test_that("los atributos llegan hasta la priorizacion", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999")
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  f <- iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba())
  a <- suppressMessages(iei_anotar(f, cache = FALSE))
  p <- iei_priorizar(iei_filtrar_frecuencia(a))
  for (n in c("ieiprio_qc", "ieiprio_parametros", "ieiprio_vcf", "ieiprio_sexo",
              "ieiprio_bitacora", "ieiprio_max_af", "ieiprio_pesos")) {
    expect_false(is.null(attr(p, n)), info = n)
  }
  expect_equal(attr(p, "vep_release"), "999")
  expect_equal(p$proteina_id[p$pos == 101370000], "ENSP_SINTETICO_1")
  expect_equal(p$pos_proteina[p$pos == 101370000], 500)
})
