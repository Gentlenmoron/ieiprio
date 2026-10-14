test_that(".predictores lee los plugins de VEP", {
  t <- list(alphamissense = list(am_pathogenicity = 0.91, am_class = "likely_pathogenic"),
            revel = 0.82, spliceai = list(DS_AG = 0.01, DS_DL = 0.73, SYMBOL = "BTK"))
  p <- .predictores(t)
  expect_equal(p$alphamissense, 0.91)
  expect_equal(p$am_clase, "likely_pathogenic")
  expect_equal(p$revel, 0.82)
  expect_equal(p$spliceai, 0.73)
  vacio <- .predictores(list())
  expect_true(all(is.na(unlist(vacio[c("alphamissense", "revel", "spliceai")]))))
})

test_that("la anotacion agrega los predictores", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999")
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  a <- suppressMessages(iei_anotar(iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba()),
                                   cache = FALSE))
  r1 <- a[a$muestra == "P01" & a$pos == 101360000, ]
  expect_equal(r1$revel, 0.82)
  expect_equal(r1$am_clase, "likely_pathogenic")
  expect_equal(a$spliceai[a$pos == 101370000], 0.73)
})

test_that(".evaluar_predictores aplica PP3 y BP4", {
  x <- tibble::tibble(revel = c(0.8, 0.1, NA, 0.5, NA),
                      alphamissense = c(NA, 0.1, NA, 0.4, NA),
                      spliceai = c(NA, 0.02, 0.6, NA, NA))
  r <- .evaluar_predictores(x)
  expect_equal(r$favor, c(TRUE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(r$contra, c(FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_match(r$texto_favor[1], "REVEL 0.80")
  expect_match(r$texto_favor[3], "SpliceAI 0.60")
})

test_that("los predictores entran al puntaje y a ACMG", {
  x <- tibble::tibble(muestra = "P1", gen = c("BTK", "BTK"), genotipo = "hemi",
                      impacto = "MODERATE", consecuencia = "missense_variant",
                      af_gnomad = NA_real_, clinvar = NA_character_,
                      revel = c(0.9, 0.05), alphamissense = c(0.95, 0.1), spliceai = NA_real_)
  p <- iei_priorizar(x)
  favor <- p[p$revel == 0.9, ]; contra <- p[p$revel == 0.05, ]
  expect_equal(favor$p_predictores, 1)
  expect_equal(contra$p_predictores, -1)
  expect_match(favor$acmg, "PP3")
  expect_match(contra$acmg, "BP4")
  expect_equal(favor$puntaje - contra$puntaje, 2)
  expect_equal(favor$p_clinvar + favor$p_impacto + favor$p_frecuencia + favor$p_predictores +
                 favor$p_herencia + favor$p_fenotipo, favor$puntaje)
})

test_that("tablas sin columnas de predictores siguen funcionando", {
  x <- tibble::tibble(muestra = "P1", gen = "BTK", genotipo = "hemi", impacto = "HIGH",
                      consecuencia = "stop_gained", af_gnomad = NA_real_, clinvar = "pathogenic")
  p <- iei_priorizar(x)
  expect_equal(p$p_predictores, 0)
  expect_equal(p$puntaje, 10)
})
