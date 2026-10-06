withr::local_options(list(rlib_message_verbosity = "quiet"))

filtrado_mini <- function() {
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba())
}
mapa_prueba <- c(
  "X 101360000 . C T . . ." = "GRCh38_X_101360000_C_T",
  "X 101370000 . A G . . ." = "GRCh38_X_101370000_A_G",
  "1 1000 . A G . . ." = "GRCh38_1_1000_A_G"
)

test_that(".vep_aplanar extrae transcritos, ClinVar, frecuencia y rsID", {
  t <- .vep_aplanar(vep_sintetico(), mapa_prueba)
  a <- t[t$clave == "GRCh38_X_101360000_C_T", ]
  expect_equal(nrow(a), 3)
  expect_equal(unique(a$rsid), "rs999001")
  expect_equal(unique(a$clinvar), "pathogenic;likely_pathogenic")
  expect_equal(unique(a$af_gnomad), 1e-05)
})

test_that("ignora transcritos de otro alelo", {
  t <- .vep_aplanar(vep_sintetico(), mapa_prueba)
  b <- t[t$clave == "GRCh38_X_101370000_A_G", ]
  expect_equal(nrow(b), 1)
  expect_equal(b$consecuencia, "stop_gained&splice_region_variant")
  expect_true(is.na(b$af_gnomad))
  expect_true(is.na(b$clinvar))
})

test_that("variantes intergenicas quedan con una fila sin transcrito", {
  t <- .vep_aplanar(vep_sintetico(), mapa_prueba)
  c1 <- t[t$clave == "GRCh38_1_1000_A_G", ]
  expect_equal(nrow(c1), 1)
  expect_true(is.na(c1$transcrito))
  expect_equal(c1$consecuencia, "intergenic_variant")
  expect_equal(c1$clinvar, "benign")
  expect_equal(c1$af_gnomad, 0.3)
})

test_that(".vep_elegir prefiere gen, luego MANE, luego canonico", {
  t <- .vep_aplanar(vep_sintetico(), mapa_prueba)
  e <- .vep_elegir("GRCh38_X_101360000_C_T", "BTK", t)
  expect_equal(e$transcrito, "ENST_SINTETICO_1")
  expect_equal(e$mane, "NM_SINTETICO.1")
  e2 <- .vep_elegir("GRCh38_X_101360000_C_T", "GENOTRO", t)
  expect_equal(e2$transcrito, "ENST_OTRO")
  e3 <- .vep_elegir("NO_EXISTE", NA, t)
  expect_true(is.na(e3$transcrito))
})

test_that("iei_anotar agrega columnas y conserva atributos", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999")
  f <- filtrado_mini()
  a <- iei_anotar(f, cache = FALSE)
  expect_equal(nrow(a), nrow(f))
  expect_true(all(c("transcrito", "mane", "consecuencia", "impacto", "hgvs_c",
                    "hgvs_p", "af_gnomad", "clinvar") %in% names(a)))
  expect_equal(a$hgvs_p[a$pos == 101360000][1], "ENSP_SINTETICO_1:p.Arg334Trp")
  expect_equal(a$impacto[a$pos == 101370000], "HIGH")
  expect_equal(a$id[a$pos == 101360000][1], "rs999001")
  expect_equal(attr(a, "vep_release"), "999")
  expect_equal(attr(a, "build"), "GRCh38")
  expect_equal(iei_bitacora(a), iei_bitacora(f))
})

test_that("la cache evita consultar de nuevo", {
  dir <- withr::local_tempdir()
  llamadas <- 0
  local_mocked_bindings(
    .vep_consultar = function(entrada, build) { llamadas <<- llamadas + 1; vep_falso(entrada, build) },
    .vep_release = function(build) "999"
  )
  f <- filtrado_mini()
  a1 <- iei_anotar(f, cache = dir)
  a2 <- iei_anotar(f, cache = dir)
  expect_equal(llamadas, 1)
  expect_equal(a1$hgvs_c, a2$hgvs_c)
  expect_equal(attr(a2, "vep_release"), "999")
  expect_equal(length(iei_limpiar_cache(dir)), 1)
  a3 <- iei_anotar(f, cache = dir)
  expect_equal(llamadas, 2)
})

test_that("divide en lotes", {
  tamanos <- integer(0)
  local_mocked_bindings(
    .vep_consultar = function(entrada, build) { tamanos <<- c(tamanos, length(entrada)); vep_falso(entrada, build) },
    .vep_release = function(build) "999"
  )
  iei_anotar(filtrado_mini(), cache = FALSE, lote = 1)
  expect_equal(tamanos, c(1L, 1L))
})

test_that("variantes sin respuesta quedan en NA y una tabla vacia funciona", {
  local_mocked_bindings(.vep_consultar = function(entrada, build) list(),
                        .vep_release = function(build) NA_character_)
  a <- iei_anotar(filtrado_mini(), cache = FALSE)
  expect_true(all(is.na(a$consecuencia)))
  vacio <- iei_anotar(filtrado_mini()[0, ], cache = FALSE)
  expect_equal(nrow(vacio), 0)
  expect_true("clinvar" %in% names(vacio))
})

test_that("valida entradas", {
  expect_error(iei_anotar(data.frame(a = 1)), "iei_leer_vcf")
  expect_error(iei_anotar(filtrado_mini(), lote = 500), "lote")
})

test_that("iei_filtrar_frecuencia usa el umbral segun herencia", {
  x <- tibble::tibble(gen = c("BTK", "BTK", "RAG1", "RAG1", "BTK"),
                      af_gnomad = c(5e-5, 5e-4, 5e-3, 0.05, NA))
  attr(x, "ieiprio_bitacora") <- tibble::tibble(paso = "Inicial", filas = 5L, descartadas = 0L)
  attr(x, "build") <- "GRCh38"
  f <- iei_filtrar_frecuencia(x)
  # BTK es XL (1e-4) y RAG1 es AR (0.01)
  expect_equal(f$af_gnomad, c(5e-5, 5e-3, NA))
  b <- iei_bitacora(f)
  expect_equal(b$paso[nrow(b)], "Frecuencia")
  expect_equal(b$descartadas[nrow(b)], 2L)
  expect_error(iei_filtrar_frecuencia(data.frame(a = 1)), "iei_anotar")
})

test_that("anotacion en vivo con Ensembl", {
  # Solo corre si activas: Sys.setenv(IEIPRIO_EN_VIVO = "true")
  skip_if_not(Sys.getenv("IEIPRIO_EN_VIVO") == "true", "test en vivo desactivado")
  # Variante comun conocida en GRCh38 (rs699, gen AGT)
  r <- .vep_consultar("1 230710048 . A G . . .", "GRCh38")
  t <- .vep_aplanar(r, c("1 230710048 . A G . . ." = "k"))
  expect_true("AGT" %in% t$gen_vep)
  expect_true(any(!is.na(t$mane)))
  expect_false(is.na(unique(t$af_gnomad)))
  expect_false(is.na(.vep_release("GRCh38")))
})
