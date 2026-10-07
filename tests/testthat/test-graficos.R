test_that("el desglose del puntaje suma el total", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999")
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  a <- suppressMessages(iei_anotar(iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba()),
                                   cache = FALSE))
  p <- iei_priorizar(a)
  suma <- p$p_clinvar + p$p_impacto + p$p_frecuencia + p$p_herencia + p$p_fenotipo
  expect_equal(suma, p$puntaje)
})

test_that("los criterios ACMG sugeridos siguen las reglas", {
  x <- tibble::tibble(
    consecuencia = c("stop_gained", "missense_variant", "frameshift_variant", "missense_variant"),
    af_gnomad = c(NA, 0.001, 5e-5, 0.08))
  r <- .acmg_sugeridos(x, gof = c(FALSE, FALSE, TRUE, FALSE),
                       compuesto = c(FALSE, TRUE, FALSE, FALSE), p_fenotipo = c(3, 0, 0, 0))
  expect_equal(r, c("PVS1, PM2_Supporting, PP4", "PM3", "PM2_Supporting", "BA1"))
})

test_that("los graficos producen SVG valido", {
  h <- .histograma(c(5, 15, 15, 120, NA), c(seq(0, 100, 10), Inf))
  expect_equal(sum(h$n), 4)
  expect_match(.svg_histograma(h, "DP", "lecturas"), "^<figure.*<svg.*</svg></figure>$")
  expect_match(.svg_histograma(.histograma(numeric(0), 0:2), "DP", "x"), "Sin datos")

  b <- tibble::tibble(paso = c("Inicial", "Panel"), filas = c(100L, 5L), descartadas = c(0L, 95L))
  e <- .svg_embudo(b)
  expect_match(e, "Inicial")
  expect_match(e, "\\(-95\\)")

  d <- tibble::tibble(gen = c("BTK", "RAG1"), hgvs_p = c("E:p.Arg1Ter", NA), puntaje = c(10, 5),
                      p_clinvar = c(4, 4), p_impacto = c(3, 3), p_frecuencia = c(1, 0),
                      p_herencia = c(2, -2), p_fenotipo = c(0, 0))
  s <- .svg_puntaje(d)
  expect_match(s, "<svg")
  expect_equal(lengths(regmatches(s, gregexpr("<rect", s))), 5 + 4 + 3)
})

test_that("diagrama de proteina con dominios y marcas", {
  j <- jsonlite::read_json(fixture("proteina_sintetica.json"))
  info <- .proteina_procesar(j$secuencia, j$dominios)
  expect_equal(info$largo, 659)
  expect_equal(nrow(info$dominios), 3)
  marcas <- data.frame(pos = c(334, 500), etiqueta = c("Arg334Trp", "Lys500Ter"),
                       categoria = c("media", "alta"))
  s <- .svg_proteina(info, marcas)
  expect_match(s, "659 aa")
  expect_match(s, "SH2 domain")
  expect_match(s, "Lys500Ter")
  expect_match(.svg_proteina(NULL, marcas), "Sin informacion")
})

test_that("la informacion de proteinas usa cache", {
  dir <- withr::local_tempdir()
  n <- 0
  local_mocked_bindings(.proteina_consultar = function(id, build) {
    n <<- n + 1
    j <- jsonlite::read_json(fixture("proteina_sintetica.json"))
    .proteina_procesar(j$secuencia, j$dominios)
  })
  r1 <- .proteinas_info(c("P1", "P1", NA), "GRCh38", dir)
  r2 <- .proteinas_info("P1", "GRCh38", dir)
  expect_equal(n, 1)
  expect_equal(r2$P1$largo, 659)
  local_mocked_bindings(.proteina_consultar = function(id, build) stop("sin red"))
  expect_length(.proteinas_info("P2", "GRCh38", FALSE), 0)
})

test_that("matriz de fenotipo", {
  hpo <- tibble::tibble(gen = c("A", "A", "B"), hpo_id = c("HP:1", "HP:2", "HP:2"), hpo_nombre = "x")
  t <- .tabla_fenotipo(tibble::tibble(gen = c("A", "B")), c("HP:1", "HP:2"), hpo)
  expect_match(t, "2 de 2")
  expect_match(t, "1 de 2")
  expect_equal(.tabla_fenotipo(tibble::tibble(gen = "A"), NULL, hpo), "")
})

test_that("reporte completo con dominios y control de calidad", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999",
                        .proteina_consultar = function(id, build) {
                          j <- jsonlite::read_json(fixture("proteina_sintetica.json"))
                          .proteina_procesar(j$secuencia, j$dominios)
                        })
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  a <- suppressMessages(iei_anotar(iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba()),
                                   cache = FALSE))
  p <- iei_priorizar(a, fenotipo = "HP:0004313")
  f <- withr::local_tempfile(fileext = ".html")
  suppressMessages(iei_reporte(p, f, paciente = "P01", cache = FALSE))
  h <- paste(readLines(f), collapse = "\n")
  expect_match(h, "Verificaci&oacute;n de sexo")
  expect_match(h, "mini.vcf")
  expect_match(h, "Profundidad \\(DP\\)")
  expect_match(h, "659 aa")
  expect_match(h, "Arg334Trp")
  expect_match(h, "Lys500Ter")
  expect_match(h, "Embudo del filtrado")
  expect_match(h, "DP &ge; 10")
})

test_that("proteina en vivo con Ensembl", {
  skip_if_not(Sys.getenv("IEIPRIO_EN_VIVO") == "true", "test en vivo desactivado")
  # BTK, proteina canonica en Ensembl
  info <- .proteina_consultar("ENSP00000308176", "GRCh38")
  expect_gt(info$largo, 600)
  expect_gt(nrow(info$dominios), 0)
})
