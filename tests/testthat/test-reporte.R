priorizado <- function() {
  x <- tibble::tibble(
    muestra = c("P01", "P01", "P01", "P02"), gen = c("BTK", "RAG1", "CYBB", "BTK"),
    chr = "X", pos = c(1L, 2L, 3L, 1L), ref = "C", alt = "T",
    id = c("rs999001", NA, NA, "rs999001"),
    genotipo = c("hemi", "het", "het", "het"), impacto = c("HIGH", "HIGH", "LOW", "HIGH"),
    consecuencia = c("stop_gained", "stop_gained", "synonymous_variant", "stop_gained"),
    af_gnomad = c(NA, NA, 0.001, NA),
    clinvar = c("pathogenic", "pathogenic", NA, "pathogenic"),
    hgvs_c = c("ENST1:c.1A>T", "ENST2:c.2A>T", NA, "ENST1:c.1A>T"),
    hgvs_p = c("ENSP1:p.Arg1Ter", "ENSP2:p.Lys2Ter", NA, "ENSP1:p.Arg1Ter")
  )
  attr(x, "ieiprio_bitacora") <- tibble::tibble(paso = c("Inicial", "Calidad"),
                                                filas = c(10L, 4L), descartadas = c(0L, 6L))
  attr(x, "build") <- "GRCh38"
  attr(x, "vep_release") <- "999"
  iei_priorizar(x, fenotipo = "HP:0004313")
}

test_that("genera un HTML con las secciones clave", {
  f <- withr::local_tempfile(fileext = ".html")
  expect_message(iei_reporte(priorizado(), f, paciente = "P01", responsable = "Lab <Inmuno>"),
                 "P01")
  h <- paste(readLines(f), collapse = "\n")
  expect_match(h, "Reporte de variantes &middot; P01")
  expect_match(h, "No es un informe diagn")
  expect_match(h, "<h2>Filtrado</h2>")
  expect_match(h, "Par&aacute;metros y trazabilidad")
  expect_match(h, "Control de calidad de la muestra")
  expect_match(h, "Coincidencia con el fenotipo")
  expect_match(h, "Desglose del puntaje")
  expect_match(h, "ACMG sugerido")
  expect_match(h, "release 999")
  expect_match(h, "p.Arg1Ter")
  expect_match(h, "gnomad.broadinstitute.org/variant/X-1-C-T\\?dataset=gnomad_r4")
  expect_match(h, "clinvar/\\?term=rs999001")
  expect_match(h, "HP:0004313")
  expect_match(h, "Portador en gen recesivo")
  expect_false(grepl("<b>CYBB</b> .", h, fixed = TRUE))
})

test_that("escapa el texto para no romper el HTML", {
  f <- withr::local_tempfile(fileext = ".html")
  suppressMessages(iei_reporte(priorizado(), f, paciente = "P01", responsable = "Lab <Inmuno>"))
  h <- paste(readLines(f), collapse = "\n")
  expect_match(h, "Lab &lt;Inmuno&gt;")
  expect_false(grepl("<Inmuno>", h, fixed = TRUE))
})

test_that("solo incluye la muestra pedida", {
  f <- withr::local_tempfile(fileext = ".html")
  suppressMessages(iei_reporte(priorizado(), f, paciente = "P02"))
  h <- paste(readLines(f), collapse = "\n")
  expect_false(grepl("RAG1", h, fixed = TRUE))
  expect_match(h, "portadora")
})

test_that("valida entradas", {
  f <- withr::local_tempfile(fileext = ".html")
  expect_error(iei_reporte(priorizado(), f, dominios = FALSE), "indica")
  expect_error(iei_reporte(priorizado(), f, paciente = "P99"), "P99")
  expect_error(iei_reporte(data.frame(a = 1), f), "iei_priorizar")
})

test_that("funciona sin variantes", {
  f <- withr::local_tempfile(fileext = ".html")
  vacio <- priorizado()[0, ]
  suppressMessages(iei_reporte(vacio, f, paciente = "P01"))
  expect_match(paste(readLines(f), collapse = ""), "No quedaron variantes")
})
