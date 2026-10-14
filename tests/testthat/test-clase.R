test_that("iei_priorizar devuelve un iei_resultado que sigue siendo tibble", {
  r <- resultado_mini()
  expect_s3_class(r, "iei_resultado")
  expect_s3_class(r, "tbl_df")
  f <- r[r$muestra == "P01", ]
  expect_s3_class(f, "iei_resultado")
  expect_equal(attr(f, "vep_release"), "999")
  skip_if_not_installed("dplyr")
  d <- dplyr::filter(r, categoria != "baja")
  expect_s3_class(d, "iei_resultado")
})

test_that("print muestra un resumen y no falla con columnas sueltas", {
  r <- resultado_mini()
  tabla <- utils::capture.output(out <- cli::cli_fmt(print(r)))
  out <- c(out, tabla)
  expect_true(any(grepl("Resultado de ieiprio", out)))
  expect_true(any(grepl("Ensembl VEP 999", out)))
  expect_output(print(r[, c("gen", "pos")]), "gen")
  expect_true(any(grepl("p.Arg334Trp", tabla, fixed = TRUE)))
})

test_that("summary resume por categoria, gen, muestra y ACMG", {
  s <- summary(resultado_mini())
  expect_s3_class(s, "summary.iei_resultado")
  expect_equal(sum(s$por_categoria), 3)
  expect_equal(s$por_gen$gen, "BTK")
  expect_setequal(s$por_muestra$muestra, c("P01", "P02"))
  expect_true("PM2_Supporting" %in% s$acmg$criterio)
  expect_false(is.null(s$bitacora))
  tablas <- utils::capture.output(lineas <- cli::cli_fmt(print(s)))
  expect_true(length(c(lineas, tablas)) > 5)
})

test_that("plot devuelve graficos de ggplot2", {
  skip_if_not_installed("ggplot2")
  r <- resultado_mini()
  for (tipo in c("puntaje", "categorias", "genes", "predictores")) {
    g <- plot(r, tipo = tipo)
    expect_s3_class(g, "ggplot")
    expect_no_error(ggplot2::ggplot_build(g))
  }
  sin_pred <- r; sin_pred$revel <- NA_real_; sin_pred$alphamissense <- NA_real_
  expect_error(plot(sin_pred, tipo = "predictores"), "REVEL")
  expect_error(plot(r[0, ]), "no tiene variantes")
})

test_that("plot interactivo usa plotly", {
  skip_if_not_installed("plotly")
  expect_s3_class(plot(resultado_mini(), interactivo = TRUE), "plotly")
})

test_that("iei_tabla e iei_metadatos", {
  r <- resultado_mini()
  t <- iei_tabla(r)
  expect_false(inherits(t, "iei_resultado"))
  expect_true(all(c("gen", "puntaje", "razones", "revel") %in% names(t)))
  expect_gt(ncol(iei_tabla(r, "todas")), ncol(t))
  m <- iei_metadatos(r)
  expect_true(all(c("build", "vep_release", "bitacora", "qc", "parametros", "vcf",
                    "fenotipo", "pesos") %in% names(m)))
  expect_equal(m$fenotipo, "HP:0004313")
})
