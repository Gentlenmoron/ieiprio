skip_if_not_installed("ggplot2")

test_that("autoplot y plot funcionan para busquedas", {
  b <- suppressMessages(iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente"), n = 8))
  for (tipo in c("ranking", "mapa")) {
    g <- ggplot2::autoplot(b, tipo = tipo)
    expect_s3_class(g, "ggplot")
    expect_no_error(ggplot2::ggplot_build(g))
  }
  expect_s3_class(plot(b), "ggplot")
})

test_that("autoplot de una ficha de gen", {
  withr::local_options(ieiprio.sin_red = FALSE)
  local_mocked_bindings(.clinvar_peticion = clinvar_falso,
                        .proteina_canonica = function(gen, build) "ENSP1",
                        .proteina_consultar = proteina_falsa)
  g <- iei_gen("BTK", cache = FALSE)
  for (tipo in c("proteina", "clinvar")) {
    expect_no_error(ggplot2::ggplot_build(ggplot2::autoplot(g, tipo = tipo)))
  }
  expect_error(ggplot2::autoplot(iei_gen("BTK", clinvar = FALSE, proteina = FALSE)), "ClinVar")
})

test_that("autoplot de un resultado incluye el mapa de cohorte", {
  r <- resultado_mini()
  for (tipo in c("puntaje", "categorias", "genes", "mapa")) {
    expect_no_error(ggplot2::ggplot_build(ggplot2::autoplot(r, tipo = tipo)))
  }
})

test_that("autoplot de un cromatograma", {
  d <- withr::local_tempdir()
  escribir_ab1(file.path(d, "a.ab1"), referencia_prueba(120, 4))
  a <- iei_leer_ab1(file.path(d, "a.ab1"))
  g <- ggplot2::autoplot(a, desde = 10, hasta = 40)
  expect_no_error(ggplot2::ggplot_build(g))
  expect_match(g$labels$title, "10 a 40")
})

test_that("autoplot de una validacion de fenotipo", {
  skip_if_not_installed("ggplot2")
  v <- suppressMessages(iei_validar_fenotipo(n = 5, metodos = c("enfermedad", "exacta")))
  g <- ggplot2::autoplot(v, hasta = 10)
  expect_s3_class(g, "ggplot")
  expect_no_error(ggplot2::ggplot_build(g))
  expect_s3_class(plot(v), "ggplot")
})
