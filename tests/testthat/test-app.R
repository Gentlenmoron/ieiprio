test_that("la aplicacion se arma y responde a la busqueda por fenotipo", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("ggplot2")
  # La aplicacion cambia opciones globales; se restauran al terminar la prueba
  # para no apagar la red en las pruebas en vivo que vienen despues.
  withr::local_options(ieiprio.sin_red = getOption("ieiprio.sin_red"),
                       ieiprio.silenciar_privacidad = getOption("ieiprio.silenciar_privacidad"))
  app <- shiny::shinyAppFile(system.file("app", "app.R", package = "ieiprio"))
  expect_s3_class(app, "shiny.appobj")
  shiny::testServer(app, {
    session$setInputs(sin_red = TRUE, fen = "agammaglobulinemia\nneumonia recurrente", aus = "",
                      sim = "enfermedad", n = 10, buscar = 1)
    b <- suppressMessages(busqueda())
    expect_s3_class(b, "iei_busqueda")
    expect_true("BTK" %in% b$gen)
    session$setInputs(val_n = 20, val_imp = 0.3, val_ruido = 1, validar = 1)
    expect_s3_class(suppressMessages(validacion()), "iei_validacion")
  })
})

test_that("iei_app avisa si falta shiny", {
  local_mocked_bindings(check_installed = function(...) stop("falta shiny"), .package = "rlang")
  expect_error(iei_app(), "falta shiny")
})
