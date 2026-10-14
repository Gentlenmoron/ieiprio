respuesta <- function(status, tipo, cuerpo) {
  httr2::response(status_code = status, headers = list(`Content-Type` = tipo),
                  body = charToRaw(cuerpo))
}
req_prueba <- function() httr2::request("https://ejemplo.invalid/api")

test_that(".ejecutar_json devuelve el JSON cuando todo va bien", {
  withr::local_options(ieiprio.sin_red = FALSE)
  httr2::local_mocked_responses(function(req) respuesta(200, "application/json", '{"a": 1}'))
  expect_equal(.ejecutar_json(req_prueba()), list(a = 1L))
})

test_that("un error de servidor da un mensaje claro con clase propia", {
  withr::local_options(ieiprio.sin_red = FALSE)
  httr2::local_mocked_responses(function(req) respuesta(500, "text/html", "<html>error</html>"))
  e <- expect_error(.ejecutar_json(req_prueba(), "Ensembl"), class = "ieiprio_servicio_caido")
  expect_match(conditionMessage(e), "Ensembl no est. respondiendo \\(HTTP 500\\)")
})

test_that("una pagina de error con estado 200 tambien se detecta", {
  withr::local_options(ieiprio.sin_red = FALSE)
  httr2::local_mocked_responses(function(req) respuesta(200, "text/html", "<html>mantenimiento</html>"))
  expect_error(.ejecutar_json(req_prueba()), "pagina de error", class = "ieiprio_servicio_caido")
})

test_that("el Variant Recoder lee un 400 como variante no reconocida", {
  withr::local_options(ieiprio.sin_red = FALSE)
  httr2::local_mocked_responses(function(req) {
    respuesta(400, "application/json", '{"error": "No se pudo interpretar BTK:c.999X>Y"}')
  })
  r <- .recoder_consultar("BTK:c.999X>Y", "GRCh38")
  expect_match(r$error, "No se pudo interpretar")
  expect_equal(nrow(.recoder_aplanar(r)), 0)
})

test_that("los errores de servidor se consideran transitorios", {
  withr::local_options(ieiprio.sin_red = FALSE)
  expect_true(.transitorio(respuesta(503, "text/html", "")))
  expect_true(.transitorio(respuesta(500, "text/html", "")))
  expect_false(.transitorio(respuesta(400, "application/json", "{}")))
})

test_that("si falta la proteina, el reporte sigue sin diagrama", {
  withr::local_options(ieiprio.sin_red = FALSE)
  local_mocked_bindings(.proteina_consultar = function(id, build) {
    cli::cli_abort("caido", class = "ieiprio_servicio_caido")
  })
  expect_length(.proteinas_info("P1", "GRCh38", FALSE), 0)
})
