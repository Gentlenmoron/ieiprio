fixture <- function(nombre) testthat::test_path("fixtures", nombre)
pa_mini <- function() .panelapp_procesar(jsonlite::read_json(fixture("panelapp_mini.json")))
