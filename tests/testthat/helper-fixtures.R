fixture <- function(nombre) testthat::test_path("fixtures", nombre)
pa_mini <- function() .panelapp_procesar(jsonlite::read_json(fixture("panelapp_mini.json")))

coords_prueba <- function() {
  tibble::tibble(
    gen = c("GENA", "BTK", "BTK"),
    build = c("GRCh38", "GRCh38", "GRCh37"),
    chr = c("1", "X", "X"),
    inicio = c(900L, 101349000L, 100604435L),
    fin = c(3100L, 101391000L, 100641212L),
    ensembl_id = NA_character_
  )
}

vep_sintetico <- function() jsonlite::read_json(fixture("vep_sintetico.json"))

# Simula la API de VEP: responde solo las entradas que estan en el fixture.
vep_falso <- function(entrada, build) {
  r <- vep_sintetico()
  Filter(function(o) o$input %in% entrada, r)
}
