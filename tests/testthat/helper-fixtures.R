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

# Las pruebas en vivo corren solo con IEIPRIO_EN_VIVO=true y con conexion.
saltar_si_no_en_vivo <- function() {
  testthat::skip_if_not(Sys.getenv("IEIPRIO_EN_VIVO") == "true", "test en vivo desactivado")
  testthat::skip_if_offline("rest.ensembl.org")
}

# En las pruebas en vivo, una caida temporal de Ensembl no es un error del
# paquete: la prueba se salta con el motivo.
en_vivo <- function(code) {
  saltar <- function(e) testthat::skip(paste("servicio externo no disponible:", conditionMessage(e)))
  tryCatch(code, ieiprio_servicio_caido = saltar, httr2_failure = saltar)
}

# Resultado priorizado pequeno, con VEP simulado.
resultado_mini <- function() {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999",
                        .env = parent.frame())
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  a <- suppressMessages(iei_anotar(iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba()),
                                   cache = FALSE))
  iei_priorizar(a, fenotipo = "HP:0004313")
}


# Las pruebas sin conexion no deben intentar servicios externos: el modo sin
# red hace que las consultas en linea se omitan. Las pruebas en vivo lo apagan
# de forma explicita, por si la sesion de R lo dejo prendido de antes.
options(ieiprio.sin_red = !identical(Sys.getenv("IEIPRIO_EN_VIVO"), "true"),
        ieiprio.silenciar_privacidad = TRUE)
