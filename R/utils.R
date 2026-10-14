.o <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.chr1 <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- x[[1]]
  if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
}

.agente <- function(req) {
  httr2::req_user_agent(req, "ieiprio (https://github.com/Gentlenmoron/ieiprio)")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Atributos que el paquete arrastra de una etapa a otra.
.attrs_ieiprio <- function(x) {
  a <- attributes(x)
  a[names(a) %in% c("build", "vep_release", "archivo") | startsWith(names(a), "ieiprio_")]
}

.poner_attrs <- function(x, a) {
  for (n in names(a)) attr(x, n) <- a[[n]]
  x
}

# Errores de servidor que suelen resolverse solos: se reintentan.
.transitorio <- function(resp) httr2::resp_status(resp) %in% c(429, 500, 502, 503, 504)

# Ejecuta una consulta a un servicio externo y devuelve el JSON. Si el
# servicio esta caido o responde algo que no es JSON (una pagina de error),
# lanza un error claro con clase `ieiprio_servicio_caido`.
.ejecutar_json <- function(req, servicio = "Ensembl") {
  .comprobar_red(servicio)
  resp <- tryCatch(httr2::req_perform(req), httr2_error = function(e) e)
  if (inherits(resp, "error")) {
    estado <- tryCatch(httr2::resp_status(resp$resp), error = function(e) NA)
    cli::cli_abort(c(
      "{servicio} no est\u00e1 respondiendo{if (!is.na(estado)) paste0(' (HTTP ', estado, ')')}.",
      "i" = "Suele ser una ca\u00edda temporal del servidor. Intenta de nuevo en unos minutos."
    ), class = "ieiprio_servicio_caido", parent = resp)
  }
  if (!grepl("json", httr2::resp_content_type(resp) %||% "", fixed = TRUE)) {
    cli::cli_abort(c(
      "{servicio} respondio con una pagina de error en lugar de datos (HTTP {httr2::resp_status(resp)}).",
      "i" = "Suele ser mantenimiento o una ca\u00edda temporal. Intenta de nuevo en unos minutos."
    ), class = "ieiprio_servicio_caido")
  }
  httr2::resp_body_json(resp, simplifyVector = FALSE)
}

# Con options(ieiprio.sin_red = TRUE) ninguna funcion consulta servicios
# externos: util cuando los datos no pueden salir de la institucion.
.comprobar_red <- function(servicio) {
  if (isTRUE(getOption("ieiprio.sin_red", FALSE))) {
    cli::cli_abort(c("Consulta a {servicio} bloqueada porque {.code options(ieiprio.sin_red = TRUE)} esta activo.",
                     "i" = "Para anotar sin conexi\u00f3n usa un VEP local con {.fn iei_vep_entrada} y {.code iei_anotar(vep_json = )}."),
                   class = "ieiprio_sin_red")
  }
  invisible(TRUE)
}

.avisos_mostrados <- new.env(parent = emptyenv())

.sin_red <- function() isTRUE(getOption("ieiprio.sin_red", FALSE))

# Una vez por sesion recuerda que las variantes salen hacia un servicio externo.
.aviso_privacidad <- function(servicio) {
  if (.sin_red()) return(invisible())
  if (isTRUE(getOption("ieiprio.silenciar_privacidad", FALSE)) || isTRUE(.avisos_mostrados[[servicio]])) {
    return(invisible())
  }
  assign(servicio, TRUE, envir = .avisos_mostrados)
  cli::cli_inform(c("i" = "Las posiciones de las variantes se enviar\u00e1n a {servicio}. Si los datos no pueden salir de tu instituci\u00f3n, usa un VEP local (ver {.fn iei_vep_entrada}).",
                    " " = "Este aviso aparece una vez por sesi\u00f3n; se apaga con {.code options(ieiprio.silenciar_privacidad = TRUE)}."))
  invisible()
}
