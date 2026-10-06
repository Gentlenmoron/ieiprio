.panelapp_base <- "https://panelapp.genomicsengland.co.uk/api/v1"

.panelapp_url_panel <- function(panel_id, version = NULL) {
  url <- sprintf("%s/panels/%s/", .panelapp_base, panel_id)
  if (!is.null(version)) url <- paste0(url, "?version=", version)
  url
}

# Descarga un panel de PanelApp. Si `destino` no es NULL guarda el JSON crudo.
# Devuelve la respuesta como lista.
.panelapp_descargar <- function(panel_id = 398, version = NULL, destino = NULL) {
  req <- httr2::request(.panelapp_base) |>
    httr2::req_url_path_append("panels", paste0(panel_id, "/")) |>
    httr2::req_url_query(version = version) |>
    .agente() |>
    httr2::req_timeout(120) |>
    httr2::req_retry(max_tries = 4)

  texto <- httr2::resp_body_string(httr2::req_perform(req))
  if (!is.null(destino)) writeLines(texto, destino, useBytes = TRUE)
  jsonlite::fromJSON(texto, simplifyVector = FALSE)
}

# Convierte la respuesta de PanelApp en tablas ordenadas.
.panelapp_procesar <- function(x) {
  genes <- x$genes
  if (is.null(genes) || length(genes) == 0) {
    cli::cli_abort("La respuesta de PanelApp no trae genes.")
  }

  gd  <- lapply(genes, function(g) .o(g$gene_data, list()))
  gen <- vapply(seq_along(genes), function(i) {
    .chr1(.o(gd[[i]]$hgnc_symbol, .o(gd[[i]]$gene_symbol, genes[[i]]$entity_name)))
  }, character(1))
  hgnc_id <- vapply(gd, function(d) .chr1(d$hgnc_id), character(1))
  herencia_pa <- vapply(genes, function(g) .chr1(g$mode_of_inheritance), character(1))
  confianza <- vapply(genes, function(g) .chr1(g$confidence_level), character(1))
  fenotipos <- vapply(genes, function(g) {
    f <- unlist(g$phenotypes)
    f <- f[!is.na(f) & nzchar(trimws(f))]
    if (length(f) == 0) NA_character_ else paste(trimws(f), collapse = "; ")
  }, character(1))

  panel <- tibble::tibble(
    gen               = gen,
    hgnc_id           = hgnc_id,
    herencia          = .mapear_herencia(herencia_pa),
    herencia_panelapp = herencia_pa,
    evidencia         = .mapear_evidencia(confianza),
    fenotipos         = fenotipos
  )
  panel <- panel[!is.na(panel$gen), ]
  panel <- panel[order(panel$gen), ]

  alias_lista <- lapply(seq_along(genes), function(i) {
    a <- unique(c(unlist(gd[[i]]$alias), unlist(gd[[i]]$hgnc_previous_symbols)))
    a <- a[!is.na(a) & nzchar(a)]
    if (length(a) == 0 || is.na(gen[i])) return(NULL)
    data.frame(alias = a, gen = gen[i], stringsAsFactors = FALSE)
  })
  alias <- do.call(rbind, alias_lista)
  if (is.null(alias)) alias <- data.frame(alias = character(0), gen = character(0))
  alias <- unique(alias[toupper(alias$alias) != toupper(alias$gen), ])
  alias <- tibble::as_tibble(alias[order(alias$alias), ])

  list(panel = panel, alias = alias,
       version = .chr1(x$version), nombre = .chr1(x$name))
}

# Une el panel de PanelApp con la curacion IUIS.
.construir_panel <- function(pa, curado) {
  panel <- pa$panel
  cats <- tapply(curado$categoria_iuis, curado$gen,
                 function(v) paste(sort(unique(v)), collapse = ";"))
  panel$categorias_iuis <- unname(cats[panel$gen])
  gof <- unique(curado$gen[as.logical(curado$ganancia_funcion)])
  panel$ganancia_funcion <- panel$gen %in% gof
  panel
}

.construir_categorias <- function(curado, genes_panel) {
  cur <- curado[curado$gen %in% genes_panel, c("gen", "categoria_iuis")]
  cur <- unique(cur)
  tibble::tibble(
    gen              = cur$gen,
    categoria_iuis   = as.integer(cur$categoria_iuis),
    categoria_nombre = unname(.categorias_iuis[as.character(cur$categoria_iuis)])
  )
}

#' Descargar una version actualizada del panel desde PanelApp
#'
#' El panel incluido en el paquete esta fijado a una version de PanelApp
#' para que los resultados sean reproducibles (ver [iei_fuentes()]). Esta
#' funcion baja otra version, o la mas reciente, sin modificar los datos
#' del paquete, y la guarda en la cache del usuario.
#'
#' @param version Version del panel en PanelApp, por ejemplo `"8.78"`.
#'   `NULL` baja la mas reciente.
#' @param panel_id Identificador del panel en PanelApp. Por defecto 398,
#'   inmunodeficiencias primarias.
#' @param guardar Si `TRUE`, guarda el resultado como RDS en la cache.
#'
#' @return Un tibble con la misma estructura que [iuis_panel], con los
#'   atributos `version` y `panel_id`.
#' @export
#' @examples
#' \dontrun{
#' nuevo <- iei_actualizar_panel()
#' attr(nuevo, "version")
#' }
iei_actualizar_panel <- function(version = NULL, panel_id = 398, guardar = TRUE) {
  pa <- .panelapp_procesar(.panelapp_descargar(panel_id, version))
  curado <- ieiprio::iuis_categorias
  curado$ganancia_funcion <- curado$gen %in%
    ieiprio::iuis_panel$gen[ieiprio::iuis_panel$ganancia_funcion]
  panel <- .construir_panel(pa, curado)
  attr(panel, "version") <- pa$version
  attr(panel, "panel_id") <- panel_id

  if (guardar) {
    dir <- tools::R_user_dir("ieiprio", "cache")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    ruta <- file.path(dir, sprintf("panelapp_%s_v%s.rds", panel_id, pa$version))
    saveRDS(panel, ruta)
    cli::cli_inform("Panel {panel_id} version {pa$version} guardado en {.path {ruta}}.")
  }

  incluida <- ieiprio::fuentes_datos$version[ieiprio::fuentes_datos$fuente == "PanelApp"]
  if (!identical(pa$version, incluida)) {
    cli::cli_inform(c("i" = "El paquete incluye la version {incluida}; bajaste la {pa$version}."))
  }
  panel
}
