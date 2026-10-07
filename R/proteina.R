# Consulta a Ensembl el largo de una proteina y sus dominios Pfam.
# Devuelve list(largo, dominios) o NULL si falla.
.proteina_consultar <- function(proteina_id, build) {
  base <- httr2::request(.vep_servidor(build)) |>
    httr2::req_headers(Accept = "application/json") |>
    .agente() |>
    httr2::req_timeout(60) |>
    httr2::req_throttle(rate = 10 / 1) |>
    httr2::req_retry(max_tries = 3)
  sec <- base |>
    httr2::req_url_path_append("sequence", "id", proteina_id) |>
    httr2::req_url_query(type = "protein") |>
    httr2::req_perform() |>
    httr2::resp_body_json()
  dom <- base |>
    httr2::req_url_path_append("overlap", "translation", proteina_id) |>
    httr2::req_url_query(type = "Pfam") |>
    httr2::req_perform() |>
    httr2::resp_body_json()
  .proteina_procesar(sec, dom)
}

.proteina_procesar <- function(sec, dom) {
  largo <- nchar(.chr1(sec$seq))
  if (is.na(largo) || largo == 0) return(NULL)
  filas <- lapply(dom, function(d) {
    ini <- suppressWarnings(as.numeric(.chr1(d$start)))
    fin <- suppressWarnings(as.numeric(.chr1(d$end)))
    if (is.na(ini) || is.na(fin)) return(NULL)
    nombre <- .chr1(d$description)
    if (is.na(nombre) || nombre == "") nombre <- .chr1(d$id)
    data.frame(inicio = ini, fin = fin, nombre = nombre, stringsAsFactors = FALSE)
  })
  dominios <- do.call(rbind, filas)
  if (is.null(dominios)) dominios <- data.frame(inicio = numeric(0), fin = numeric(0),
                                                nombre = character(0))
  dominios <- unique(dominios[order(dominios$inicio), ])
  list(largo = largo, dominios = dominios)
}

# Obtiene la informacion de varias proteinas usando una cache en disco.
.proteinas_info <- function(ids, build, cache) {
  ids <- unique(ids[!is.na(ids)])
  if (length(ids) == 0) return(list())
  ruta <- if (isFALSE(cache)) NULL else file.path(cache, sprintf("proteinas_v1_%s.rds", build))
  guardado <- if (!is.null(ruta) && file.exists(ruta)) readRDS(ruta) else list()
  faltan <- setdiff(ids, names(guardado))
  if (length(faltan) > 0) {
    for (id in faltan) {
      guardado[[id]] <- tryCatch(.proteina_consultar(id, build), error = function(e) NULL)
    }
    guardado <- guardado[!vapply(guardado, is.null, logical(1))]
    if (!is.null(ruta)) {
      dir.create(cache, recursive = TRUE, showWarnings = FALSE)
      saveRDS(guardado, ruta)
    }
  }
  guardado[intersect(ids, names(guardado))]
}
