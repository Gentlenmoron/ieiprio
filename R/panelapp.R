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

  coords <- do.call(rbind, lapply(seq_along(genes), function(i) {
    .coordenadas_gen(gd[[i]], gen[i])
  }))
  if (is.null(coords)) {
    coords <- data.frame(gen = character(0), build = character(0), chr = character(0),
                         inicio = integer(0), fin = integer(0), ensembl_id = character(0))
  }
  coords <- tibble::as_tibble(coords[order(coords$build, coords$gen), ])

  list(panel = panel, alias = alias, coordenadas = coords,
       version = .chr1(x$version), nombre = .chr1(x$name))
}

# Extrae las coordenadas de un gen desde gene_data$ensembl_genes, que tiene la
# forma list(GRch37 = list(`82` = list(location = "X:1-2", ensembl_id = ...)),
#            GRch38 = list(`90` = list(...))).
.coordenadas_gen <- function(gd, gen) {
  eg <- gd$ensembl_genes
  if (is.null(eg) || length(eg) == 0 || is.na(gen)) return(NULL)
  filas <- list()
  for (nombre in names(eg)) {
    build <- if (grepl("37", nombre)) "GRCh37" else if (grepl("38", nombre)) "GRCh38" else NA
    if (is.na(build)) next
    for (version in eg[[nombre]]) {
      loc <- .chr1(version$location)
      partes <- regmatches(loc, regexec("^([^:]+):([0-9]+)-([0-9]+)$", loc))[[1]]
      if (length(partes) != 4) next
      filas[[length(filas) + 1]] <- data.frame(
        gen = gen, build = build, chr = .normalizar_chr(partes[2]),
        inicio = as.integer(partes[3]), fin = as.integer(partes[4]),
        ensembl_id = .chr1(version$ensembl_id), stringsAsFactors = FALSE
      )
      break
    }
  }
  if (length(filas) == 0) NULL else do.call(rbind, filas)
}

# Une el panel de PanelApp con la curacion IUIS.
# Nombres de las diez tablas IUIS tal como aparecen en los fenotipos de
# PanelApp, y palabras clave de fenotipo que indican cada categoria cuando el
# nombre de la tabla no esta.
.iuis_tablas <- c(
  "1" = "cellular and humoral immunity",
  "2" = "combined immunodeficiencies with associated or syndromic",
  "3" = "predominantly antibody deficienc",
  "4" = "diseases of immune dysregulation",
  "5" = "defects of phagocyte",
  "6" = "defects in intrinsic and innate immunity",
  "7" = "autoinflammatory disorders",
  "8" = "complement deficienc",
  "9" = "bone marrow failure",
  "10" = "phenocopies"
)
.iuis_palabras <- c(
  "1" = "severe combined immunodeficiency|\\bSCID\\b|combined immunodeficiency|combined immune deficiency|omenn",
  "2" = "hyper[- ]?ige|DNA repair|ataxia[- ]telangiectasia|nijmegen|wiskott|immuno-?osseous|cartilage[- ]hair|digeorge|thymic (aplasia|hypoplasia)|ectodermal dysplasia.*immunodeficiency",
  "3" = "agammaglobulin|hypogammaglobulin|common variable|\\bCVID\\b|IgA deficiency|antibody deficiency|class[- ]switch",
  "4" = "hemophagocytic|\\bHLH\\b|lymphoproliferative|\\bALPS\\b|\\bIPEX\\b|regulatory T|immune dysregulation",
  "5" = "neutropeni|granulomatous disease|leukocyte adhesion",
  "6" = "mycobacter|\\bMSMD\\b|candidiasis|herpes simplex encephalitis|epidermodysplasia|susceptibility to (viral|invasive|severe)|WHIM",
  "7" = "autoinflamm|interferonopath|aicardi|periodic fever|muckle|cold.{0,20}inflammatory|\\bCINCA\\b|mediterranean fever|mevalonate|hyper[- ]?IgD",
  "8" = "complement|properdin|angioedema",
  "9" = "bone marrow failure|dyskeratosis|aplastic anemia|fanconi"
)

# Categorias IUIS de cada gen en formato largo, con la fuente de cada una.
.categorias_iuis_derivar <- function(genes, fenotipos, curado = NULL) {
  f <- ifelse(is.na(fenotipos), "", fenotipos)
  buscar <- function(patrones, fuente) {
    do.call(rbind, lapply(names(patrones), function(k) {
      g <- genes[grepl(patrones[[k]], f, ignore.case = TRUE, perl = TRUE)]
      if (length(g) == 0) return(NULL)
      data.frame(gen = g, categoria_iuis = as.integer(k), fuente = fuente)
    }))
  }
  tabla <- buscar(.iuis_tablas, "PanelApp, tabla IUIS")
  palabras <- buscar(.iuis_palabras, "PanelApp, fenotipo")
  # Las palabras clave solo completan genes que no traen la tabla IUIS.
  palabras <- palabras[!palabras$gen %in% tabla$gen, ]
  cur <- NULL
  if (!is.null(curado) && nrow(curado) > 0) {
    cur <- data.frame(gen = curado$gen, categoria_iuis = as.integer(curado$categoria_iuis),
                      fuente = "curado")
    cur <- cur[!is.na(cur$categoria_iuis) & cur$gen %in% genes, ]
  }
  todo <- rbind(tabla, palabras, cur)
  todo <- todo[!duplicated(todo[c("gen", "categoria_iuis")]), ]
  todo <- todo[order(todo$gen, todo$categoria_iuis), ]
  tibble::tibble(
    gen = todo$gen,
    categoria_iuis = todo$categoria_iuis,
    categoria_nombre = unname(.categorias_iuis[as.character(todo$categoria_iuis)]),
    fuente = todo$fuente
  )
}

# `mecanismo` tiene columnas gen y mecanismo ("GOF" o "LOF y GOF").
.construir_panel <- function(pa, curado, mecanismo = NULL) {
  panel <- pa$panel
  cats <- .categorias_iuis_derivar(panel$gen, panel$fenotipos, curado)
  agrupadas <- tapply(cats$categoria_iuis, cats$gen, function(v) paste(sort(unique(v)), collapse = ";"))
  panel$categorias_iuis <- as.character(unname(agrupadas[panel$gen]))
  if (is.null(mecanismo) && !is.null(curado$ganancia_funcion)) {
    gof <- unique(curado$gen[as.logical(curado$ganancia_funcion)])
    mecanismo <- data.frame(gen = gof, mecanismo = rep("GOF", length(gof)))
  }
  panel$mecanismo <- if (is.null(mecanismo)) NA_character_ else mecanismo$mecanismo[match(panel$gen, mecanismo$gen)]
  panel$ganancia_funcion <- !is.na(panel$mecanismo)
  panel
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
  curado <- ieiprio::iuis_categorias[ieiprio::iuis_categorias$fuente == "curado", ]
  mec <- ieiprio::iuis_panel[!is.na(ieiprio::iuis_panel$mecanismo), c("gen", "mecanismo")]
  panel <- .construir_panel(pa, curado, mec)
  attr(panel, "version") <- pa$version
  attr(panel, "panel_id") <- panel_id

  if (guardar) {
    dir <- tools::R_user_dir("ieiprio", "cache")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    ruta <- file.path(dir, sprintf("panelapp_%s_v%s.rds", panel_id, pa$version))
    saveRDS(panel, ruta)
    cli::cli_inform("Panel {panel_id} versi\u00f3n {pa$version} guardado en {.path {ruta}}.")
  }

  incluida <- ieiprio::fuentes_datos$version[ieiprio::fuentes_datos$fuente == "PanelApp"]
  if (!identical(pa$version, incluida)) {
    cli::cli_inform(c("i" = "El paquete incluye la versi\u00f3n {incluida}; bajaste la {pa$version}."))
  }
  panel
}
