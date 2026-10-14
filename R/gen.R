# Ficha de un gen: lo que el paquete sabe de el y lo que se conoce en ClinVar.

# Transcrito canonico de Ensembl y su proteina.
.proteina_canonica <- function(gen, build) {
  r <- httr2::request(.vep_servidor(build)) |>
    httr2::req_url_path_append("lookup", "symbol", "homo_sapiens", gen) |>
    httr2::req_url_query(expand = 1) |>
    httr2::req_headers(Accept = "application/json") |>
    .agente() |>
    httr2::req_timeout(60) |>
    httr2::req_retry(max_tries = 3, is_transient = .transitorio) |>
    .ejecutar_json("Ensembl")
  tx <- r$Transcript %||% list()
  canon <- Filter(function(t) isTRUE(as.logical(t$is_canonical)), tx)
  if (length(canon) == 0) return(NA_character_)
  .chr1(canon[[1]]$Translation$id)
}

# Guarda en la cache por dia el resultado de una consulta en linea. Si falla,
# avisa y devuelve NULL.
.en_cache <- function(clave, simbolo, cache, f) {
  ruta <- if (isFALSE(cache)) NULL else file.path(cache, sprintf("%s_%s_%s.rds", clave, simbolo, Sys.Date()))
  if (!is.null(ruta) && file.exists(ruta)) return(readRDS(ruta))
  if (.sin_red()) {
    cli::cli_inform(c("i" = "Consulta en l\u00ednea omitida para {simbolo} (modo sin red)."))
    return(NULL)
  }
  r <- tryCatch(f(), error = function(e) {
    cli::cli_warn(c("No se pudo completar la consulta en l\u00ednea para {simbolo}.", "i" = conditionMessage(e)))
    NULL
  })
  if (!is.null(r) && !is.null(ruta)) {
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    saveRDS(r, ruta)
  }
  r
}

# Estrellas de ClinVar segun el estado de revision.
.estrellas_clinvar <- function(revision) {
  r <- tolower(ifelse(is.na(revision), "", revision))
  ifelse(grepl("practice guideline", r), 4L,
  ifelse(grepl("expert panel", r), 3L,
  ifelse(grepl("multiple submitters, no conflicts", r), 2L,
  ifelse(grepl("criteria provided", r), 1L, 0L))))
}

# Posicion en la proteina a partir del titulo de ClinVar, por ejemplo
# "NM_000061.3(BTK):c.1630A>G (p.Arg544Gly)", o del cambio en una letra.
.pos_proteina_clinvar <- function(titulo, cambio) {
  p <- suppressWarnings(as.integer(sub(".*\\(p\\.[A-Z][a-z]{2}([0-9]+).*", "\\1", titulo)))
  p[!grepl("\\(p\\.[A-Z][a-z]{2}[0-9]+", titulo %||% "")] <- NA
  sin <- is.na(p) & !is.na(cambio)
  p[sin] <- suppressWarnings(as.integer(sub("^[A-Z]([0-9]+).*", "\\1", cambio[sin])))
  p
}

.grupo_consecuencia <- function(cons, titulo) {
  x <- tolower(paste(cons, titulo))
  ifelse(grepl("frameshift|fs\\b|fs\\*|del|dup|ins", x) & !grepl("missense", x), "frameshift o indel",
  ifelse(grepl("nonsense|stop gained|ter\\)|\\*\\)", x), "sin sentido",
  ifelse(grepl("splice|[+-][0-9]+[acgt]>", x), "splicing",
  ifelse(grepl("missense", x), "missense",
  ifelse(grepl("synonymous", x), "sinonima", "otra")))))
}

# Variantes de ClinVar del gen: conteo por clasificacion y detalle de las
# patogenicas y probablemente patogenicas.
.clinvar_gen <- function(gen, max_variantes = 500) {
  clases <- c("Patogenica" = "clinsig_pathogenic", "Probablemente patogenica" = "clinsig_likely_pathogenic",
              "Significado incierto" = "clinsig_vus", "Probablemente benigna" = "clinsig_likely_benign",
              "Benigna" = "clinsig_benign")
  conteo <- vapply(clases, function(cl) {
    r <- .clinvar_peticion("esearch.fcgi", term = sprintf("%s[gene] AND %s[Properties]", gen, cl), retmax = 0)
    as.integer(r$esearchresult$count %||% 0)
  }, integer(1))
  ids <- .clinvar_buscar(sprintf("%s[gene] AND (clinsig_pathogenic[Properties] OR clinsig_likely_pathogenic[Properties])", gen),
                         retmax = max_variantes)
  res <- .clinvar_resumen(ids)
  filas <- lapply(seq_along(res), function(i) {
    r <- res[[i]]
    genes <- vapply(r$genes %||% list(), function(g) .chr1(g$symbol), character(1))
    if (length(genes) > 0 && !gen %in% genes) return(NULL)
    gc <- r$germline_classification %||% r$clinical_significance %||% list()
    titulo <- .chr1(r$title)
    cambio <- .chr1(strsplit(.chr1(r$protein_change) %||% "", ",\\s*")[[1]])
    cons <- paste(unlist(r$molecular_consequence_list %||% list()), collapse = ", ")
    id <- .chr1(r$uid)
    if (is.na(id)) id <- names(res)[i]
    data.frame(clinvar_id = as.integer(id),
               variante = titulo,
               clasificacion = .chr1(gc$description),
               revision = .chr1(gc$review_status),
               consecuencia = if (nzchar(cons)) cons else NA_character_,
               cambio_proteina = cambio, stringsAsFactors = FALSE)
  })
  v <- do.call(rbind, filas)
  if (is.null(v)) {
    v <- data.frame(clinvar_id = integer(0), variante = character(0), clasificacion = character(0),
                    revision = character(0), consecuencia = character(0), cambio_proteina = character(0))
  }
  v$estrellas <- .estrellas_clinvar(v$revision)
  v$pos_proteina <- .pos_proteina_clinvar(v$variante, v$cambio_proteina)
  v$tipo <- .grupo_consecuencia(v$consecuencia, v$variante)
  list(fecha = format(Sys.Date()),
       conteo = tibble::tibble(clasificacion = names(conteo), variantes = unname(conteo)),
       variantes = tibble::as_tibble(v[order(-v$estrellas, v$pos_proteina), ]),
       consultadas = length(ids))
}

#' Ficha de un gen
#'
#' Reune todo lo que se sabe de un gen del panel: herencia, evidencia en
#' PanelApp, categorias IUIS (con la fuente de cada una), mecanismo de
#' enfermedad, alias, ubicacion, fenotipos de PanelApp y terminos HPO. Si hay
#' conexion agrega las variantes conocidas en ClinVar y la proteina con sus
#' dominios Pfam, para ver donde se concentran las variantes patogenicas.
#'
#' Las consultas en linea se guardan en la cache por dia, asi que repetir la
#' misma ficha no vuelve a consultar. Sin conexion la ficha sale igual con la
#' informacion incluida en el paquete.
#'
#' @param gen Simbolo del gen o un alias (`"BTK"`, `"TACI"`).
#' @param clinvar Si `TRUE`, consulta las variantes de ClinVar.
#' @param proteina Si `TRUE`, consulta la proteina y sus dominios en Ensembl.
#' @param build `"GRCh38"` o `"GRCh37"`.
#' @param cache Carpeta de la cache, o `FALSE` para no usarla.
#' @return Una lista de clase `iei_gen`. Se imprime como ficha y se grafica
#'   con [plot()] o `autoplot()`; ver [autoplot.iei_gen()].
#' @export
#' @examples
#' g <- iei_gen("BTK", clinvar = FALSE, proteina = FALSE)
#' g
#' \dontrun{
#' g <- iei_gen("BTK")
#' g$clinvar$variantes
#' plot(g)
#' }
iei_gen <- function(gen, clinvar = TRUE, proteina = TRUE, build = c("GRCh38", "GRCh37"),
                    cache = tools::R_user_dir("ieiprio", "cache")) {
  build <- match.arg(build)
  if (length(gen) != 1) cli::cli_abort("{.fn iei_gen} recibe un solo gen. Para varios usa {.fn iei_panel} o {.fn iei_buscar_gen}.")
  oficial <- .normalizar_simbolo(gen, ieiprio::iuis_panel$gen, ieiprio::iuis_alias)
  en_panel <- !is.na(oficial)
  simbolo <- if (en_panel) oficial else toupper(gen)
  if (!en_panel) cli::cli_warn("{.val {gen}} no est\u00e1 en el panel de errores innatos de la inmunidad; la ficha solo tendr\u00e1 lo que se consulte en l\u00ednea.")

  p <- ieiprio::iuis_panel[ieiprio::iuis_panel$gen == simbolo, ]
  cats <- ieiprio::iuis_categorias[ieiprio::iuis_categorias$gen == simbolo, ]
  hpo <- ieiprio::hpo_genes[ieiprio::hpo_genes$gen == simbolo, ]
  if (nrow(hpo) > 0) {
    o <- .ontologia()
    hpo$ic <- unname(o$ic[o$idx[hpo$hpo_id]])
    hpo <- hpo[order(-hpo$ic, hpo$hpo_nombre), ]
  }
  coord <- ieiprio::genes_coordenadas[ieiprio::genes_coordenadas$gen == simbolo, ]

  en_cache <- function(clave, f) .en_cache(clave, simbolo, cache, f)
  cv <- if (clinvar) en_cache("clinvar_gen", function() .clinvar_gen(simbolo)) else NULL
  prot <- if (proteina) en_cache(paste0("proteina_gen_", build), function() {
    id <- .proteina_canonica(simbolo, build)
    if (is.na(id)) return(NULL)
    c(list(id = id), .proteina_consultar(id, build))
  }) else NULL

  structure(list(
    gen = simbolo,
    en_panel = en_panel,
    herencia = if (nrow(p)) p$herencia else NA_character_,
    herencia_panelapp = if (nrow(p)) p$herencia_panelapp else NA_character_,
    evidencia = if (nrow(p)) p$evidencia else NA_character_,
    mecanismo = if (nrow(p) && "mecanismo" %in% names(p)) p$mecanismo else NA_character_,
    fenotipos = if (nrow(p)) p$fenotipos else NA_character_,
    categorias = tibble::as_tibble(cats),
    alias = ieiprio::iuis_alias$alias[ieiprio::iuis_alias$gen == simbolo],
    coordenadas = tibble::as_tibble(coord),
    hpo = tibble::as_tibble(hpo),
    clinvar = cv,
    proteina = prot,
    build = build
  ), class = "iei_gen")
}

#' @export
print.iei_gen <- function(x, ...) {
  b <- cli::symbol$bullet
  cli::cli_h2(x$gen)
  if (!x$en_panel) cli::cli_text("{.emph Fuera del panel de errores innatos de la inmunidad.}")
  if (x$en_panel) {
    cli::cli_text("{b} Herencia: {x$herencia %|NA|% 'sin dato'} {cli::symbol$line} evidencia en PanelApp: {x$evidencia}")
    mec <- switch(x$mecanismo %|NA|% "", "GOF" = "ganancia de funci\u00f3n", "LOF y GOF" = "p\u00e9rdida y ganancia de funci\u00f3n seg\u00fan la enfermedad", "p\u00e9rdida de funci\u00f3n o no especificado")
    cli::cli_text("{b} Mecanismo: {mec}")
    if (nrow(x$categorias) > 0) {
      cats <- unique(sprintf("%d %s (%s)", x$categorias$categoria_iuis, x$categorias$categoria_nombre, x$categorias$fuente))
      cli::cli_text("{b} Categor\u00eda IUIS: {paste(cats, collapse = '; ')}")
    } else {
      cli::cli_text("{b} Categor\u00eda IUIS: sin asignar")
    }
    if (length(x$alias) > 0) cli::cli_text("{b} Alias: {paste(utils::head(x$alias, 8), collapse = ', ')}")
  }
  co <- x$coordenadas[x$coordenadas$build == x$build, ]
  if (nrow(co) > 0) cli::cli_text("{b} Ubicaci\u00f3n: {co$chr[1]}:{format(co$inicio[1], big.mark = ',')}-{format(co$fin[1], big.mark = ',')} ({x$build})")
  if (!is.na(x$fenotipos)) {
    fen <- .resumir_fenotipos(x$fenotipos)
    cli::cli_text("{b} Fenotipos en PanelApp: {fen}")
  }
  if (nrow(x$hpo) > 0) {
    cli::cli_text("{b} HPO: {nrow(x$hpo)} t\u00e9rminos; los m\u00e1s espec\u00edficos: {paste(utils::head(x$hpo$hpo_nombre, 6), collapse = ', ')}")
  }
  if (!is.null(x$proteina)) {
    d <- x$proteina$dominios
    cli::cli_text("{b} Prote\u00edna {x$proteina$id}: {x$proteina$largo} aa{if (nrow(d)) paste0(', dominios ', paste(unique(d$nombre), collapse = ', ')) else ''}")
  }
  if (!is.null(x$clinvar)) {
    k <- x$clinvar$conteo
    cli::cli_text("{b} ClinVar ({x$clinvar$fecha}): {paste(sprintf('%s %d', tolower(k$clasificacion), k$variantes), collapse = ', ')}")
    v <- x$clinvar$variantes
    if (nrow(v) > 0) {
      t <- sort(table(v$tipo), decreasing = TRUE)
      cli::cli_text("{b} Patog\u00e9nicas y probablemente patog\u00e9nicas por tipo: {paste(sprintf('%s %d', names(t), as.integer(t)), collapse = ', ')}")
      if (!is.null(x$proteina) && nrow(x$proteina$dominios) > 0 && any(!is.na(v$pos_proteina))) {
        dom <- .dominio_de(v$pos_proteina, x$proteina$dominios)
        td <- sort(table(dom[!is.na(dom)]), decreasing = TRUE)
        if (length(td)) cli::cli_text("{b} Donde se concentran en la prote\u00edna: {paste(sprintf('%s %d', names(td), as.integer(td)), collapse = ', ')}")
      }
    }
  }
  cli::cli_text("{.emph Grafico: plot(x) o autoplot(x); variantes de ClinVar en x$clinvar$variantes.}")
  invisible(x)
}

`%|NA|%` <- function(a, b) if (length(a) == 0 || is.na(a)) b else a

.dominio_de <- function(pos, dominios) {
  vapply(pos, function(p) {
    if (is.na(p)) return(NA_character_)
    k <- which(dominios$inicio <= p & dominios$fin >= p)
    if (length(k)) dominios$nombre[k[1]] else "fuera de dominios"
  }, character(1))
}
