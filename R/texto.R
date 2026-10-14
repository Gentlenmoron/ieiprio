# Entrada de variantes desde texto libre, PDF o Word.

#' Leer el texto de un archivo
#'
#' Extrae el texto de un informe para buscar variantes en el. Acepta texto
#' plano (`.txt`, `.md`, `.csv`, `.tsv`), PDF con texto (no escaneado) y
#' Word (`.docx`).
#'
#' @param archivo Ruta al archivo.
#' @return Un vector de texto, una linea por elemento.
#' @export
#' @examples
#' \dontrun{
#' texto <- iei_leer_texto("informe_secuenciacion.pdf")
#' iei_extraer_variantes(texto)
#' }
iei_leer_texto <- function(archivo) {
  if (!file.exists(archivo)) cli::cli_abort("No encuentro {.path {archivo}}.")
  ext <- tolower(tools::file_ext(archivo))
  lineas <- switch(ext,
    pdf = {
      rlang::check_installed("pdftools", reason = "para leer PDF.")
      unlist(strsplit(pdftools::pdf_text(archivo), "\n", fixed = TRUE))
    },
    docx = {
      rlang::check_installed("xml2", reason = "para leer Word.")
      tmp <- tempfile()
      on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
      utils::unzip(archivo, files = "word/document.xml", exdir = tmp)
      doc <- xml2::read_xml(file.path(tmp, "word", "document.xml"))
      parrafos <- xml2::xml_find_all(doc, "//w:p", xml2::xml_ns(doc))
      vapply(parrafos, function(p) {
        paste(xml2::xml_text(xml2::xml_find_all(p, ".//w:t", xml2::xml_ns(doc))), collapse = "")
      }, character(1))
    },
    doc = cli::cli_abort("El formato {.file .doc} antiguo no se puede leer. Guardalo como {.file .docx} o PDF."),
    readLines(archivo, warn = FALSE, encoding = "UTF-8")
  )
  if (ext == "pdf" && sum(nchar(trimws(lineas))) == 0) {
    cli::cli_warn("El PDF no tiene texto. Probablemente es una imagen escaneada y necesita OCR.")
  }
  lineas
}

# Patrones de nomenclatura. El orden importa: primero lo mas especifico.
.patron_transcrito <- "\\b(?:NM|NR|XM|ENST)_?[0-9]+(?:\\.[0-9]+)?(?:\\([A-Za-z0-9-]+\\))?:[cnr]\\.[^\\s,;]+"
.patron_genomico <- "\\bN[CG]_[0-9]+\\.[0-9]+:g\\.[^\\s,;]+"
.patron_rsid <- "\\brs[0-9]{3,}\\b"
.patron_gen <- paste0("\\b([A-Z][A-Z0-9]{1,9})(?:\\s*:\\s*|\\s+)",
                      "(c\\.[0-9*+_-][^\\s,;]*|p\\.\\(?[A-Z][a-z]{2}[0-9]+[A-Za-z*=]+\\)?)")

.palabras_genotipo <- list(
  hemi = "hemic[i\u00ed]go|hemizyg",
  hom_alt = "homoc[i\u00ed]go|homozyg",
  het = "heteroc[i\u00ed]go|heterozyg"
)

#' Encontrar variantes escritas en un texto
#'
#' Busca descripciones de variantes en texto libre: nomenclatura HGVS de
#' transcrito (`NM_000061.3:c.1632-2A>G`, tambien con el gen entre parentesis
#' como aparece en ClinVar), HGVS genomica (`NC_000023.11:g.101353990T>C`),
#' identificadores rsID y el gen seguido del cambio (`BTK c.1632-2A>G` o
#' `BTK p.Arg562Trp`). Si en la misma linea aparece heterocigoto,
#' homocigoto o hemicigoto, tambien lo registra.
#'
#' @param texto Vector de texto, por ejemplo la salida de [iei_leer_texto()].
#' @return Un tibble con `entrada` (lo encontrado, normalizado para
#'   Ensembl), `tipo`, `genotipo` y `linea`.
#' @export
#' @examples
#' iei_extraer_variantes(c(
#'   "Se identifico la variante NM_000061.3(BTK):c.1632-2A>G en hemicigosis.",
#'   "Tambien rs121908711 y STAT1 p.Arg274Gln en heterocigosis."
#' ))
iei_extraer_variantes <- function(texto) {
  texto <- .normalizar_texto(as.character(texto))
  filas <- list()
  for (i in seq_along(texto)) {
    linea <- texto[i]
    if (is.na(linea) || !nzchar(trimws(linea))) next
    resto <- linea
    agregar <- function(m, tipo) {
      for (e in m) filas[[length(filas) + 1]] <<- data.frame(
        entrada = e, tipo = tipo, linea = i, stringsAsFactors = FALSE)
    }
    for (pt in list(c(.patron_transcrito, "hgvs_transcrito"), c(.patron_genomico, "hgvs_genomico"))) {
      m <- regmatches(resto, gregexpr(pt[1], resto, perl = TRUE))[[1]]
      agregar(.limpiar_hgvs(m), pt[2])
      for (e in m) resto <- sub(e, " ", resto, fixed = TRUE)
    }
    agregar(regmatches(resto, gregexpr(.patron_rsid, resto, perl = TRUE))[[1]], "rsid")
    gm <- regmatches(resto, gregexpr(.patron_gen, resto, perl = TRUE))[[1]]
    if (length(gm) > 0) {
      partes <- regmatches(gm, regexec(.patron_gen, gm, perl = TRUE))
      agregar(vapply(partes, function(p) paste0(p[2], ":", .limpiar_hgvs(p[3])), character(1)),
              "gen_y_cambio")
    }
  }
  out <- do.call(rbind, filas)
  if (is.null(out)) {
    return(tibble::tibble(entrada = character(0), tipo = character(0),
                          genotipo = character(0), linea = integer(0)))
  }
  out$genotipo <- vapply(out$linea, function(i) .genotipo_en_texto(texto[i]), character(1))
  out <- out[!duplicated(out$entrada), ]
  tibble::as_tibble(out[, c("entrada", "tipo", "genotipo", "linea")])
}

# Los informes suelen traer guiones tipograficos, signos de mayor de ancho
# completo y espacios no separables que rompen la nomenclatura.
.normalizar_texto <- function(x) {
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x, perl = TRUE)
  x <- gsub("[\uFF1E\u203A]", ">", x, perl = TRUE)
  x <- gsub("[\u00A0\u2007\u202F]", " ", x, perl = TRUE)
  x <- gsub("&gt;", ">", x, fixed = TRUE)
  x
}

.limpiar_hgvs <- function(x) {
  x <- sub("\\(([A-Za-z0-9-]+)\\):", ":", x)       # NM_000061.3(BTK):c. -> NM_000061.3:c.
  x <- sub("[.)\\]]+$", "", x)                       # punto final o parentesis de cierre
  x <- gsub("p\\.\\((.*)\\)$", "p.\\1", x)           # p.(Arg1Ter) -> p.Arg1Ter
  x
}

.genotipo_en_texto <- function(linea) {
  for (g in names(.palabras_genotipo)) {
    if (grepl(.palabras_genotipo[[g]], linea, ignore.case = TRUE, perl = TRUE)) return(g)
  }
  NA_character_
}

# Consulta el Variant Recoder de Ensembl. Unica funcion que habla con la API.
.recoder_consultar <- function(entradas, build) {
  httr2::request(.vep_servidor(build)) |>
    httr2::req_url_path_append("variant_recoder", "homo_sapiens") |>
    httr2::req_url_query(vcf_string = 1) |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_body_json(list(ids = as.list(entradas))) |>
    .agente() |>
    httr2::req_timeout(120) |>
    httr2::req_throttle(rate = 10 / 1) |>
    httr2::req_retry(max_tries = 4, is_transient = .transitorio) |>
    # Un 400 trae en JSON el motivo (nomenclatura invalida): se lee como
    # variante no reconocida. Solo los errores de servidor detienen el analisis.
    httr2::req_error(is_error = function(resp) httr2::resp_status(resp) >= 500) |>
    .ejecutar_json("Ensembl Variant Recoder")
}

# La respuesta es una lista con un elemento por entrada; cada elemento tiene
# un bloque por alelo alternativo con input, vcf_string, hgvsg, etc.
.recoder_aplanar <- function(resp) {
  if (!is.null(resp$error)) return(.recoder_vacio())
  filas <- list()
  for (el in resp) {
    bloques <- if (!is.null(el$input)) list(el) else Filter(is.list, el)
    for (b in bloques) {
      entrada <- .chr1(b$input)
      vcfs <- unlist(b$vcf_string)
      for (v in vcfs) {
        p <- strsplit(v, "-", fixed = TRUE)[[1]]
        if (length(p) != 4) next
        filas[[length(filas) + 1]] <- data.frame(
          entrada = entrada, chr = .normalizar_chr(p[1]), pos = as.integer(p[2]),
          ref = p[3], alt = p[4], ambigua = FALSE, stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, filas)
  if (is.null(out)) return(.recoder_vacio())
  # Ensembl tambien devuelve coordenadas sobre secuencias LRG y parches; solo
  # interesan los cromosomas del ensamblaje primario.
  out <- unique(out[out$chr %in% c(as.character(1:22), "X", "Y", "MT"), ])
  if (nrow(out) == 0) return(.recoder_vacio())
  n_pos <- table(out$entrada)
  out$ambigua <- as.vector(n_pos[out$entrada] > 1)
  out
}

.recoder_vacio <- function() {
  data.frame(entrada = character(0), chr = character(0), pos = integer(0),
             ref = character(0), alt = character(0), ambigua = logical(0))
}

#' Crear una tabla de variantes desde texto
#'
#' Convierte descripciones de variantes en una tabla con la misma estructura
#' que [iei_leer_vcf()], lista para [iei_anotar()] e [iei_priorizar()]. Sirve
#' cuando no hay VCF, por ejemplo con un informe de otro centro de secuenciacion o una
#' lista de variantes escrita a mano. Las coordenadas se obtienen del
#' Variant Recoder de Ensembl.
#'
#' @param x Un vector con variantes (`"BTK c.1632-2A>G"`, `"rs699"`,
#'   `"NM_000061.3:c.1632-2A>G"`), un texto libre, o la ruta a un archivo
#'   de texto, PDF o Word.
#' @param muestra Codigo de la muestra.
#' @param genotipo Genotipo de las variantes (`"het"`, `"hom_alt"`, `"hemi"`).
#'   Puede ser un valor o uno por variante. Si es `NA` se usa el que diga el
#'   texto, y si no dice nada queda sin indicar.
#' @param sexo `"M"` o `"F"`. En un varon, una variante homocigota en el X
#'   fuera de la region pseudoautosomica pasa a ser hemicigota.
#' @param build `"GRCh38"` o `"GRCh37"`.
#'
#' @return Una tabla como la de [iei_leer_vcf()] con las columnas extra
#'   `gen` (si cae en un gen del panel), `entrada` y `ambigua` (una
#'   descripcion de proteina que corresponde a varias posiciones posibles).
#' @export
#' @examples
#' \dontrun{
#' v <- iei_variantes(c("BTK c.1632-2A>G", "RAG1 p.Arg396Cys"), sexo = "M")
#' iei_priorizar(iei_anotar(v), fenotipo = "HP:0004313")
#'
#' v <- iei_variantes("informe.pdf")
#' }
iei_variantes <- function(x, muestra = "consulta", genotipo = NA, sexo = NULL,
                          build = c("GRCh38", "GRCh37")) {
  build <- match.arg(build)
  if (!is.character(x) || length(x) == 0) cli::cli_abort("{.arg x} debe ser texto.")

  origen <- "texto"
  if (length(x) == 1 && file.exists(x)) {
    origen <- basename(x)
    x <- iei_leer_texto(x)
  }
  encontradas <- iei_extraer_variantes(x)
  if (nrow(encontradas) == 0) {
    cli::cli_abort(c("No encontre variantes en el texto.",
                     "i" = "Usa nomenclatura HGVS ({.val NM_000061.3:c.1632-2A>G}), un rsID o gen y cambio ({.val BTK p.Arg562Trp})."))
  }
  if (!all(is.na(genotipo))) {
    if (length(genotipo) == 1) genotipo <- rep(genotipo, nrow(encontradas))
    if (length(genotipo) != nrow(encontradas)) {
      cli::cli_abort("Di un {.arg genotipo} o uno por variante ({nrow(encontradas)}).")
    }
    encontradas$genotipo <- ifelse(is.na(genotipo), encontradas$genotipo, genotipo)
  }

  # Una consulta por variante: en grupo, Ensembl solo responde la primera.
  coords <- do.call(rbind, lapply(encontradas$entrada, function(e) {
    .recoder_aplanar(.recoder_consultar(e, build))
  }))
  if (is.null(coords)) coords <- .recoder_vacio()
  sin_coord <- setdiff(encontradas$entrada, coords$entrada)
  if (length(sin_coord) > 0) {
    cli::cli_warn("Ensembl no reconocio: {.val {sin_coord}}. Revisa la nomenclatura.")
  }
  if (any(coords$ambigua)) {
    cli::cli_inform(c("i" = "Sin transcrito, {.val {unique(coords$entrada[coords$ambigua])}} corresponde{?n} a varias posiciones. Al anotar se conserva la que coincide con esa nomenclatura en el transcrito MANE.",
                      "i" = "Para evitarlo, indica el transcrito: {.val NM_000061.3:c.1632-2A>G}."))
  }
  if (nrow(coords) == 0) cli::cli_abort("Ninguna variante pudo ubicarse en el genoma.")

  g <- encontradas$genotipo[match(coords$entrada, encontradas$entrada)]
  sexo_m <- identical(toupper(sexo %||% ""), "M")
  g[sexo_m & coords$chr == "X" & !.en_par(coords$pos, build) & g %in% "hom_alt"] <- "hemi"

  tabla <- tibble::tibble(
    muestra = muestra, chr = coords$chr, pos = coords$pos, id = coords$entrada,
    ref = coords$ref, alt = coords$alt, genotipo = g, qual = NA_real_, dp = NA_real_,
    gq = NA_real_, ab = NA_real_, filtro = "PASS", build = build,
    entrada = coords$entrada, ambigua = coords$ambigua
  )
  tabla <- .genes_de_variantes(tabla, build)

  attr(tabla, "build") <- build
  attr(tabla, "ieiprio_sexo") <- stats::setNames(if (is.null(sexo)) NA_character_ else toupper(sexo), muestra)
  attr(tabla, "ieiprio_vcf") <- list(archivo = origen, md5 = .md5_texto(encontradas$entrada))
  attr(tabla, "ieiprio_texto") <- TRUE
  attr(tabla, "ieiprio_bitacora") <- tibble::tibble(
    paso = c("Encontradas en el texto", "Ubicadas en el genoma"),
    filas = c(nrow(encontradas), nrow(tabla)),
    descartadas = c(0L, length(sin_coord)))
  tabla
}

# md5 de un vector de texto, para registrar exactamente que se analizo.
.md5_texto <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(x, f)
  unname(tools::md5sum(f))
}

# Asigna el gen del panel a cada variante; las que no caen en ninguno se
# conservan con gen NA.
.genes_de_variantes <- function(tabla, build, margen = 20) {
  co <- ieiprio::genes_coordenadas[ieiprio::genes_coordenadas$build == build, ]
  tabla$gen <- vapply(seq_len(nrow(tabla)), function(i) {
    k <- which(co$chr == tabla$chr[i] & tabla$pos[i] >= co$inicio - margen &
                 tabla$pos[i] <= co$fin + margen)
    if (length(k) == 0) NA_character_ else paste(sort(unique(co$gen[k])), collapse = ";")
  }, character(1))
  fuera <- is.na(tabla$gen)
  if (any(fuera)) {
    cli::cli_inform(c("i" = "{sum(fuera)} variante{?s} fuera de los genes del panel; se anota{?n} igual."))
  }
  tabla[, c("muestra", "gen", setdiff(names(tabla), c("muestra", "gen")))]
}

#' Consultar una variante al vuelo
#'
#' Atajo para revisar una o pocas variantes sin VCF: las ubica, las anota con
#' Ensembl VEP y las prioriza. Imprime una ficha por variante.
#'
#' @inheritParams iei_variantes
#' @param fenotipo Terminos HPO del paciente, opcional.
#' @return Un `iei_resultado`, de forma invisible.
#' @export
#' @examples
#' \dontrun{
#' iei_variante("BTK c.1632-2A>G", genotipo = "hemi")
#' }
iei_variante <- function(x, genotipo = NA, sexo = NULL, fenotipo = NULL,
                         build = c("GRCh38", "GRCh37")) {
  build <- match.arg(build)
  v <- iei_variantes(x, genotipo = genotipo, sexo = sexo, build = build)
  r <- iei_priorizar(suppressMessages(iei_anotar(v)), fenotipo = fenotipo)
  for (i in seq_len(nrow(r))) .imprimir_ficha(r[i, ])
  invisible(r)
}

.imprimir_ficha <- function(r) {
  panel <- ieiprio::iuis_panel
  info <- panel[match(r$gen, panel$gen), ]
  sin <- function(v, otro = "-") if (length(v) == 0 || is.na(v)) otro else as.character(v)
  titulo <- paste(sin(r$gen, "fuera del panel"), sin(sub("^[^:]+:", "", r$hgvs_c), sin(r$entrada)))
  lineas <- c(
    paste("Proteina:", sin(sub("^[^:]+:", "", r$hgvs_p))),
    paste0("Consecuencia: ", sin(gsub("_", " ", r$consecuencia)), " (", sin(r$impacto), ")"),
    paste0("gnomAD: ", .formato_af(r$af_gnomad), " | ClinVar: ",
           sin(gsub("_", " ", r$clinvar), "sin datos")),
    paste("Predictores:", .texto_predictores(r)),
    paste0("Herencia del gen: ", sin(info$herencia), " | Genotipo: ", sin(r$genotipo, "no indicado")),
    paste("Fenotipos del gen:", .resumir_fenotipos(info$fenotipos))
  )
  cli::cli_h2(titulo)
  cli::cli_bullets(stats::setNames(lineas, rep("*", length(lineas))))
  cli::cli_bullets(c(">" = paste0("Puntaje ", r$puntaje, " (", r$categoria, "): ", r$razones)))
  if (!is.na(r$acmg)) cli::cli_bullets(c(">" = paste("ACMG sugerido:", r$acmg)))
  if (!is.na(r$nota)) cli::cli_bullets(c("!" = r$nota))
}

# PanelApp repite fenotipos con variantes de nombre; se muestran los primeros
# distintos, sin codigos OMIM.
.resumir_fenotipos <- function(f, n = 3) {
  if (length(f) == 0 || is.na(f)) return("-")
  partes <- trimws(strsplit(f, ";", fixed = TRUE)[[1]])
  partes <- trimws(sub(",?\\s*[0-9]{6}$", "", partes))
  partes <- partes[nzchar(partes) & !duplicated(tolower(partes))]
  resto <- length(partes) - n
  paste0(paste(utils::head(partes, n), collapse = "; "),
         if (resto > 0) sprintf(" (y %d mas)", resto) else "")
}
