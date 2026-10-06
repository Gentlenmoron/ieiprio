.vep_servidor <- function(build) {
  if (identical(build, "GRCh37")) "https://grch37.rest.ensembl.org" else "https://rest.ensembl.org"
}

.rango_impacto <- c(HIGH = 4, MODERATE = 3, LOW = 2, MODIFIER = 1)

#' Anotar variantes con Ensembl VEP
#'
#' Consulta la API REST de Ensembl VEP para cada variante unica y agrega su
#' consecuencia, nomenclatura HGVS, frecuencia en gnomAD y significancia
#' clinica en ClinVar. Usa el servidor de GRCh37 o GRCh38 segun el build
#' del VCF. Las respuestas se guardan en una cache local, asi que volver a
#' anotar las mismas variantes no consulta la API de nuevo.
#'
#' Para cada variante se elige un transcrito, en este orden de preferencia.
#' Primero uno del mismo gen de la columna `gen` (si existe), luego el MANE
#' Select, luego el canonico de Ensembl y por ultimo el de mayor impacto.
#'
#' @param x Tabla devuelta por [iei_filtrar()] o [iei_leer_vcf()].
#' @param cache Carpeta de la cache, o `FALSE` para no usarla.
#' @param lote Variantes por consulta. La API acepta hasta 200.
#'
#' @return La tabla `x` con las columnas `transcrito`, `mane`,
#'   `consecuencia`, `impacto`, `hgvs_c`, `hgvs_p`, `af_gnomad` y
#'   `clinvar`. Si la columna `id` estaba vacia, se completa con el rsID de
#'   Ensembl. Conserva la bitacora y agrega el atributo `vep_release`.
#' @export
#' @examples
#' \dontrun{
#' anot <- iei_anotar(iei_filtrar(iei_leer_vcf("paciente01.vcf.gz", sexo = "M")))
#' }
iei_anotar <- function(x, cache = tools::R_user_dir("ieiprio", "cache"), lote = 200) {
  necesarias <- c("chr", "pos", "ref", "alt")
  if (!is.data.frame(x) || !all(necesarias %in% names(x))) {
    cli::cli_abort("{.arg x} debe venir de {.fn iei_leer_vcf} o {.fn iei_filtrar}.")
  }
  if (lote < 1 || lote > 200) cli::cli_abort("{.arg lote} debe estar entre 1 y 200.")
  build <- attr(x, "build") %||% unique(x$build)[1]
  if (is.null(build) || is.na(build)) cli::cli_abort("No se sabe el build de {.arg x}.")

  bitacora <- attr(x, "ieiprio_bitacora")
  vacias <- c("transcrito", "mane", "consecuencia", "impacto", "hgvs_c",
              "hgvs_p", "af_gnomad", "clinvar")
  if (nrow(x) == 0) {
    for (v in vacias) x[[v]] <- if (v == "af_gnomad") numeric(0) else character(0)
    return(.restaurar_atributos(x, build, bitacora, NA_character_))
  }

  clave <- paste(build, x$chr, x$pos, x$ref, x$alt, sep = "_")
  unicas <- unique(clave)

  tabla_cache <- .cache_leer(cache, build)
  faltan <- setdiff(unicas, tabla_cache$clave)
  release <- attr(tabla_cache, "vep_release") %||% NA_character_

  if (length(faltan) > 0) {
    cli::cli_inform("Consultando Ensembl VEP ({build}) para {length(faltan)} variante{?s}.")
    partes <- strsplit(faltan, "_", fixed = TRUE)
    entrada <- vapply(partes, function(p) paste(p[2], p[3], ".", p[4], p[5], ". . ."), character(1))
    grupos <- split(seq_along(faltan), ceiling(seq_along(faltan) / lote))
    nuevas <- lapply(grupos, function(idx) {
      resp <- .vep_consultar(entrada[idx], build)
      .vep_aplanar(resp, stats::setNames(faltan[idx], entrada[idx]))
    })
    nuevas <- do.call(rbind, nuevas)
    sin_resp <- setdiff(faltan, nuevas$clave)
    if (length(sin_resp) > 0) {
      nuevas <- rbind(nuevas, .vep_fila_vacia(sin_resp))
    }
    release <- .vep_release(build)
    tabla_cache <- rbind(tabla_cache, nuevas)
    attr(tabla_cache, "vep_release") <- release
    .cache_guardar(tabla_cache, cache, build)
  }

  anot <- .vep_elegir(clave, if ("gen" %in% names(x)) x$gen else rep(NA_character_, nrow(x)),
                      tabla_cache)
  for (v in vacias) x[[v]] <- anot[[v]]
  if ("id" %in% names(x)) x$id <- ifelse(is.na(x$id), anot$rsid, x$id)

  .restaurar_atributos(x, build, bitacora, release)
}

.restaurar_atributos <- function(x, build, bitacora, release) {
  attr(x, "build") <- build
  if (!is.null(bitacora)) attr(x, "ieiprio_bitacora") <- bitacora
  attr(x, "vep_release") <- release
  x
}

# Unica funcion que habla con la API. En los tests se reemplaza.
.vep_consultar <- function(entrada, build) {
  req <- httr2::request(.vep_servidor(build)) |>
    httr2::req_url_path_append("vep", "homo_sapiens", "region") |>
    httr2::req_url_query(canonical = 1, hgvs = 1, mane = 1,
                         af_gnomade = 1, af_gnomadg = 1) |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_body_json(list(variants = as.list(entrada))) |>
    .agente() |>
    httr2::req_timeout(180) |>
    httr2::req_throttle(rate = 10 / 1) |>
    httr2::req_retry(max_tries = 5,
                     is_transient = function(resp) httr2::resp_status(resp) %in% c(429, 503))
  httr2::resp_body_json(httr2::req_perform(req), simplifyVector = FALSE)
}

.vep_release <- function(build) {
  tryCatch({
    r <- httr2::request(.vep_servidor(build)) |>
      httr2::req_url_path_append("info", "software") |>
      httr2::req_headers(Accept = "application/json") |>
      .agente() |>
      httr2::req_timeout(30) |>
      httr2::req_perform() |>
      httr2::resp_body_json()
    as.character(r$release)
  }, error = function(e) NA_character_)
}

# Convierte la respuesta de VEP en una tabla larga: una fila por variante y
# transcrito, con los datos de la variante repetidos.
# `mapa` es un vector con nombres: entrada enviada -> clave interna.
.vep_aplanar <- function(resp, mapa) {
  filas <- lapply(resp, function(o) {
    entrada <- .chr1(o$input)
    clave <- unname(mapa[entrada])
    if (is.na(clave)) return(NULL)
    alt <- strsplit(entrada, " ", fixed = TRUE)[[1]][5]

    co <- .o(o$colocated_variants, list())
    rsid <- NA_character_
    clinvar <- NA_character_
    af <- NA_real_
    for (cv in co) {
      id <- .chr1(cv$id)
      if (is.na(rsid) && !is.na(id) && startsWith(id, "rs")) rsid <- id
      cs <- .clinvar_alelo(cv, alt)
      if (is.na(clinvar) && !is.na(cs)) clinvar <- cs
      f <- cv$frequencies[[alt]]
      if (!is.null(f)) {
        valores <- suppressWarnings(as.numeric(unlist(f[c("gnomade", "gnomadg")])))
        valores <- valores[!is.na(valores)]
        if (length(valores) > 0) af <- max(c(af, valores), na.rm = TRUE)
      }
    }

    tc <- .o(o$transcript_consequences, list())
    tc <- Filter(function(t) is.null(t$variant_allele) || identical(t$variant_allele, alt), tc)
    if (length(tc) == 0) {
      return(data.frame(
        clave = clave, gen_vep = NA_character_, transcrito = NA_character_,
        mane = NA_character_, canonico = FALSE,
        consecuencia = .chr1(o$most_severe_consequence), impacto = NA_character_,
        hgvs_c = NA_character_, hgvs_p = NA_character_,
        af_gnomad = af, clinvar = clinvar, rsid = rsid, stringsAsFactors = FALSE
      ))
    }
    do.call(rbind, lapply(tc, function(t) data.frame(
      clave = clave,
      gen_vep = .chr1(t$gene_symbol),
      transcrito = .chr1(t$transcript_id),
      mane = .chr1(t$mane_select),
      canonico = isTRUE(as.logical(.o(t$canonical, 0))),
      consecuencia = paste(unlist(t$consequence_terms), collapse = "&"),
      impacto = .chr1(t$impact),
      hgvs_c = .chr1(t$hgvsc), hgvs_p = .chr1(t$hgvsp),
      af_gnomad = af, clinvar = clinvar, rsid = rsid, stringsAsFactors = FALSE
    )))
  })
  out <- do.call(rbind, filas)
  if (is.null(out)) out <- .vep_fila_vacia(character(0))
  out
}

# Significancia ClinVar para el alelo alternativo, si VEP la separa por alelo.
.clinvar_alelo <- function(cv, alt) {
  por_alelo <- .chr1(cv$clin_sig_allele)
  if (!is.na(por_alelo)) {
    pares <- strsplit(strsplit(por_alelo, ";", fixed = TRUE)[[1]], ":", fixed = TRUE)
    sig <- unique(vapply(Filter(function(p) length(p) == 2 && p[1] == alt, pares),
                         function(p) p[2], character(1)))
    return(if (length(sig) == 0) NA_character_ else paste(sig, collapse = ";"))
  }
  cs <- unlist(cv$clin_sig)
  if (length(cs) == 0) NA_character_ else paste(unique(cs), collapse = ";")
}

.vep_fila_vacia <- function(claves) {
  n <- length(claves)
  data.frame(clave = claves, gen_vep = rep(NA_character_, n), transcrito = rep(NA_character_, n),
             mane = rep(NA_character_, n), canonico = rep(FALSE, n),
             consecuencia = rep(NA_character_, n), impacto = rep(NA_character_, n),
             hgvs_c = rep(NA_character_, n), hgvs_p = rep(NA_character_, n),
             af_gnomad = rep(NA_real_, n), clinvar = rep(NA_character_, n),
             rsid = rep(NA_character_, n), stringsAsFactors = FALSE)
}

# Elige un transcrito por fila: mismo gen, MANE, canonico, mayor impacto.
.vep_elegir <- function(clave, gen, tabla) {
  por_clave <- split(seq_len(nrow(tabla)), tabla$clave)
  filas <- vapply(seq_along(clave), function(i) {
    idx <- por_clave[[clave[i]]]
    if (is.null(idx)) return(NA_integer_)
    t <- tabla[idx, , drop = FALSE]
    mismo_gen <- !is.na(gen[i]) & !is.na(t$gen_vep) & t$gen_vep == gen[i]
    con_mane <- !is.na(t$mane)
    impacto <- ifelse(is.na(t$impacto), 0, .rango_impacto[t$impacto])
    puntaje <- 1000 * mismo_gen + 100 * con_mane + 10 * t$canonico + impacto
    idx[which.max(puntaje)]
  }, integer(1))
  tabla[filas, , drop = FALSE]
}

.cache_ruta <- function(cache, build) file.path(cache, sprintf("vep_%s.rds", build))

.cache_leer <- function(cache, build) {
  if (isFALSE(cache)) return(.vep_fila_vacia(character(0)))
  ruta <- .cache_ruta(cache, build)
  if (!file.exists(ruta)) return(.vep_fila_vacia(character(0)))
  readRDS(ruta)
}

.cache_guardar <- function(tabla, cache, build) {
  if (isFALSE(cache)) return(invisible(NULL))
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  saveRDS(tabla, .cache_ruta(cache, build))
}

#' Borrar la cache de anotaciones
#'
#' @param cache Carpeta de la cache.
#' @export
iei_limpiar_cache <- function(cache = tools::R_user_dir("ieiprio", "cache")) {
  rutas <- list.files(cache, pattern = "^vep_.*\\.rds$", full.names = TRUE)
  unlink(rutas)
  cli::cli_inform("Se borraron {length(rutas)} archivo{?s} de cache.")
  invisible(rutas)
}

#' Filtrar por frecuencia poblacional segun el modo de herencia
#'
#' Descarta variantes demasiado frecuentes en gnomAD para causar una
#' enfermedad rara. El umbral depende de la herencia del gen, porque una
#' variante recesiva puede ser mas frecuente que una dominante. Las
#' variantes ausentes de gnomAD siempre pasan.
#'
#' @param x Tabla devuelta por [iei_anotar()], con columna `gen`.
#' @param max_af Umbrales por modo de herencia. Los genes sin herencia
#'   conocida usan el de `AR`.
#'
#' @return La tabla filtrada, con el paso "Frecuencia" agregado a la bitacora.
#' @export
iei_filtrar_frecuencia <- function(x, max_af = c(AD = 1e-4, XL = 1e-4, AR = 0.01,
                                                 `AR/AD` = 0.01, MT = 0.01)) {
  if (!all(c("gen", "af_gnomad") %in% names(x))) {
    cli::cli_abort("{.arg x} debe venir de {.fn iei_anotar} sobre una tabla de {.fn iei_filtrar}.")
  }
  herencia <- ieiprio::iuis_panel$herencia[match(x$gen, ieiprio::iuis_panel$gen)]
  umbral <- unname(max_af[herencia])
  umbral[is.na(umbral)] <- max_af[["AR"]]
  conservar <- is.na(x$af_gnomad) | x$af_gnomad <= umbral

  atributos <- attributes(x)[c("build", "ieiprio_bitacora", "vep_release")]
  out <- x[conservar, ]
  b <- atributos$ieiprio_bitacora
  if (!is.null(b)) {
    b <- rbind(b, tibble::tibble(paso = "Frecuencia", filas = nrow(out),
                                 descartadas = nrow(x) - nrow(out)))
  }
  .restaurar_atributos(out, atributos$build, b, atributos$vep_release %||% NA_character_)
}
