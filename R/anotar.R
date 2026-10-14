.vep_servidor <- function(build) {
  if (identical(build, "GRCh37")) "https://grch37.rest.ensembl.org" else "https://rest.ensembl.org"
}

.rango_impacto <- c(HIGH = 4, MODERATE = 3, LOW = 2, MODIFIER = 1)

#' Anotar variantes con Ensembl VEP
#'
#' Consulta Ensembl VEP para cada variante unica y agrega su consecuencia,
#' nomenclatura HGVS, exon, frecuencias en gnomAD (global y del grupo mas
#' alto), significancia en ClinVar y predictores. Usa el servidor de GRCh37 o
#' GRCh38 segun el build. Las respuestas se guardan en una cache local por
#' release de Ensembl, asi que repetir el analisis no vuelve a consultar y
#' nunca se mezclan anotaciones de releases distintos.
#'
#' Antes de anotar, las inserciones y deleciones se normalizan (recorte de
#' bases comunes y alineamiento a la izquierda contra la referencia), para
#' que coincidan con gnomAD y ClinVar aunque el VCF las escriba de otra forma.
#'
#' Para cada variante se elige un transcrito, en este orden de preferencia.
#' Primero uno del mismo gen de la columna `gen` (si existe), luego el MANE
#' Select, luego el canonico de Ensembl y por ultimo el de mayor impacto.
#'
#' **Privacidad.** En modo en linea cada variante se envia a Ensembl (Reino
#' Unido) y, con `clinvar = TRUE`, el gen se consulta en NCBI (Estados
#' Unidos). Para no enviar nada, corre VEP en tu computadora con el archivo
#' de [iei_vep_entrada()] y pasa su salida en `vep_json`, con
#' `clinvar = FALSE`. Con `options(ieiprio.sin_red = TRUE)` el paquete se
#' niega a hacer cualquier consulta externa.
#'
#' @param x Tabla devuelta por [iei_filtrar()], [iei_sanger()],
#'   [iei_variantes()] o [iei_leer_vcf()].
#' @param cache Carpeta de la cache, o `FALSE` para no usarla.
#' @param lote Variantes por consulta. La API acepta hasta 200.
#' @param normalizar Si `TRUE`, normaliza inserciones y deleciones.
#' @param clinvar Si `TRUE`, descarga las variantes patogenicas de ClinVar de
#'   cada gen afectado para obtener las estrellas de revision y los criterios
#'   PS1, PM5 y PM1 (ver [iei_priorizar()]).
#' @param gnomad Si `TRUE`, consulta en gnomAD los homocigotos, hemicigotos,
#'   la frecuencia filtrada al 95 % y la restriccion del gen. Es lento (gnomAD
#'   limita las consultas); conviene usarlo despues de filtrar, o aplicar
#'   [iei_gnomad()] solo a las candidatas.
#' @param vep_json Ruta a la salida JSON de un VEP local, en lugar de la API.
#'
#' @return La tabla `x` con las columnas `transcrito`, `mane`, `exon`,
#'   `consecuencia`, `impacto`, `hgvs_c`, `hgvs_p`, `af_gnomad` (global),
#'   `af_grpmax` y `grupo_grpmax` (grupo poblacional mas alto), `clinvar`,
#'   `proteina_id`, `pos_proteina`, los predictores `alphamissense`,
#'   `am_clase`, `revel` y `spliceai`, y `vep_release`. Con `clinvar = TRUE`
#'   agrega `clinvar_estrellas`, `ps1`, `pm5` y `pm1`; con `gnomad = TRUE`,
#'   las columnas de [iei_gnomad()]. REVEL es de uso no comercial.
#' @export
#' @examples
#' \dontrun{
#' anot <- iei_anotar(iei_filtrar(iei_leer_vcf("paciente01.vcf.gz", sexo = "M")))
#' }
iei_anotar <- function(x, cache = tools::R_user_dir("ieiprio", "cache"), lote = 200,
                       normalizar = TRUE, clinvar = TRUE, gnomad = FALSE, vep_json = NULL) {
  necesarias <- c("chr", "pos", "ref", "alt")
  if (!is.data.frame(x) || !all(necesarias %in% names(x))) {
    cli::cli_abort("{.arg x} debe venir de {.fn iei_leer_vcf} o {.fn iei_filtrar}.")
  }
  if (lote < 1 || lote > 200) cli::cli_abort("{.arg lote} debe estar entre 1 y 200.")
  build <- attr(x, "build") %||% unique(x$build)[1]
  if (is.null(build) || is.na(build)) cli::cli_abort("No se sabe el build de {.arg x}.")

  a <- .attrs_ieiprio(x)
  vacias <- c("transcrito", "mane", "exon", "consecuencia", "impacto", "hgvs_c",
              "hgvs_p", "af_gnomad", "af_grpmax", "grupo_grpmax", "clinvar", "proteina_id",
              "pos_proteina", "alphamissense", "am_clase", "revel", "spliceai")
  numericas <- c("af_gnomad", "af_grpmax", "pos_proteina", "alphamissense", "revel", "spliceai")
  if (nrow(x) == 0) {
    for (v in vacias) x[[v]] <- if (v %in% numericas) numeric(0) else character(0)
    x$vep_release <- character(0)
    x <- .poner_attrs(x, a)
    attr(x, "vep_release") <- NA_character_
    return(x)
  }

  if (normalizar && is.null(vep_json)) {
    x <- .normalizar_variantes(x, build)
    a <- .attrs_ieiprio(x)
  }

  clave <- paste(build, x$chr, x$pos, x$ref, x$alt, sep = "_")
  unicas <- unique(clave)

  if (!is.null(vep_json)) {
    release <- "local"
    tabla <- .vep_local(vep_json, unicas)
  } else {
    .aviso_privacidad("Ensembl VEP")
    release <- .vep_release(build)
    if (is.na(release)) release <- .release_en_cache(cache, build)
    tabla <- .cache_leer(cache, build, release)
    faltan <- setdiff(unicas, tabla$clave)
    if (length(faltan) > 0) {
      cli::cli_inform("Consultando Ensembl VEP ({build}) para {length(faltan)} variante{?s}.")
      entrada <- .entrada_vep(faltan)
      grupos <- split(seq_along(faltan), ceiling(seq_along(faltan) / lote))
      nuevas <- lapply(grupos, function(idx) {
        resp <- .vep_consultar(entrada[idx], build)
        .vep_aplanar(resp, stats::setNames(faltan[idx], entrada[idx]))
      })
      nuevas <- do.call(rbind, nuevas)
      sin_resp <- setdiff(faltan, nuevas$clave)
      if (length(sin_resp) > 0) nuevas <- rbind(nuevas, .vep_fila_vacia(sin_resp))
      tabla <- rbind(tabla, nuevas)
      .cache_guardar(tabla, cache, build, release)
    }
  }

  anot <- .vep_elegir(clave, if ("gen" %in% names(x)) x$gen else rep(NA_character_, nrow(x)), tabla)
  for (v in vacias) x[[v]] <- anot[[v]]
  x$vep_release <- rep(release, nrow(x))
  if ("id" %in% names(x)) x$id <- ifelse(is.na(x$id), anot$rsid, x$id)

  x <- .poner_attrs(x, a)
  x <- .resolver_ambiguas(x)
  attr(x, "build") <- build
  attr(x, "vep_release") <- release
  if (clinvar) x <- .anotar_clinvar(x, cache)
  if (gnomad) x <- iei_gnomad(x, cache = cache)
  x
}

# Entrada de VEP en formato por defecto: "chr pos . ref alt . . .".
.entrada_vep <- function(claves) {
  partes <- strsplit(claves, "_", fixed = TRUE)
  vapply(partes, function(p) paste(p[2], p[3], ".", p[4], p[5], ". . ."), character(1))
}

#' Archivo de entrada para un VEP local
#'
#' Escribe las variantes en el formato por defecto de VEP, para anotarlas en
#' tu computadora sin enviarlas a internet. Despues corre VEP con salida JSON
#' y los mismos complementos que usa el paquete, y pasa el resultado a
#' [iei_anotar()] con `vep_json`.
#'
#' Un comando de ejemplo (requiere VEP con su cache y los complementos
#' instalados):
#'
#' `vep -i entrada.txt --json -o salida.json --cache --offline --assembly GRCh38
#' --hgvs --mane --canonical --protein --numbers --af_gnomade --af_gnomadg
#' --check_existing --plugin AlphaMissense,file=... --plugin REVEL,file=...
#' --plugin SpliceAI,snv=...,indel=...`
#'
#' @param x Tabla con `chr`, `pos`, `ref` y `alt`.
#' @param archivo Ruta del archivo a escribir.
#' @return La ruta, de forma invisible.
#' @export
iei_vep_entrada <- function(x, archivo) {
  build <- attr(x, "build") %||% unique(x$build)[1] %||% "GRCh38"
  claves <- unique(paste(build, x$chr, x$pos, x$ref, x$alt, sep = "_"))
  writeLines(sub(" \\. \\. \\.$", "", .entrada_vep(claves)), archivo)
  cli::cli_inform("{length(claves)} variante{?s} escritas en {.path {archivo}}.")
  invisible(archivo)
}

# Lee la salida JSON de un VEP local (un objeto por linea o un arreglo).
.vep_local <- function(ruta, claves) {
  if (!file.exists(ruta)) cli::cli_abort("No encuentro {.path {ruta}}.")
  lineas <- readLines(ruta, warn = FALSE)
  lineas <- lineas[nzchar(trimws(lineas))]
  resp <- if (length(lineas) > 0 && startsWith(trimws(lineas[1]), "[")) {
    jsonlite::fromJSON(paste(lineas, collapse = "\n"), simplifyVector = FALSE)
  } else {
    lapply(lineas, jsonlite::fromJSON, simplifyVector = FALSE)
  }
  entrada <- .entrada_vep(claves)
  # VEP local devuelve la entrada sin los puntos finales; se aceptan ambas formas.
  mapa <- c(stats::setNames(claves, entrada), stats::setNames(claves, sub(" \\. \\. \\.$", "", entrada)))
  tabla <- .vep_aplanar(resp, mapa)
  sin <- setdiff(claves, tabla$clave)
  if (length(sin) > 0) {
    cli::cli_warn("La salida de VEP no trae {length(sin)} variante{?s}; quedan sin anotar.")
    tabla <- rbind(tabla, .vep_fila_vacia(sin))
  }
  tabla
}

# Una descripcion sin transcrito ("BTK c.1632-2A>G") puede corresponder a
# varias posiciones del genoma, una por transcrito. Despues de anotar se
# conserva la posicion cuya nomenclatura en el transcrito elegido (MANE si
# existe) coincide con lo escrito. Si ninguna coincide, se conservan todas.
.resolver_ambiguas <- function(x) {
  if (!all(c("entrada", "ambigua") %in% names(x)) || !any(x$ambigua %in% TRUE)) return(x)
  cambio <- sub("^.*:", "", x$entrada)
  hc <- sub("^[^:]+:", "", .col(x, "hgvs_c", NA_character_))
  hp <- sub("^[^:]+:", "", .col(x, "hgvs_p", NA_character_))
  coincide <- (!is.na(hc) & hc == cambio) | (!is.na(hp) & hp == cambio)
  conservar <- rep(TRUE, nrow(x))
  for (e in unique(x$entrada[x$ambigua %in% TRUE])) {
    idx <- which(x$entrada == e)
    if (any(coincide[idx])) conservar[idx] <- coincide[idx]
  }
  a <- .attrs_ieiprio(x)
  out <- .poner_attrs(x[conservar, ], a)
  out$ambigua[out$entrada %in% x$entrada[!conservar]] <- FALSE
  descartadas <- sum(!conservar)
  if (descartadas > 0 && !is.null(a$ieiprio_bitacora)) {
    attr(out, "ieiprio_bitacora") <- rbind(a$ieiprio_bitacora, tibble::tibble(
      paso = "Posiciones ambiguas resueltas", filas = nrow(out), descartadas = descartadas))
  }
  out
}

# Unica funcion que habla con la API. En los tests se reemplaza.
.vep_consultar <- function(entrada, build) {
  req <- httr2::request(.vep_servidor(build)) |>
    httr2::req_url_path_append("vep", "homo_sapiens", "region") |>
    httr2::req_url_query(canonical = 1, hgvs = 1, mane = 1, protein = 1, numbers = 1,
                         af_gnomade = 1, af_gnomadg = 1,
                         AlphaMissense = 1, REVEL = 1, SpliceAI = 2) |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_body_json(list(variants = as.list(entrada))) |>
    .agente() |>
    httr2::req_timeout(180) |>
    httr2::req_throttle(rate = 10 / 1) |>
    httr2::req_retry(max_tries = 5, is_transient = .transitorio)
  .ejecutar_json(req, "Ensembl VEP")
}

.vep_release <- function(build) {
  if (.sin_red()) return(NA_character_)
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
    grpmax <- list(af = NA_real_, grupo = NA_character_)
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
        g <- .grpmax(f)
        if (!is.na(g$af) && (is.na(grpmax$af) || g$af > grpmax$af)) grpmax <- g
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
        af_gnomad = af, af_grpmax = grpmax$af, grupo_grpmax = grpmax$grupo, exon = NA_character_,
        clinvar = clinvar, rsid = rsid,
        proteina_id = NA_character_, pos_proteina = NA_real_,
        alphamissense = NA_real_, am_clase = NA_character_, revel = NA_real_,
        spliceai = NA_real_, stringsAsFactors = FALSE
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
      af_gnomad = af, af_grpmax = grpmax$af, grupo_grpmax = grpmax$grupo, exon = .chr1(t$exon),
      clinvar = clinvar, rsid = rsid,
      proteina_id = .chr1(t$protein_id),
      pos_proteina = suppressWarnings(as.numeric(.chr1(t$protein_start))),
      alphamissense = .predictores(t)$alphamissense, am_clase = .predictores(t)$am_clase,
      revel = .predictores(t)$revel, spliceai = .predictores(t)$spliceai,
      stringsAsFactors = FALSE
    )))
  })
  out <- do.call(rbind, filas)
  if (is.null(out)) out <- .vep_fila_vacia(character(0))
  out
}

# Frecuencia del grupo poblacional mas alto (grpmax), como recomienda
# ClinGen: se excluyen los grupos con efecto fundador o mezclados (judios
# asquenazies, fineses, amish, Medio Oriente y "otros" o "resto").
.grupos_excluidos_grpmax <- c("asj", "fin", "ami", "mid", "oth", "remaining")
.grpmax <- function(f) {
  nombres <- names(f)
  sub_ <- nombres[grepl("^gnomad[eg]_", nombres)]
  grupo <- sub("^gnomad[eg]_", "", sub_)
  sub_ <- sub_[!grupo %in% .grupos_excluidos_grpmax]
  if (length(sub_) == 0) return(list(af = NA_real_, grupo = NA_character_))
  v <- suppressWarnings(as.numeric(unlist(f[sub_])))
  if (all(is.na(v))) return(list(af = NA_real_, grupo = NA_character_))
  k <- which.max(v)
  list(af = v[k], grupo = sub("^gnomad[eg]_", "", sub_[k]))
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
             af_gnomad = rep(NA_real_, n), af_grpmax = rep(NA_real_, n),
             grupo_grpmax = rep(NA_character_, n), exon = rep(NA_character_, n),
             clinvar = rep(NA_character_, n),
             rsid = rep(NA_character_, n), proteina_id = rep(NA_character_, n),
             pos_proteina = rep(NA_real_, n), alphamissense = rep(NA_real_, n),
             am_clase = rep(NA_character_, n), revel = rep(NA_real_, n),
             spliceai = rep(NA_real_, n), stringsAsFactors = FALSE)
}

.num1 <- function(x) suppressWarnings(as.numeric(.chr1(x)))

# Lee los predictores de los plugins de VEP de un transcrito. Acepta las
# variantes de nombre que usa la API segun la version.
.predictores <- function(t) {
  am <- t$alphamissense %||% list()
  am_score <- if (is.list(am)) .num1(am$am_pathogenicity) else .num1(am)
  am_clase <- if (is.list(am)) .chr1(am$am_class) else NA_character_
  revel <- .num1(t$revel %||% t$revel_score)
  sa <- t$spliceai %||% list()
  ds <- suppressWarnings(as.numeric(unlist(sa[intersect(names(sa), c("DS_AG", "DS_AL", "DS_DG", "DS_DL"))])))
  ds <- ds[!is.na(ds)]
  list(alphamissense = am_score, am_clase = am_clase, revel = revel,
       spliceai = if (length(ds) == 0) NA_real_ else max(ds))
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

# La cache se separa por build y por release de Ensembl: una actualizacion
# de Ensembl no mezcla anotaciones nuevas con viejas. El numero de version del
# formato cambia cuando cambian las columnas guardadas.
.cache_ruta <- function(cache, build, release) {
  file.path(cache, sprintf("vep_v4_%s_%s.rds", build, ifelse(is.na(release), "desconocido", release)))
}

# Sin conexion no se sabe el release actual: se usa la cache mas reciente.
.release_en_cache <- function(cache, build) {
  if (isFALSE(cache)) return(NA_character_)
  f <- list.files(cache, pattern = sprintf("^vep_v4_%s_.*\\.rds$", build))
  if (length(f) == 0) return(NA_character_)
  r <- sub(sprintf("^vep_v4_%s_(.*)\\.rds$", build), "\\1", f)
  r <- r[r != "desconocido"]
  if (length(r) == 0) return(NA_character_)
  as.character(max(suppressWarnings(as.numeric(r)), na.rm = TRUE))
}

.cache_leer <- function(cache, build, release) {
  if (isFALSE(cache)) return(.vep_fila_vacia(character(0)))
  ruta <- .cache_ruta(cache, build, release)
  if (!file.exists(ruta)) return(.vep_fila_vacia(character(0)))
  readRDS(ruta)
}

.cache_guardar <- function(tabla, cache, build, release) {
  if (isFALSE(cache)) return(invisible(NULL))
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  saveRDS(tabla, .cache_ruta(cache, build, release))
}

#' Borrar la cache de anotaciones
#'
#' @param cache Carpeta de la cache.
#' @export
iei_limpiar_cache <- function(cache = tools::R_user_dir("ieiprio", "cache")) {
  rutas <- list.files(cache, pattern = "^(vep|clinvar|proteina|gnomad|ref)_.*\\.rds$", full.names = TRUE)
  unlink(rutas)
  cli::cli_inform("Se borraron {length(rutas)} archivo{?s} de cach\u00e9.")
  invisible(rutas)
}

# Frecuencia que se compara con los umbrales: FAF95 > grpmax > global.
.af_referencia <- function(x) {
  af <- .col(x, "af_gnomad")
  g <- .col(x, "af_grpmax")
  af <- ifelse(!is.na(g), g, af)
  f <- .col(x, "faf95")
  ifelse(!is.na(f), f, af)
}

#' Filtrar por frecuencia poblacional segun el modo de herencia
#'
#' Descarta variantes demasiado frecuentes en gnomAD para causar una
#' enfermedad rara. El umbral depende de la herencia del gen, porque una
#' variante recesiva puede ser mas frecuente que una dominante. Se compara la
#' frecuencia filtrada al 95 % si se consulto gnomAD, si no la del grupo
#' poblacional mas alto (grpmax) y, si tampoco esta, la global. Las variantes
#' ausentes de gnomAD siempre pasan.
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
  af <- .af_referencia(x)
  conservar <- is.na(af) | af <= umbral

  a <- .attrs_ieiprio(x)
  out <- .poner_attrs(x[conservar, ], a)
  if (!is.null(a$ieiprio_bitacora)) {
    attr(out, "ieiprio_bitacora") <- rbind(
      a$ieiprio_bitacora,
      tibble::tibble(paso = "Frecuencia", filas = nrow(out), descartadas = nrow(x) - nrow(out)))
  }
  attr(out, "ieiprio_max_af") <- max_af
  out
}
