# Regiones pseudoautosomicas del cromosoma X (1-based, inclusivas).
.par_x <- list(
  GRCh37 = list(c(60001, 2699520), c(154931044, 155260560)),
  GRCh38 = list(c(10001, 2781479), c(155701383, 156030895))
)

# Longitud del cromosoma 1 en cada build, para detectar el build sin encabezado.
.largo_chr1 <- c(GRCh37 = 249250621, GRCh38 = 248956422)

#' Leer un archivo VCF
#'
#' Lee un VCF (comprimido o no) y lo convierte en una tabla larga con una
#' fila por variante, alelo alternativo y muestra. Separa los sitios
#' multialelicos, normaliza los nombres de cromosoma, detecta el build del
#' genoma y marca como hemicigotas las variantes del X en varones.
#'
#' @param ruta Ruta al archivo `.vcf` o `.vcf.gz`.
#' @param muestras Nombres de las muestras a leer. `NULL` lee todas.
#' @param build `"auto"` detecta el build desde el encabezado. Tambien se
#'   puede indicar `"GRCh37"` o `"GRCh38"`.
#' @param sexo Sexo de las muestras, `"M"` o `"F"`. Puede ser un solo valor
#'   para todas o un vector con nombres de muestra, por ejemplo
#'   `c(P01 = "M", P02 = "F")`. `NULL` si no se conoce; en ese caso no se
#'   marcan hemicigotos salvo que el genotipo ya sea haploide.
#'
#' @return Un tibble con las columnas `muestra`, `chr`, `pos`, `id`, `ref`,
#'   `alt`, `genotipo` (`het`, `hom_alt`, `hemi`, `hom_ref` o `NA` si no
#'   se llamo), `qual`, `dp`, `gq`, `ab` (fraccion de lecturas con el alelo
#'   alternativo), `filtro` y `build`. Lleva los atributos `build` y
#'   `archivo`.
#' @export
#' @examples
#' \dontrun{
#' vcf <- iei_leer_vcf("paciente01.vcf.gz", sexo = "M")
#' }
iei_leer_vcf <- function(ruta, muestras = NULL,
                         build = c("auto", "GRCh37", "GRCh38"), sexo = NULL) {
  build <- match.arg(build)
  if (!file.exists(ruta)) {
    cli::cli_abort("No encuentro el archivo {.path {ruta}}.")
  }

  con <- gzfile(ruta, "r")
  lineas <- readLines(con, warn = FALSE)
  close(con)

  meta <- lineas[startsWith(lineas, "##")]
  i_enc <- which(startsWith(lineas, "#CHROM"))
  if (length(i_enc) != 1) {
    cli::cli_abort("{.path {ruta}} no tiene una linea de encabezado {.code #CHROM}; no parece un VCF.")
  }
  encabezado <- strsplit(sub("^#", "", lineas[i_enc]), "\t", fixed = TRUE)[[1]]
  datos <- lineas[-seq_len(i_enc)]
  datos <- datos[nzchar(datos)]

  if (build == "auto") build <- .detectar_build(meta)

  todas <- if (length(encabezado) > 9) encabezado[-(1:9)] else character(0)
  if (length(todas) == 0) {
    cli::cli_abort("El VCF no tiene columnas de muestras, solo sitios.")
  }
  if (is.null(muestras)) muestras <- todas
  faltan <- setdiff(muestras, todas)
  if (length(faltan) > 0) {
    cli::cli_abort(c("Muestras que no estan en el VCF: {.val {faltan}}.",
                     "i" = "Disponibles: {.val {todas}}."))
  }
  sexo <- .preparar_sexo(sexo, muestras)

  if (length(datos) == 0) {
    cli::cli_warn("El VCF no tiene variantes.")
    return(.marcar_vcf(.vcf_vacio(build, ruta), build, ruta, sexo))
  }

  campos <- strsplit(datos, "\t", fixed = TRUE)
  n_col <- lengths(campos)
  if (any(n_col != length(encabezado))) {
    malas <- which(n_col != length(encabezado))[1]
    cli::cli_abort("La variante {malas} tiene {n_col[malas]} columnas y el encabezado {length(encabezado)}.")
  }
  m <- matrix(unlist(campos, use.names = FALSE), ncol = length(encabezado),
              byrow = TRUE, dimnames = list(NULL, encabezado))

  # Expandir sitios multialelicos: una fila por alelo alternativo
  alts <- strsplit(m[, "ALT"], ",", fixed = TRUE)
  fila <- rep(seq_len(nrow(m)), lengths(alts))
  k <- unlist(lapply(lengths(alts), seq_len), use.names = FALSE)
  alt <- unlist(alts, use.names = FALSE)
  util <- !alt %in% c("*", ".", "<*>", "<NON_REF>")
  fila <- fila[util]; k <- k[util]; alt <- alt[util]

  chr <- .normalizar_chr(m[fila, "CHROM"])
  pos <- as.integer(m[fila, "POS"])
  qual <- suppressWarnings(as.numeric(m[fila, "QUAL"]))
  formato <- m[fila, "FORMAT"]
  avisos <- character(0)

  tablas <- lapply(muestras, function(s) {
    valores <- m[fila, s]
    gt <- .extraer_campo(formato, valores, "GT")
    ad <- .extraer_campo(formato, valores, "AD")
    dp <- suppressWarnings(as.numeric(.extraer_campo(formato, valores, "DP")))
    gq <- suppressWarnings(as.numeric(.extraer_campo(formato, valores, "GQ")))

    geno <- .clasificar_genotipo(gt, k)
    ab <- .fraccion_alelica(ad, k)

    es_hemi_x <- identical(sexo[[s]], "M") & chr == "X" & !.en_par(pos, build) &
      !is.na(geno) & geno == "hom_alt"
    geno[es_hemi_x] <- "hemi"

    tibble::tibble(
      muestra = s, chr = chr, pos = pos,
      id = ifelse(m[fila, "ID"] == ".", NA_character_, m[fila, "ID"]),
      ref = m[fila, "REF"], alt = alt, genotipo = geno,
      qual = qual, dp = dp, gq = gq, ab = ab,
      filtro = m[fila, "FILTER"], build = build
    )
  })
  out <- do.call(rbind, tablas)

  for (campo in c("gq", "dp", "ab")) {
    if (all(is.na(out[[campo]]))) avisos <- c(avisos, campo)
  }
  if (length(avisos) > 0) {
    cli::cli_warn("El VCF no trae informacion para {.field {avisos}}; esas columnas quedan en NA.")
  }

  .marcar_vcf(out, build, ruta, sexo)
}

.marcar_vcf <- function(out, build, ruta, sexo) {
  attr(out, "build") <- build
  attr(out, "archivo") <- normalizePath(ruta)
  attr(out, "ieiprio_sexo") <- vapply(sexo, function(v) as.character(v), character(1))
  attr(out, "ieiprio_vcf") <- list(archivo = basename(ruta),
                                   md5 = unname(tools::md5sum(ruta)))
  out
}

.vcf_vacio <- function(build, ruta) {
  out <- tibble::tibble(
    muestra = character(0), chr = character(0), pos = integer(0),
    id = character(0), ref = character(0), alt = character(0),
    genotipo = character(0), qual = numeric(0), dp = numeric(0),
    gq = numeric(0), ab = numeric(0), filtro = character(0),
    build = character(0)
  )
  out
}

.detectar_build <- function(meta) {
  ref <- meta[startsWith(meta, "##reference")]
  if (length(ref) > 0) {
    if (any(grepl("GRCh38|hg38|hs38", ref, ignore.case = TRUE))) return("GRCh38")
    if (any(grepl("GRCh37|hg19|b37|hs37", ref, ignore.case = TRUE))) return("GRCh37")
  }
  contig1 <- meta[grepl("^##contig=<ID=(chr)?1,", meta)]
  if (length(contig1) > 0) {
    largo <- suppressWarnings(as.numeric(sub(".*length=([0-9]+).*", "\\1", contig1[1])))
    b <- names(.largo_chr1)[match(largo, .largo_chr1)]
    if (!is.na(b)) return(b)
  }
  cli::cli_abort(c(
    "No pude detectar el build del genoma desde el encabezado del VCF.",
    "i" = "Indicalo a mano con {.code build = \"GRCh38\"} o {.code build = \"GRCh37\"}."
  ))
}

.normalizar_chr <- function(x) {
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x[toupper(x) %in% c("M", "MT")] <- "MT"
  x[toupper(x) %in% c("X", "Y")] <- toupper(x[toupper(x) %in% c("X", "Y")])
  x
}

.en_par <- function(pos, build) {
  par <- .par_x[[build]]
  (pos >= par[[1]][1] & pos <= par[[1]][2]) | (pos >= par[[2]][1] & pos <= par[[2]][2])
}

.preparar_sexo <- function(sexo, muestras) {
  if (is.null(sexo)) return(stats::setNames(as.list(rep(NA_character_, length(muestras))), muestras))
  sexo <- toupper(sexo)
  if (!all(sexo %in% c("M", "F", NA))) {
    cli::cli_abort("{.arg sexo} solo acepta {.val M} o {.val F}.")
  }
  if (is.null(names(sexo))) {
    if (length(sexo) != 1) {
      cli::cli_abort("Si das mas de un {.arg sexo}, ponle nombres de muestra: {.code c(P01 = \"M\")}.")
    }
    sexo <- stats::setNames(rep(sexo, length(muestras)), muestras)
  }
  faltan <- setdiff(names(sexo), muestras)
  if (length(faltan) > 0) cli::cli_abort("{.arg sexo} nombra muestras que no se leyeron: {.val {faltan}}.")
  out <- stats::setNames(as.list(rep(NA_character_, length(muestras))), muestras)
  out[names(sexo)] <- as.list(unname(sexo))
  out
}

# Extrae un campo del FORMAT (GT, AD, DP...) para cada fila.
# Agrupa por FORMAT unico porque el orden puede cambiar entre filas.
.extraer_campo <- function(formato, valores, campo) {
  out <- rep(NA_character_, length(valores))
  for (f in unique(formato)) {
    idx <- which(formato == f)
    pos <- match(campo, strsplit(f, ":", fixed = TRUE)[[1]])
    if (is.na(pos)) next
    partes <- strsplit(valores[idx], ":", fixed = TRUE)
    out[idx] <- vapply(partes, function(p) if (length(p) >= pos) p[pos] else NA_character_,
                       character(1))
  }
  out[out %in% c(".", "")] <- NA_character_
  out
}

# Clasifica el genotipo respecto al alelo alternativo numero k.
.clasificar_genotipo <- function(gt, k) {
  alelos <- strsplit(ifelse(is.na(gt), ".", gt), "[/|]")
  vapply(seq_along(alelos), function(i) {
    a <- alelos[[i]]
    a <- a[a != "."]
    if (length(a) == 0) return(NA_character_)
    n_k <- sum(a == as.character(k[i]))
    if (length(a) == 1) return(if (n_k == 1) "hemi" else "hom_ref")
    if (n_k == length(a)) return("hom_alt")
    if (n_k >= 1) return("het")
    "hom_ref"
  }, character(1))
}

# Fraccion de lecturas que apoyan al alelo k segun el campo AD.
.fraccion_alelica <- function(ad, k) {
  vapply(seq_along(ad), function(i) {
    if (is.na(ad[i])) return(NA_real_)
    n <- suppressWarnings(as.numeric(strsplit(ad[i], ",", fixed = TRUE)[[1]]))
    total <- sum(n, na.rm = TRUE)
    if (length(n) < k[i] + 1 || total == 0) return(NA_real_)
    n[k[i] + 1] / total
  }, numeric(1))
}
