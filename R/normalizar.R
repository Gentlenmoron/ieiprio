# Normalizacion de variantes: recorte de bases comunes y alineamiento a la
# izquierda (como bcftools norm o vt normalize). Sin esto, una misma indel
# escrita de otra forma no coincide con gnomAD ni con ClinVar.

# Recorta el sufijo comun y luego el prefijo comun, dejando una base de anclaje.
.recortar_alelos <- function(pos, ref, alt) {
  for (k in seq_along(ref)) {
    r <- ref[k]; a <- alt[k]; p <- pos[k]
    while (nchar(r) > 1 && nchar(a) > 1 && substring(r, nchar(r)) == substring(a, nchar(a))) {
      r <- substr(r, 1, nchar(r) - 1); a <- substr(a, 1, nchar(a) - 1)
    }
    while (nchar(r) > 1 && nchar(a) > 1 && substr(r, 1, 1) == substr(a, 1, 1)) {
      r <- substring(r, 2); a <- substring(a, 2); p <- p + 1L
    }
    ref[k] <- r; alt[k] <- a; pos[k] <- p
  }
  list(pos = pos, ref = ref, alt = alt)
}

.es_indel <- function(ref, alt) {
  nchar(ref) != nchar(alt) & !grepl("[^ACGTNacgtn]", paste0(ref, alt))
}

# Normaliza la tabla. El alineamiento a la izquierda necesita la secuencia de
# referencia alrededor de cada indel, que se pide a Ensembl; si no hay
# conexion solo se recortan las bases comunes.
.normalizar_variantes <- function(x, build, margen = 100L) {
  a <- .attrs_ieiprio(x)
  pos0 <- x$pos; ref0 <- x$ref; alt0 <- x$alt
  rec <- .recortar_alelos(x$pos, x$ref, x$alt)
  x$pos <- as.integer(rec$pos); x$ref <- rec$ref; x$alt <- rec$alt
  indel <- which(.es_indel(x$ref, x$alt))
  claves <- unique(paste(x$chr[indel], x$pos[indel], x$ref[indel], x$alt[indel]))
  hechas <- list()
  if (.sin_red()) claves <- character(0)
  for (cl in claves) {
    p <- strsplit(cl, " ", fixed = TRUE)[[1]]
    chr <- p[1]; pos <- as.integer(p[2])
    largo <- max(nchar(p[3]), nchar(p[4]))
    ini <- max(1L, pos - margen)
    ref <- tryCatch(.referencia_consultar(chr, ini, pos + largo + 5L, build),
                    error = function(e) NULL)
    if (is.null(ref$seq)) next
    vars <- data.frame(pos = pos, ref = p[3], alt = p[4])
    hechas[[cl]] <- .normalizar_indels(vars, list(inicio = ini, secuencia = toupper(ref$seq)))
  }
  if (length(hechas) > 0) {
    k <- paste(x$chr, x$pos, x$ref, x$alt)
    for (cl in names(hechas)) {
      i <- which(k == cl)
      x$pos[i] <- as.integer(hechas[[cl]]$pos); x$ref[i] <- hechas[[cl]]$ref; x$alt[i] <- hechas[[cl]]$alt
    }
  }
  cambio <- x$pos != pos0 | x$ref != ref0 | x$alt != alt0
  if (any(cambio)) {
    x$variante_original <- ifelse(cambio, paste0(x$chr, ":", pos0, ":", ref0, ">", alt0), NA_character_)
    cli::cli_inform("{sum(cambio)} variante{?s} normalizada{?s} (recorte o alineamiento a la izquierda).")
  }
  x <- .poner_attrs(x, a)
  if (length(indel) > 0 && .sin_red()) {
    cli::cli_inform(c("i" = "Indels recortadas pero no alineadas a la izquierda (modo sin red)."))
  } else if (length(indel) > 0 && length(hechas) < length(claves)) {
    cli::cli_warn("No se pudo obtener la referencia para alinear {length(claves) - length(hechas)} indel{?s}; solo se recortaron.")
  }
  x
}
