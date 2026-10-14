# Usa las variantes patogenicas de ClinVar de cada gen para completar la
# anotacion: estrellas de revision de la propia variante y los criterios
# ACMG que dependen de otras variantes del mismo gen.
#
# PS1  mismo cambio de aminoacido que una patogenica conocida, por otro
#      cambio de nucleotido.
# PM5  otro cambio missense en un codon donde ya hay una missense patogenica.
# PM1  (de apoyo) zona caliente: al menos 3 posiciones con missense
#      patogenicas a 5 aminoacidos o menos.

.aa3 <- "Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Sec|Ter"

# Separa un cambio p. de tres letras en aminoacido original, posicion y nuevo.
.partir_p <- function(p) {
  p <- sub("^.*p\\.", "", ifelse(is.na(p), "", p))
  p <- gsub("[()]", "", p)
  patron <- sprintf("^(%s)([0-9]+)(%s|=)$", .aa3, .aa3)
  ok <- grepl(patron, p)
  data.frame(ref = ifelse(ok, sub(patron, "\\1", p), NA_character_),
             pos = suppressWarnings(as.integer(ifelse(ok, sub(patron, "\\2", p), NA))),
             alt = ifelse(ok, sub(patron, "\\3", p), NA_character_),
             stringsAsFactors = FALSE)
}

.c_de_titulo <- function(titulo) {
  ok <- grepl("c\\.[^ ]+", titulo %||% "")
  ifelse(ok, sub("^.*?(c\\.[^ ]+).*$", "\\1", titulo, perl = TRUE), NA_character_)
}

.p_de_titulo <- function(titulo) {
  ok <- grepl("\\(p\\.[^)]+\\)", titulo %||% "")
  ifelse(ok, sub("^.*\\(p\\.([^)]+)\\).*$", "p.\\1", titulo), NA_character_)
}

# Criterios para una tabla anotada a partir de las listas de ClinVar por gen.
.criterios_clinvar <- function(x, listas) {
  n <- nrow(x)
  estrellas <- rep(NA_integer_, n); ps1 <- pm5 <- pm1 <- rep(FALSE, n)
  c_prop <- sub("^[^:]+:", "", .col(x, "hgvs_c", NA_character_))
  p_prop <- .partir_p(.col(x, "hgvs_p", NA_character_))
  missense <- grepl("missense", .col(x, "consecuencia", NA_character_) %||% "")
  patog <- grepl("pathogenic", tolower(.col(x, "clinvar", NA_character_))) &
    !grepl("conflict", tolower(.col(x, "clinvar", NA_character_)))
  for (i in seq_len(n)) {
    g <- x$gen[i]
    if (is.na(g)) next
    L <- listas[[g]]
    if (is.null(L) || nrow(L) == 0) next
    lc <- .c_de_titulo(L$variante)
    lp <- .partir_p(.p_de_titulo(L$variante))
    if (patog[i] && !is.na(c_prop[i])) {
      m <- which(lc == c_prop[i])
      if (length(m) > 0) estrellas[i] <- max(L$estrellas[m])
    }
    if (missense[i] && !is.na(p_prop$pos[i])) {
      mis <- !is.na(lp$pos) & lp$alt != lp$ref & !lp$alt %in% c("Ter", "=")
      ps1[i] <- any(mis & lp$pos == p_prop$pos[i] & lp$alt == p_prop$alt[i] &
                      (is.na(lc) | lc != c_prop[i]))
      pm5[i] <- any(mis & lp$pos == p_prop$pos[i] & lp$alt != p_prop$alt[i])
      cerca <- unique(lp$pos[mis & abs(lp$pos - p_prop$pos[i]) <= 5 & lp$pos != p_prop$pos[i]])
      pm1[i] <- length(cerca) >= 3
    }
  }
  list(estrellas = estrellas, ps1 = ps1, pm5 = pm5, pm1 = pm1)
}

.anotar_clinvar <- function(x, cache) {
  genes <- unique(stats::na.omit(x$gen %||% character(0)))
  interes <- grepl("pathogenic", tolower(.col(x, "clinvar", NA_character_))) |
    grepl("missense", .col(x, "consecuencia", NA_character_) %||% "")
  genes <- intersect(genes, x$gen[interes])
  listas <- list()
  if (length(genes) > 0 && .sin_red()) {
    cli::cli_inform(c("i" = "Sin consultar ClinVar (modo sin red): no hay estrellas ni PS1, PM5 o PM1."))
    genes <- character(0)
  }
  if (length(genes) > 0) {
    .aviso_privacidad("NCBI ClinVar")
    for (g in genes) {
      r <- .en_cache("clinvar_gen", g, cache, function() .clinvar_gen(g))
      if (!is.null(r)) listas[[g]] <- r$variantes
    }
  }
  cr <- .criterios_clinvar(x, listas)
  a <- .attrs_ieiprio(x)
  x$clinvar_estrellas <- cr$estrellas
  x$ps1 <- cr$ps1; x$pm5 <- cr$pm5; x$pm1 <- cr$pm1
  x <- .poner_attrs(x, a)
  attr(x, "ieiprio_clinvar_fecha") <- format(Sys.Date())
  x
}
