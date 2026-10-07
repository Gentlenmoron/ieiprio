.histograma <- function(x, cortes) {
  x <- x[!is.na(x)]
  n <- if (length(x) == 0) rep(0L, length(cortes) - 1) else
    as.integer(table(cut(x, cortes, right = FALSE, include.lowest = TRUE)))
  list(cortes = cortes, n = n)
}

# Estima el sexo con la fraccion de heterocigotos en el X fuera de la PAR.
.sexo_desde_x <- function(d, build, minimo = 20) {
  x <- d[d$chr == "X" & !.en_par(d$pos, build) & !is.na(d$genotipo) &
           d$genotipo != "hom_ref", ]
  n <- nrow(x)
  f <- if (n > 0) mean(x$genotipo == "het") else NA_real_
  sexo <- if (n < minimo) NA_character_ else if (f < 0.1) "M" else if (f > 0.25) "F" else "indeterminado"
  list(sexo = sexo, frac_het = f, n = n)
}

.calcular_qc <- function(v, build) {
  declarado <- attr(v, "ieiprio_sexo")
  muestras <- unique(v$muestra)
  out <- lapply(muestras, function(s) {
    d <- v[v$muestra == s, ]
    llamado <- !is.na(d$genotipo)
    var <- d[llamado & d$genotipo != "hom_ref", ]
    snv <- nchar(var$ref) == 1 & nchar(var$alt) == 1
    ti <- paste0(var$ref, var$alt)[snv] %in% c("AG", "GA", "CT", "TC")
    n_tv <- sum(!ti)
    sx <- .sexo_desde_x(d, build)
    list(
      n_sitios = nrow(d), sin_llamar = sum(!llamado), n_variantes = nrow(var),
      n_het = sum(var$genotipo == "het"),
      n_hom = sum(var$genotipo %in% c("hom_alt", "hemi")),
      titv = if (n_tv > 0) sum(ti) / n_tv else NA_real_,
      mediana_dp = stats::median(var$dp, na.rm = TRUE),
      mediana_gq = stats::median(var$gq, na.rm = TRUE),
      hist_dp = .histograma(pmin(var$dp, 100), c(seq(0, 100, 10), Inf)),
      hist_gq = .histograma(var$gq, seq(0, 100, 10)),
      hist_ab = .histograma(var$ab[var$genotipo == "het"], seq(0, 1, 0.05)),
      sexo_declarado = if (is.null(declarado)) NA_character_ else unname(declarado[s]),
      sexo_inferido = sx$sexo, frac_het_x = sx$frac_het, n_x = sx$n
    )
  })
  stats::setNames(out, muestras)
}

#' Verificar el sexo de las muestras desde el cromosoma X
#'
#' Estima el sexo de cada muestra con la proporcion de variantes
#' heterocigotas en el cromosoma X fuera de las regiones pseudoautosomicas.
#' Un varon casi no tiene heterocigotos ahi. Lo compara con el sexo
#' declarado en [iei_leer_vcf()], porque un error en ese dato cambia la
#' interpretacion de todos los genes ligados al X.
#'
#' @param vcf Tabla devuelta por [iei_leer_vcf()], antes de filtrar.
#' @param minimo Numero minimo de variantes en el X para estimar.
#'
#' @return Un tibble con `muestra`, `sexo_declarado`, `sexo_inferido` (`M`,
#'   `F`, `indeterminado` o `NA` si hay pocas variantes), `frac_het_x`,
#'   `n_variantes_x` y `concordante`.
#' @export
iei_verificar_sexo <- function(vcf, minimo = 20) {
  build <- attr(vcf, "build") %||% unique(vcf$build)[1]
  declarado <- attr(vcf, "ieiprio_sexo")
  muestras <- unique(vcf$muestra)
  filas <- lapply(muestras, function(s) {
    sx <- .sexo_desde_x(vcf[vcf$muestra == s, ], build, minimo)
    dec <- if (is.null(declarado)) NA_character_ else unname(declarado[s])
    tibble::tibble(muestra = s, sexo_declarado = dec, sexo_inferido = sx$sexo,
                   frac_het_x = sx$frac_het, n_variantes_x = sx$n,
                   concordante = if (is.na(dec) || is.na(sx$sexo) || sx$sexo == "indeterminado")
                     NA else dec == sx$sexo)
  })
  out <- do.call(rbind, filas)
  if (any(out$concordante %in% FALSE)) {
    cli::cli_warn("El sexo declarado no coincide con el inferido en: {.val {out$muestra[out$concordante %in% FALSE]}}.")
  }
  out
}
