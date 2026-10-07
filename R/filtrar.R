#' Filtrar variantes por genotipo, panel y calidad
#'
#' Reduce la tabla de [iei_leer_vcf()] a las variantes presentes en genes
#' del panel que pasan los filtros de calidad. Cada paso queda registrado en
#' una bitacora que se consulta con [iei_bitacora()].
#'
#' Los pasos, en orden, son estos.
#' 1. **Genotipo.** Quita `hom_ref` y genotipos sin llamar.
#' 2. **Panel.** Deja las variantes dentro de un gen del panel, con un margen
#'    a cada lado para no perder variantes de splicing. Agrega la columna
#'    `gen`. Una variante en dos genes solapados aparece una vez por gen.
#' 3. **Calidad.** Exige profundidad, calidad de genotipo, FILTER igual a
#'    PASS y una fraccion alelica coherente con el genotipo. Un valor
#'    faltante no descarta la variante.
#'
#' El filtro por frecuencia poblacional se aplica despues de anotar.
#'
#' @param vcf Tabla devuelta por [iei_leer_vcf()].
#' @param genes Genes a conservar. Por defecto, los de [iei_panel()] con el
#'   nivel de `evidencia` indicado.
#' @param evidencia Niveles de evidencia de PanelApp para el panel por
#'   defecto. Se ignora si das `genes`.
#' @param min_dp Profundidad minima.
#' @param min_gq Calidad de genotipo minima.
#' @param min_ab Fraccion alelica minima para heterocigotos. El maximo es
#'   `1 - min_ab`, y homocigotos y hemicigotos deben tener al menos
#'   `1 - min_ab`.
#' @param solo_pass Si `TRUE`, solo variantes con FILTER igual a `PASS` o `.`.
#' @param margen Pares de bases que se agregan a cada lado del gen.
#' @param coordenadas Tabla de coordenadas de genes. Por defecto
#'   [genes_coordenadas].
#'
#' @return La tabla filtrada con la columna `gen` agregada despues de
#'   `muestra`, y el atributo `ieiprio_bitacora`.
#' @export
#' @examples
#' \dontrun{
#' vcf <- iei_leer_vcf("paciente01.vcf.gz", sexo = "M")
#' filt <- iei_filtrar(vcf)
#' iei_bitacora(filt)
#' }
iei_filtrar <- function(vcf, genes = NULL, evidencia = "verde",
                        min_dp = 10, min_gq = 20, min_ab = 0.2,
                        solo_pass = TRUE, margen = 20,
                        coordenadas = ieiprio::genes_coordenadas) {
  necesarias <- c("muestra", "chr", "pos", "genotipo", "dp", "gq", "ab", "filtro")
  if (!is.data.frame(vcf) || !all(necesarias %in% names(vcf))) {
    cli::cli_abort("{.arg vcf} debe ser la tabla que devuelve {.fn iei_leer_vcf}.")
  }
  if (min_ab < 0 || min_ab >= 0.5) {
    cli::cli_abort("{.arg min_ab} debe estar entre 0 y 0.5.")
  }
  build <- attr(vcf, "build") %||% unique(vcf$build)[1]
  if (is.null(build) || is.na(build)) {
    cli::cli_abort("No se sabe el build del VCF. Vuelve a leerlo con {.fn iei_leer_vcf}.")
  }
  genes_por_defecto <- is.null(genes)
  if (genes_por_defecto) genes <- iei_panel(evidencia = evidencia)$gen
  a <- .attrs_ieiprio(vcf)
  qc <- .calcular_qc(vcf, build)

  bitacora <- list(c("Inicial", nrow(vcf)))
  x <- vcf

  # 1. Genotipo
  x <- x[!is.na(x$genotipo) & x$genotipo != "hom_ref", ]
  bitacora[[length(bitacora) + 1]] <- c("Genotipo", nrow(x))

  # 2. Panel
  co <- coordenadas[coordenadas$build == build & coordenadas$gen %in% genes, ]
  sin_coord <- setdiff(genes, co$gen)
  if (length(sin_coord) > 0) {
    cli::cli_inform(c("i" = "{length(sin_coord)} gen{?es} sin coordenadas en {build} no se pueden buscar."))
  }
  x <- .asignar_genes(x, co, margen)
  bitacora[[length(bitacora) + 1]] <- c("Panel", nrow(x))

  # 3. Calidad
  ok_dp <- is.na(x$dp) | x$dp >= min_dp
  ok_gq <- is.na(x$gq) | x$gq >= min_gq
  ok_filtro <- !solo_pass | x$filtro %in% c("PASS", ".")
  ok_ab <- is.na(x$ab) |
    (x$genotipo == "het" & x$ab >= min_ab & x$ab <= 1 - min_ab) |
    (x$genotipo %in% c("hom_alt", "hemi") & x$ab >= 1 - min_ab)
  x <- x[ok_dp & ok_gq & ok_filtro & ok_ab, ]
  bitacora[[length(bitacora) + 1]] <- c("Calidad", nrow(x))

  b <- do.call(rbind, bitacora)
  n <- as.integer(b[, 2])
  x <- .poner_attrs(x, a)
  attr(x, "ieiprio_bitacora") <- tibble::tibble(
    paso = b[, 1], filas = n, descartadas = c(0L, -diff(n))
  )
  attr(x, "build") <- build
  attr(x, "ieiprio_qc") <- qc
  attr(x, "ieiprio_parametros") <- list(
    genes = if (genes_por_defecto) sprintf("panel PanelApp, evidencia %s (%d genes)",
                                           paste(evidencia, collapse = "/"), length(genes))
            else sprintf("lista propia (%d genes)", length(genes)),
    min_dp = min_dp, min_gq = min_gq, min_ab = min_ab, solo_pass = solo_pass, margen = margen
  )
  x
}

# Cruza variantes con intervalos de genes. Devuelve una fila por variante y gen.
.asignar_genes <- function(x, co, margen) {
  if (nrow(x) == 0 || nrow(co) == 0) {
    out <- x[0, ]
    out <- tibble::add_column(out, gen = character(0), .after = "muestra")
    return(out)
  }
  pares <- lapply(seq_len(nrow(co)), function(j) {
    idx <- which(x$chr == co$chr[j] &
                 x$pos >= co$inicio[j] - margen &
                 x$pos <= co$fin[j] + margen)
    if (length(idx) == 0) NULL else cbind(idx, j)
  })
  pares <- do.call(rbind, pares)
  if (is.null(pares)) return(.asignar_genes(x[0, ], co[0, ], margen))
  pares <- pares[order(pares[, 1], co$gen[pares[, 2]]), , drop = FALSE]
  out <- x[pares[, 1], ]
  tibble::add_column(out, gen = co$gen[pares[, 2]], .after = "muestra")
}

#' Ver la bitacora del filtrado
#'
#' @param x Tabla devuelta por [iei_filtrar()].
#' @return Un tibble con cada paso, las filas que quedaron y las descartadas.
#'   Cada fila es una variante en una muestra.
#' @export
iei_bitacora <- function(x) {
  b <- attr(x, "ieiprio_bitacora")
  if (is.null(b)) cli::cli_abort("{.arg x} no tiene bitacora. Pasalo antes por {.fn iei_filtrar}.")
  b
}
