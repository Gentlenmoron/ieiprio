# CNV y variantes estructurales leidas del VCF.

.info_campo <- function(info, campo) {
  patron <- sprintf("^(.*;)?%s=([^;]+).*$", campo)
  ifelse(grepl(patron, info), sub(patron, "\\2", info), NA_character_)
}

.leer_cnv <- function(m, muestras, build) {
  tipo <- .info_campo(m[, "INFO"], "SVTYPE")
  alt_tipo <- sub("^<([A-Z]+).*$", "\\1", m[, "ALT"])
  tipo <- ifelse(is.na(tipo), ifelse(grepl("^<", m[, "ALT"]), alt_tipo, NA_character_), tipo)
  inicio <- as.integer(m[, "POS"])
  fin <- suppressWarnings(as.integer(.info_campo(m[, "INFO"], "END")))
  largo <- suppressWarnings(abs(as.integer(.info_campo(m[, "INFO"], "SVLEN"))))
  fin <- ifelse(is.na(fin), inicio + ifelse(is.na(largo), 0L, largo), fin)
  chr <- .normalizar_chr(m[, "CHROM"])
  formato <- m[, "FORMAT"]
  tablas <- lapply(muestras, function(s) {
    gt <- .extraer_campo(formato, m[, s], "GT")
    cn <- suppressWarnings(as.numeric(.extraer_campo(formato, m[, s], "CN")))
    presente <- !is.na(gt) & grepl("[1-9]", gt) | (!is.na(cn) & cn != 2)
    tibble::tibble(muestra = s, chr = chr, inicio = inicio, fin = as.integer(fin), tipo = tipo,
                   gt = gt, cn = cn, filtro = m[, "FILTER"], build = build)[presente, ]
  })
  do.call(rbind, tablas)
}

#' CNV y variantes estructurales en genes del panel
#'
#' Devuelve los CNV y variantes estructurales que traia el VCF (registros con
#' ALT simbolico como `<DEL>` o `<DUP>`, o con `SVTYPE` en INFO), con los
#' genes del panel que tocan. ieiprio no llama CNV; esta funcion muestra los
#' que ya llamo otro programa (por ejemplo un llamador de CNV de exoma), para
#' no perder las deleciones grandes de genes como BTK, CYBB o la region
#' 22q11.2.
#'
#' @param x Tabla de [iei_leer_vcf()] o cualquier tabla derivada que conserve
#'   sus atributos.
#' @param genes Genes a cruzar. Por defecto, todo el panel.
#' @return Un tibble con `muestra`, `chr`, `inicio`, `fin`, `tamano` (pb),
#'   `tipo`, `cn` (numero de copias si el VCF lo trae), `genes` del panel que
#'   toca y `nota`.
#' @export
#' @examples
#' \dontrun{
#' iei_cnv(iei_leer_vcf("paciente_con_cnv.vcf.gz"))
#' }
iei_cnv <- function(x, genes = NULL) {
  cnv <- attr(x, "ieiprio_cnv")
  vacio <- tibble::tibble(muestra = character(0), chr = character(0), inicio = integer(0),
                          fin = integer(0), tamano = integer(0), tipo = character(0),
                          cn = numeric(0), genes = character(0), nota = character(0))
  if (is.null(cnv) || nrow(cnv) == 0) {
    cli::cli_inform("El VCF no trae CNV ni variantes estructurales.")
    return(vacio)
  }
  build <- unique(cnv$build)[1]
  co <- ieiprio::genes_coordenadas[ieiprio::genes_coordenadas$build == build, ]
  if (!is.null(genes)) co <- co[co$gen %in% genes, ]
  toca <- vapply(seq_len(nrow(cnv)), function(i) {
    g <- co$gen[co$chr == cnv$chr[i] & co$inicio <= cnv$fin[i] & co$fin >= cnv$inicio[i]]
    if (length(g) == 0) NA_character_ else paste(sort(unique(g)), collapse = ";")
  }, character(1))
  out <- tibble::tibble(muestra = cnv$muestra, chr = cnv$chr, inicio = cnv$inicio, fin = cnv$fin,
                        tamano = as.integer(cnv$fin - cnv$inicio + 1L), tipo = cnv$tipo, cn = cnv$cn,
                        genes = toca)
  dif <- ieiprio::genes_dificiles
  out$nota <- vapply(out$genes, function(g) {
    if (is.na(g)) return(NA_character_)
    gs <- strsplit(g, ";", fixed = TRUE)[[1]]
    if ("TBX1" %in% gs) return("Compatible con la regi\u00f3n 22q11.2; confirmar con MLPA o microarreglo.")
    if (any(gs %in% dif$gen[dif$tipo %in% c("numero de copias", "duplicacion segmentaria")]))
      return("Gen de n\u00famero de copias variable; confirmar con MLPA o un m\u00e9todo dirigido.")
    "Confirmar con un m\u00e9todo independiente (MLPA, microarreglo o PCR cuantitativa)."
  }, character(1))
  out <- out[order(is.na(out$genes), out$muestra, out$chr, out$inicio), ]
  out
}
