# Funciones internas para armar el VCF de ejemplo desde ClinVar (NCBI E-utilities).
# Las usa data-raw/ejemplo_vcf.R; no forman parte del flujo de analisis.

.eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

.clinvar_peticion <- function(endpoint, ...) {
  httr2::request(.eutils) |>
    httr2::req_url_path_append(endpoint) |>
    httr2::req_url_query(db = "clinvar", retmode = "json", ...) |>
    .agente() |>
    httr2::req_timeout(60) |>
    httr2::req_throttle(rate = 2 / 1) |>
    httr2::req_retry(max_tries = 4) |>
    httr2::req_perform() |>
    httr2::resp_body_json(simplifyVector = FALSE)
}

.clinvar_buscar <- function(termino, retmax = 300) {
  r <- .clinvar_peticion("esearch.fcgi", term = termino, retmax = retmax)
  ids <- unlist(r$esearchresult$idlist)
  if (is.null(ids)) character(0) else sort(as.integer(ids))
}

.clinvar_resumen <- function(ids) {
  if (length(ids) == 0) return(list())
  grupos <- split(ids, ceiling(seq_along(ids) / 150))
  unlist(unname(lapply(grupos, function(g) {
    r <- .clinvar_peticion("esummary.fcgi", id = paste(g, collapse = ","))$result
    r[setdiff(names(r), "uids")]
  })), recursive = FALSE)
}

# Convierte SPDI (posicion 0-based) a coordenadas VCF. Solo variantes de un nucleotido.
.spdi_a_vcf <- function(spdi) {
  p <- strsplit(spdi, ":", fixed = TRUE)
  out <- lapply(p, function(x) {
    if (length(x) != 4 || nchar(x[3]) != 1 || nchar(x[4]) != 1 ||
        !grepl("^NC_0+[0-9]+\\.", x[1])) return(NULL)
    n <- as.integer(sub("^NC_0+([0-9]+)\\..*", "\\1", x[1]))
    chr <- if (n == 23) "X" else if (n == 24) "Y" else if (n == 12920) "MT" else as.character(n)
    data.frame(chr = chr, pos = as.integer(x[2]) + 1L, ref = x[3], alt = x[4])
  })
  out
}

# Tabla con una fila por variante SNV de un solo gen.
.clinvar_procesar <- function(resumen) {
  filas <- lapply(seq_along(resumen), function(i) {
    r <- resumen[[i]]
    id <- .chr1(r$uid)
    if (is.na(id)) id <- names(resumen)[i]
    vs <- r$variation_set
    if (length(vs) != 1) return(NULL)
    spdi <- .chr1(vs[[1]]$canonical_spdi)
    if (is.na(spdi)) return(NULL)
    coord <- .spdi_a_vcf(spdi)[[1]]
    if (is.null(coord)) return(NULL)
    genes <- vapply(r$genes %||% list(), function(g) .chr1(g$symbol), character(1))
    if (length(genes) != 1) return(NULL)
    clas <- .chr1(r$germline_classification$description %||% r$clinical_significance$description)
    data.frame(clinvar_id = as.integer(id), gen = genes, clasificacion = clas, coord,
               nombre = .chr1(r$title), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, filas)
  if (is.null(out)) {
    out <- data.frame(clinvar_id = integer(0), gen = character(0), clasificacion = character(0),
                      chr = character(0), pos = integer(0), ref = character(0),
                      alt = character(0), nombre = character(0))
  }
  tibble::as_tibble(out[order(out$clinvar_id), ])
}

# Elige la primera variante (menor ID) de un gen con la clasificacion pedida.
.clinvar_elegir <- function(gen, clasificacion, excluir = integer(0)) {
  propiedad <- c(Pathogenic = "clinsig_pathogenic", Benign = "clinsig_benign",
                 `Uncertain significance` = "clinsig_vus")[[clasificacion]]
  ids <- .clinvar_buscar(sprintf("%s[gene] AND %s[Properties]", gen, propiedad))
  if (length(ids) == 0) ids <- .clinvar_buscar(sprintf("%s[gene]", gen), retmax = 500)
  ids <- setdiff(ids, excluir)
  t <- .clinvar_procesar(.clinvar_resumen(utils::head(ids, 150)))
  t <- t[t$gen == gen & t$clasificacion %in% clasificacion, ]
  if (nrow(t) == 0) cli::cli_abort("No encontre una variante {clasificacion} de un nucleotido en {gen}.")
  t[1, ]
}

#' Ruta del VCF de ejemplo
#'
#' El paquete incluye un VCF de un paciente ficticio (EJEMPLO01, varon) con
#' cinco variantes reales tomadas de ClinVar en genes del panel y unas 4000
#' variantes de fondo sinteticas fuera de esos genes. Se genera con
#' `data-raw/ejemplo_vcf.R`; las variantes de ClinVar estan fijadas por su
#' identificador en `data-raw/ejemplo_clinvar.csv`.
#'
#' @return La ruta al archivo `ejemplo.vcf.gz`.
#' @export
#' @examples
#' \dontrun{
#' v <- iei_leer_vcf(iei_ejemplo_vcf(), sexo = "M")
#' }
iei_ejemplo_vcf <- function() {
  ruta <- system.file("extdata", "ejemplo.vcf.gz", package = "ieiprio")
  if (!nzchar(ruta)) {
    cli::cli_abort("El VCF de ejemplo no esta instalado. Generalo con {.file data-raw/ejemplo_vcf.R}.")
  }
  ruta
}
