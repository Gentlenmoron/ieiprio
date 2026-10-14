# Consultas a la API de gnomAD (GraphQL): homocigotos, hemicigotos, frecuencia
# filtrada al 95 % y restriccion del gen. VEP no da estos datos.

.gnomad_api <- "https://gnomad.broadinstitute.org/api"

.gnomad_consulta <- function(query, variables) {
  .comprobar_red("gnomAD")
  httr2::request(.gnomad_api) |>
    httr2::req_body_json(list(query = query, variables = variables)) |>
    .agente() |>
    httr2::req_timeout(60) |>
    httr2::req_throttle(rate = 10 / 60) |>
    httr2::req_retry(max_tries = 4, is_transient = .transitorio) |>
    .ejecutar_json("gnomAD")
}

.gnomad_dataset <- function(build) if (identical(build, "GRCh37")) "gnomad_r2_1" else "gnomad_r4"

.gnomad_q_variante <- "
query Variante($id: String!, $dataset: DatasetId!) {
  variant(variantId: $id, dataset: $dataset) {
    exome { ac an homozygote_count hemizygote_count faf95 { popmax popmax_population } }
    genome { ac an homozygote_count hemizygote_count faf95 { popmax popmax_population } }
  }
}"

.gnomad_q_gen <- "
query Gen($gen: String!, $ref: ReferenceGenomeId!) {
  gene(gene_symbol: $gen, reference_genome: $ref) {
    gnomad_constraint { oe_lof_upper mis_z pli }
  }
}"

# Junta exomas y genomas: suma conteos y toma la mayor FAF95.
.gnomad_variante_procesar <- function(r) {
  v <- r$data$variant
  if (is.null(v)) return(list(hom = 0L, hemi = 0L, faf95 = NA_real_, grupo_faf = NA_character_, encontrada = FALSE))
  partes <- Filter(Negate(is.null), list(v$exome, v$genome))
  num <- function(campo) sum(vapply(partes, function(p) as.numeric(p[[campo]] %||% 0), numeric(1)))
  faf <- vapply(partes, function(p) as.numeric(p$faf95$popmax %||% NA), numeric(1))
  grupo <- vapply(partes, function(p) as.character(p$faf95$popmax_population %||% NA), character(1))
  k <- if (all(is.na(faf))) NA_integer_ else which.max(faf)
  list(hom = as.integer(num("homozygote_count")), hemi = as.integer(num("hemizygote_count")),
       faf95 = if (is.na(k)) NA_real_ else faf[k], grupo_faf = if (is.na(k)) NA_character_ else grupo[k],
       encontrada = TRUE)
}

#' Consultar gnomAD para las variantes candidatas
#'
#' Agrega datos de gnomAD que Ensembl VEP no entrega: cuantas personas sanas
#' son homocigotas o hemicigotas para la variante (la prueba mas fuerte contra
#' una variante en un gen recesivo o ligado al X, criterio BS2), la frecuencia
#' filtrada al 95 % que recomienda ClinGen para BA1 y BS1, y la restriccion
#' del gen (intolerancia a perdida de funcion y a missense, criterio PP2).
#'
#' gnomAD limita las consultas a unas diez por minuto, asi que conviene
#' aplicarlo solo a las candidatas: despues de [iei_filtrar_frecuencia()] o a
#' las primeras filas de un resultado, y luego volver a priorizar.
#'
#' @param x Tabla anotada con [iei_anotar()].
#' @param max_variantes Maximo de variantes a consultar; las demas quedan en NA.
#' @param cache Carpeta de la cache, o `FALSE`.
#' @return `x` con `gnomad_hom`, `gnomad_hemi`, `faf95`, `grupo_faf95`,
#'   `loeuf` (limite superior de observado/esperado para perdida de funcion;
#'   menor a 0.6 indica intolerancia) y `mis_z` (mayor a 3.09 indica
#'   intolerancia a missense).
#' @export
#' @examples
#' \dontrun{
#' anot <- iei_gnomad(iei_filtrar_frecuencia(iei_anotar(filtrado)))
#' }
iei_gnomad <- function(x, max_variantes = 50, cache = tools::R_user_dir("ieiprio", "cache")) {
  build <- attr(x, "build") %||% "GRCh38"
  a <- .attrs_ieiprio(x)
  n <- nrow(x)
  hom <- hemi <- rep(NA_integer_, n); faf <- loeuf <- misz <- rep(NA_real_, n); gfaf <- rep(NA_character_, n)
  if (n > 0) {
    .aviso_privacidad("gnomAD")
    ids <- paste(x$chr, x$pos, x$ref, x$alt, sep = "-")
    unicos <- unique(ids)
    if (length(unicos) > max_variantes) {
      cli::cli_warn("Se consultan solo las primeras {max_variantes} de {length(unicos)} variantes (gnomAD limita las consultas).")
      unicos <- unicos[seq_len(max_variantes)]
    }
    cli::cli_inform("Consultando gnomAD para {length(unicos)} variante{?s}; puede tardar unos minutos.")
    res <- lapply(unicos, function(id) {
      .en_cache(paste0("gnomad_", build), gsub("[^A-Za-z0-9]", "-", id), cache, function() {
        .gnomad_variante_procesar(.gnomad_consulta(.gnomad_q_variante,
                                                   list(id = id, dataset = .gnomad_dataset(build))))
      })
    })
    names(res) <- unicos
    for (i in seq_len(n)) {
      r <- res[[ids[i]]]
      if (is.null(r)) next
      hom[i] <- r$hom; hemi[i] <- r$hemi; faf[i] <- r$faf95; gfaf[i] <- r$grupo_faf
    }
    genes <- unique(stats::na.omit(x$gen %||% character(0)))
    rg <- lapply(genes, function(g) {
      .en_cache(paste0("gnomad_gen_", build), g, cache, function() {
        r <- .gnomad_consulta(.gnomad_q_gen, list(gen = g, ref = build))
        cst <- r$data$gene$gnomad_constraint
        list(loeuf = as.numeric(cst$oe_lof_upper %||% NA), mis_z = as.numeric(cst$mis_z %||% NA))
      })
    })
    names(rg) <- genes
    for (i in seq_len(n)) {
      if (is.na(x$gen[i])) next
      r <- rg[[x$gen[i]]]
      if (is.null(r)) next
      loeuf[i] <- r$loeuf; misz[i] <- r$mis_z
    }
  }
  x$gnomad_hom <- hom; x$gnomad_hemi <- hemi; x$faf95 <- faf; x$grupo_faf95 <- gfaf
  x$loeuf <- loeuf; x$mis_z <- misz
  .poner_attrs(x, a)
}
