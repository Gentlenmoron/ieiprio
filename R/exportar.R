# Exportacion de resultados a otros formatos y programas.

#' Exportar resultados
#'
#' Guarda una tabla de ieiprio (de cualquier etapa: lectura, anotacion o
#' priorizacion) en el formato que necesite el siguiente paso.
#'
#' * `"tsv"`: texto separado por tabuladores, para Excel, R o Python.
#' * `"xlsx"`: libro de Excel con una hoja de variantes, una de parametros y
#'   una de fuentes de datos (necesita el paquete writexl).
#' * `"vcf"`: VCF 4.2 con una columna por muestra y la anotacion de ieiprio
#'   en INFO (gen, consecuencia, HGVS, puntaje, prioridad y ACMG), para abrir
#'   en IGV o pasar a otro programa.
#' * `"bed"`: coordenadas en formato BED (base 0), para IGV o bedtools.
#' * `"phenopacket"`: un Phenopacket de GA4GH (version 2) en JSON con el
#'   fenotipo del caso, los signos ausentes y las variantes de prioridad alta
#'   y media como interpretaciones genomicas. Es el formato que leen
#'   Exomiser, LIRICAL y los registros de enfermedades raras.
#'
#' @param x Tabla de ieiprio.
#' @param archivo Ruta del archivo a escribir.
#' @param formato `"tsv"`, `"xlsx"`, `"vcf"`, `"bed"` o `"phenopacket"`. Por
#'   defecto se deduce de la extension (`.json` es phenopacket).
#' @param muestra Para `"phenopacket"`, la muestra a exportar si hay mas de
#'   una.
#' @return La ruta del archivo, invisible.
#' @export
#' @examples
#' \dontrun{
#' iei_exportar(prior, "resultado.xlsx")
#' iei_exportar(prior, "candidatas.vcf")
#' iei_exportar(prior, "P01.json", muestra = "P01")
#' }
iei_exportar <- function(x, archivo, formato = NULL, muestra = NULL) {
  if (is.null(formato)) {
    ext <- tolower(tools::file_ext(archivo))
    formato <- switch(ext, txt = , tsv = "tsv", xlsx = "xlsx", vcf = "vcf", bed = "bed",
                      json = "phenopacket",
                      cli::cli_abort("No reconozco la extension {.val {ext}}; indica {.arg formato}."))
  }
  formato <- match.arg(formato, c("tsv", "xlsx", "vcf", "bed", "phenopacket"))
  switch(formato,
    tsv = utils::write.table(.tabla_plana(x), archivo, sep = "\t", row.names = FALSE, quote = FALSE, na = ""),
    xlsx = .exportar_xlsx(x, archivo),
    vcf = writeLines(.vcf_lineas(x), archivo),
    bed = .exportar_bed(x, archivo),
    phenopacket = writeLines(jsonlite::toJSON(iei_phenopacket(x, muestra), auto_unbox = TRUE,
                                              pretty = TRUE, null = "null"), archivo)
  )
  cli::cli_inform("Guardado en {.path {archivo}} ({formato}).")
  invisible(archivo)
}

# Columnas de lista o con saltos de linea pasan a texto plano.
.tabla_plana <- function(x) {
  d <- as.data.frame(tibble::as_tibble(unclass(x)[names(x)]), stringsAsFactors = FALSE)
  for (n in names(d)) {
    if (is.list(d[[n]])) d[[n]] <- vapply(d[[n]], function(v) paste(v, collapse = ";"), "")
    if (is.character(d[[n]])) d[[n]] <- gsub("[\t\r\n]+", " ", d[[n]])
  }
  d
}

.exportar_xlsx <- function(x, archivo) {
  rlang::check_installed("writexl", reason = "para exportar a Excel.")
  a <- iei_metadatos(x)
  param <- data.frame(
    parametro = c("ieiprio", "fecha", "genoma", "Ensembl VEP", "fenotipo", "fenotipo ausente", "similitud"),
    valor = c(.version_paquete(), format(Sys.Date()), a$build %||% NA, a$vep_release %||% NA,
              paste(a$fenotipo %||% character(0), collapse = "; "),
              paste(a$fenotipo_ausente %||% character(0), collapse = "; "),
              a$similitud %||% NA),
    stringsAsFactors = FALSE)
  hojas <- list(variantes = .tabla_plana(x), parametros = param,
                fuentes = as.data.frame(ieiprio::fuentes_datos))
  writexl::write_xlsx(hojas, archivo)
}

.gt_vcf <- function(g) {
  out <- c(het = "0/1", hom_alt = "1/1", hemi = "1", hom_ref = "0/0")[g]
  ifelse(is.na(out), "./.", out)
}

.info_valor <- function(v) {
  v <- gsub("[;=, \t]", "_", as.character(v))
  v <- gsub("_+", "_", v)
  ifelse(is.na(v) | v == "" | v == "NA", NA_character_, v)
}

.vcf_lineas <- function(x) {
  if (!all(c("chr", "pos", "ref", "alt") %in% names(x))) cli::cli_abort("{.arg x} necesita las columnas chr, pos, ref y alt.")
  muestras <- unique(x$muestra %||% "MUESTRA")
  muestra <- x$muestra %||% rep("MUESTRA", nrow(x))
  clave <- paste(x$chr, x$pos, x$ref, x$alt, sep = "\t")
  sitios <- !duplicated(clave)
  campos <- c(GENE = "gen", CSQ = "consecuencia", HGVSC = "hgvs_c", HGVSP = "hgvs_p",
              CLINVAR = "clinvar", AF = "af_gnomad", SCORE = "puntaje", PRIORITY = "categoria",
              ACMG = "acmg", ACMG_CLASS = "acmg_clase")
  campos <- campos[campos %in% names(x)]
  desc <- c(GENE = "Gen", CSQ = "Consecuencia (Ensembl VEP)", HGVSC = "HGVS c.", HGVSP = "HGVS p.",
            CLINVAR = "Significancia en ClinVar", AF = "Frecuencia en gnomAD", SCORE = "Puntaje de ieiprio",
            PRIORITY = "Prioridad de ieiprio", ACMG = "Criterios ACMG sugeridos",
            ACMG_CLASS = "Clase ACMG por puntos")
  tipo <- c(AF = "Float", SCORE = "Float")
  build <- attr(x, "build") %||% (x$build[1] %||% "GRCh38")
  cab <- c("##fileformat=VCFv4.2", sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
           sprintf("##source=ieiprio_%s", .version_paquete()), sprintf("##reference=%s", build),
           sprintf('##INFO=<ID=%s,Number=1,Type=%s,Description="%s">', names(campos),
                   ifelse(names(campos) %in% names(tipo), tipo[names(campos)], "String"), desc[names(campos)]),
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotipo">',
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", muestras), collapse = "\t"))
  s <- which(sitios)
  info <- vapply(s, function(i) {
    v <- vapply(campos, function(cc) .info_valor(x[[cc]][i]), "")
    v <- v[!is.na(v)]
    if (length(v) == 0) "." else paste(sprintf("%s=%s", names(v), v), collapse = ";")
  }, "")
  gts <- vapply(muestras, function(m) {
    g <- rep(NA_character_, length(s))
    k <- match(clave[s], clave[muestra == m])
    g[!is.na(k)] <- x$genotipo[muestra == m][k[!is.na(k)]]
    .gt_vcf(g)
  }, character(length(s)))
  gts <- matrix(gts, nrow = length(s))
  base <- sub("^chr", "", x$chr[s])
  chr <- if (grepl("38", build)) paste0("chr", ifelse(base == "MT", "M", base)) else base
  rango <- suppressWarnings(as.numeric(base))
  rango[base == "X"] <- 23; rango[base == "Y"] <- 24; rango[base %in% c("M", "MT")] <- 25
  orden <- order(rango, x$pos[s])
  id <- .col(x, "id", NA_character_)[s]
  filas <- paste(chr, x$pos[s], ifelse(is.na(id), ".", id), x$ref[s], x$alt[s],
                 ".", "PASS", info, "GT", apply(gts, 1, paste, collapse = "\t"), sep = "\t")
  c(cab, filas[orden])
}

.exportar_bed <- function(x, archivo) {
  nombre <- paste0(x$gen %||% "", ifelse(is.na(.col(x, "hgvs_p", NA_character_)), "",
                                          paste0(":", sub("^[^:]+:", "", .col(x, "hgvs_p", NA_character_)))))
  bed <- data.frame(chrom = x$chr, inicio = as.integer(x$pos) - 1L,
                    fin = as.integer(x$pos) - 1L + pmax(nchar(x$ref), 1L),
                    nombre = gsub("[\t ]", "_", nombre),
                    puntaje = round(pmin(1000, pmax(0, 100 * .col(x, "puntaje", 0))), 0))
  bed$puntaje[is.na(bed$puntaje)] <- 0
  utils::write.table(bed, archivo, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

#' Convertir a GRanges
#'
#' Devuelve la tabla como un objeto `GRanges` de Bioconductor, con todas las
#' columnas como metadatos, para seguir el analisis con GenomicRanges,
#' VariantAnnotation o plyranges.
#'
#' @param x Tabla de ieiprio con `chr`, `pos` y `ref`.
#' @return Un objeto `GRanges` de GenomicRanges, que se instala aparte con
#'   `BiocManager::install("GenomicRanges")`.
#' @export
#' @examples
#' \dontrun{
#' gr <- iei_granges(prior)
#' }
iei_granges <- function(x) {
  rlang::check_installed("GenomicRanges", reason = "para convertir a GRanges (instalalo con BiocManager).")
  d <- .tabla_plana(x)
  d$fin <- as.integer(d$pos) + pmax(nchar(d$ref), 1L) - 1L
  # Se llama sin declarar GenomicRanges como dependencia, para no obligar a
  # instalar Bioconductor a quien no lo usa.
  hacer <- getExportedValue("GenomicRanges", "makeGRangesFromDataFrame")
  gr <- hacer(d, seqnames.field = "chr", start.field = "pos", end.field = "fin", keep.extra.columns = TRUE)
  gr
}

# Enfermedad mas parecida al fenotipo entre las del gen (o la primera del
# gen si no hay fenotipo), con su identificador OMIM u Orphanet.
.enfermedad_candidata <- function(gen, fenotipo) {
  if (is.null(gen) || is.na(gen)) return(NULL)
  enf <- ieiprio::hpo_enfermedades[ieiprio::hpo_enfermedades$gen == gen, ]
  if (nrow(enf) == 0) return(NULL)
  if (length(fenotipo) > 0) {
    r <- .fenotipo_por_gen(fenotipo, gen, "enfermedad")
    if (!is.na(r$enfermedad_id)) return(list(id = r$enfermedad_id, label = trimws(r$enfermedad)))
  }
  list(id = enf$enfermedad_id[1], label = trimws(enf$enfermedad[1]))
}

#' Phenopacket de un caso
#'
#' Arma un Phenopacket de GA4GH (esquema 2.0) como lista de R: el sujeto, el
#' fenotipo usado en la priorizacion (con los signos ausentes marcados como
#' `excluded`) y las variantes de prioridad alta y media como
#' interpretaciones genomicas candidatas, con su clase ACMG sugerida. El
#' diagnostico propuesto es la enfermedad del gen de la primera variante que
#' mas se parece al fenotipo. Se guarda como JSON con [iei_exportar()].
#'
#' @param x Resultado de [iei_priorizar()].
#' @param muestra Muestra a exportar; obligatoria si hay mas de una.
#' @return Una lista con la estructura del Phenopacket.
#' @export
#' @examples
#' \dontrun{
#' pp <- iei_phenopacket(prior, muestra = "P01")
#' }
iei_phenopacket <- function(x, muestra = NULL) {
  muestras <- unique(x$muestra)
  if (is.null(muestra)) {
    if (length(muestras) > 1) cli::cli_abort("Hay {length(muestras)} muestras; indica {.arg muestra}.")
    muestra <- muestras[1]
  }
  if (!muestra %in% muestras) cli::cli_abort("La muestra {.val {muestra}} no esta en {.arg x}.")
  a <- iei_metadatos(x)
  d <- x[x$muestra == muestra, ]
  if ("categoria" %in% names(d)) d <- d[d$categoria %in% c("alta", "media"), ]
  sexo <- a$sexo[[muestra]] %||% NA
  sexo <- if (identical(sexo, "M")) "MALE" else if (identical(sexo, "F")) "FEMALE" else "UNKNOWN_SEX"
  rasgo <- function(id, excluido) {
    r <- list(type = list(id = id, label = .hpo_nombre(id) %|% id))
    if (excluido) r$excluded <- TRUE
    r
  }
  rasgos <- c(lapply(a$fenotipo %||% character(0), rasgo, excluido = FALSE),
              lapply(a$fenotipo_ausente %||% character(0), rasgo, excluido = TRUE))
  acmg_enum <- c("Patog\u00e9nica" = "PATHOGENIC", "Probablemente patog\u00e9nica" = "LIKELY_PATHOGENIC",
                 "Significado incierto" = "UNCERTAIN_SIGNIFICANCE",
                 "Probablemente benigna" = "LIKELY_BENIGN", "Benigna" = "BENIGN")
  cigosidad <- list(het = list(id = "GENO:0000135", label = "heterozygous"),
                    hom_alt = list(id = "GENO:0000136", label = "homozygous"),
                    hemi = list(id = "GENO:0000134", label = "hemizygous"))
  build <- a$build %||% "GRCh38"
  interp <- lapply(seq_len(nrow(d)), function(i) {
    r <- d[i, ]
    expr <- Filter(Negate(is.null), list(
      if (!is.na(.celda(r, "hgvs_c"))) list(syntax = "hgvs.c", value = .celda(r, "hgvs_c")),
      if (!is.na(.celda(r, "hgvs_p"))) list(syntax = "hgvs.p", value = .celda(r, "hgvs_p"))))
    vd <- list(id = sprintf("%s-%s-%s-%s", r$chr, r$pos, r$ref, r$alt),
               geneContext = list(valueId = .celda(r, "gen"), symbol = .celda(r, "gen")),
               expressions = expr,
               vcfRecord = list(genomeAssembly = build, chrom = as.character(r$chr), pos = as.character(r$pos),
                                ref = r$ref, alt = r$alt))
    z <- cigosidad[[.celda(r, "genotipo") %|% ""]]
    if (!is.null(z)) vd$allelicState <- z
    clase <- acmg_enum[.celda(r, "acmg_clase") %|% ""]
    list(subjectOrBiosampleId = muestra, interpretationStatus = "CANDIDATE",
         variantInterpretation = list(
           acmgPathogenicityClassification = if (is.na(clase)) "NOT_PROVIDED" else unname(clase),
           variationDescriptor = vd))
  })
  pp <- list(
    id = paste0("ieiprio-", muestra),
    subject = list(id = muestra, sex = sexo),
    phenotypicFeatures = rasgos,
    interpretations = if (length(interp) == 0) list() else list(list(
      id = paste0("ieiprio-", muestra, "-1"), progressStatus = "IN_PROGRESS",
      diagnosis = Filter(Negate(is.null), list(disease = .enfermedad_candidata(d$gen[1], a$fenotipo),
                                               genomicInterpretations = interp)))),
    metaData = list(
      created = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      createdBy = paste("ieiprio", .version_paquete()),
      resources = list(list(id = "hp", name = "human phenotype ontology", url = "http://purl.obolibrary.org/obo/hp.owl",
                            version = ieiprio::fuentes_datos$version[ieiprio::fuentes_datos$fuente == "HPO ontologia"][1],
                            namespacePrefix = "HP", iriPrefix = "http://purl.obolibrary.org/obo/HP_"),
                       list(id = "geno", name = "Genotype Ontology", url = "http://purl.obolibrary.org/obo/geno.owl",
                            version = "2023-10-08", namespacePrefix = "GENO", iriPrefix = "http://purl.obolibrary.org/obo/GENO_")),
      phenopacketSchemaVersion = "2.0")
  )
  pp
}
