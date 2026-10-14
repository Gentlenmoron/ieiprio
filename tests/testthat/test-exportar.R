priorizado_exp <- function() {
  x <- tibble::tibble(
    muestra = c("P01", "P01", "P02"), gen = c("BTK", "RAG1", "BTK"), chr = c("X", "11", "X"),
    pos = c(101360000L, 36573000L, 101360000L), id = c("rs1", NA, "rs1"), ref = "C", alt = "T",
    genotipo = c("hemi", "het", "het"), impacto = "HIGH", consecuencia = "stop_gained",
    af_gnomad = NA_real_, clinvar = "pathogenic",
    hgvs_c = c("ENST1:c.1A>T", "ENST2:c.2A>T", "ENST1:c.1A>T"),
    hgvs_p = c("ENSP1:p.Arg1Ter", "ENSP2:p.Lys2Ter", "ENSP1:p.Arg1Ter"))
  attr(x, "build") <- "GRCh38"
  attr(x, "ieiprio_sexo") <- c(P01 = "M", P02 = "F")
  suppressMessages(iei_priorizar(x, fenotipo = c("HP:0004432", "HP:0006532"), fenotipo_ausente = "HP:0000964"))
}

test_that("se exporta a TSV y la extension decide el formato", {
  p <- priorizado_exp()
  f <- withr::local_tempfile(fileext = ".tsv")
  expect_message(iei_exportar(p, f), "tsv")
  t <- utils::read.delim(f, stringsAsFactors = FALSE)
  expect_equal(nrow(t), 3)
  expect_true(all(c("gen", "puntaje", "acmg_clase") %in% names(t)))
  expect_error(iei_exportar(p, "x.doc"), "extension")
})

test_that("se exporta a Excel con hojas de parametros y fuentes", {
  skip_if_not_installed("writexl")
  skip_if_not_installed("readxl")
  f <- withr::local_tempfile(fileext = ".xlsx")
  suppressMessages(iei_exportar(priorizado_exp(), f))
  expect_setequal(readxl::excel_sheets(f), c("variantes", "parametros", "fuentes"))
  par <- readxl::read_excel(f, "parametros")
  expect_match(par$valor[par$parametro == "fenotipo ausente"], "HP:0000964")
})

test_that("el VCF exportado tiene una columna por muestra y la anotacion en INFO", {
  f <- withr::local_tempfile(fileext = ".vcf")
  suppressMessages(iei_exportar(priorizado_exp(), f))
  l <- readLines(f)
  expect_equal(l[1], "##fileformat=VCFv4.2")
  expect_true(any(grepl("^##INFO=<ID=ACMG_CLASS", l)))
  cab <- strsplit(l[grepl("^#CHROM", l)], "\t")[[1]]
  expect_equal(cab[10:11], c("P01", "P02"))
  datos <- l[!startsWith(l, "#")]
  expect_length(datos, 2)
  expect_match(datos[1], "^chr11\t36573000")
  btk <- strsplit(datos[2], "\t")[[1]]
  expect_equal(btk[c(1, 3, 10, 11)], c("chrX", "rs1", "1", "0/1"))
  expect_match(btk[8], "GENE=BTK")
  # Vuelve a leerse con iei_leer_vcf.
  v <- suppressWarnings(suppressMessages(iei_leer_vcf(f)))
  expect_equal(sort(unique(v$muestra)), c("P01", "P02"))
})

test_that("el BED usa coordenadas desde cero", {
  f <- withr::local_tempfile(fileext = ".bed")
  suppressMessages(iei_exportar(priorizado_exp(), f))
  b <- utils::read.delim(f, header = FALSE)
  expect_equal(b$V2[1], 101360000L - 1L)
  expect_equal(b$V3[1] - b$V2[1], 1L)
  expect_match(b$V4[1], "BTK:p.Arg1Ter")
})

test_that("el Phenopacket lleva fenotipo, ausentes y variantes candidatas", {
  p <- priorizado_exp()
  expect_error(iei_phenopacket(p), "indica")
  pp <- iei_phenopacket(p, muestra = "P01")
  expect_equal(pp$subject$sex, "MALE")
  expect_equal(pp$metaData$phenopacketSchemaVersion, "2.0")
  ids <- vapply(pp$phenotypicFeatures, function(r) r$type$id, "")
  expect_equal(ids, c("HP:0004432", "HP:0006532", "HP:0000964"))
  expect_true(isTRUE(pp$phenotypicFeatures[[3]]$excluded))
  gi <- pp$interpretations[[1]]$diagnosis$genomicInterpretations
  expect_gte(length(gi), 1)
  vd <- gi[[1]]$variantInterpretation$variationDescriptor
  expect_equal(vd$geneContext$symbol, "BTK")
  expect_equal(vd$allelicState$id, "GENO:0000134")
  expect_match(pp$interpretations[[1]]$diagnosis$disease$id, "^(OMIM|ORPHA):")
  f <- withr::local_tempfile(fileext = ".json")
  suppressMessages(iei_exportar(p, f, muestra = "P01"))
  j <- jsonlite::read_json(f)
  expect_equal(j$subject$id, "P01")
  expect_error(iei_phenopacket(p, muestra = "P09"), "no esta")
})

test_that("iei_granges pide GenomicRanges", {
  skip_if(nzchar(system.file(package = "GenomicRanges")))
  expect_error(iei_granges(priorizado_exp()), "GenomicRanges")
})

test_that("iei_granges devuelve un GRanges", {
  skip_if_not_installed("GenomicRanges")
  gr <- iei_granges(priorizado_exp())
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
})
