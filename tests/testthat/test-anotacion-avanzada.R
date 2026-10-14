test_that("el recorte de alelos deja una base de anclaje", {
  r <- .recortar_alelos(c(100L, 200L, 300L), c("GCAT", "AC", "A"), c("GCT", "GC", "G"))
  expect_equal(r$pos, c(101L, 200L, 300L))
  expect_equal(r$ref, c("CA", "A", "A"))
  expect_equal(r$alt, c("C", "G", "G"))
  expect_equal(.es_indel(c("CA", "A", "<DEL>"), c("C", "G", "A")), c(TRUE, FALSE, FALSE))
})

test_that("las indels se alinean a la izquierda con la referencia", {
  withr::local_options(ieiprio.sin_red = FALSE)
  # Referencia desde la posicion 1: GCTAAAAGCGTACGT. Borrar una A del tramo
  # AAAA escrito al final (pos 6, AA>A) debe quedar como pos 3, TA>T.
  local_mocked_bindings(.referencia_consultar = function(chr, inicio, fin, build) {
    list(seq = substr("GCTAAAAGCGTACGTACGTACGT", inicio, fin))
  })
  x <- tibble::tibble(muestra = "P01", chr = "1", pos = c(6L, 10L), ref = c("AA", "G"), alt = c("A", "T"))
  msj <- testthat::capture_messages(n <- .normalizar_variantes(x, "GRCh38"))
  expect_true(any(grepl("normalizada", msj)))
  expect_equal(n$pos, c(3L, 10L))
  expect_equal(n$ref, c("TA", "G"))
  expect_equal(n$alt, c("T", "T"))
  expect_equal(n$variante_original, c("1:6:AA>A", NA))
})

test_that("sin red solo se recorta y se avisa", {
  withr::local_options(ieiprio.sin_red = TRUE)
  x <- tibble::tibble(muestra = "P01", chr = "1", pos = 6L, ref = "AAT", alt = "AT")
  expect_message(n <- .normalizar_variantes(x, "GRCh38"), "sin red")
  expect_equal(n$ref, "AA")
  expect_equal(n$alt, "A")
})

test_that("grpmax excluye los grupos con efecto fundador o sin ascendencia definida", {
  f <- list(gnomade = 5e-4, gnomade_afr = 1e-3, gnomade_fin = 0.5, gnomadg_nfe = 2e-3,
            gnomadg_asj = 0.3, gnomadg_remaining = 0.2)
  g <- .grpmax(f)
  expect_equal(g$af, 2e-3)
  expect_equal(g$grupo, "nfe")
  expect_true(is.na(.grpmax(list(gnomade = 1e-3))$af))
  expect_true(is.na(.grpmax(list(gnomade_fin = 0.1))$grupo))
})

test_that("la frecuencia de referencia prefiere FAF95, luego grpmax y luego la global", {
  x <- tibble::tibble(af_gnomad = c(1e-5, 1e-5, 1e-5, NA),
                      af_grpmax = c(NA, 3e-4, 3e-4, NA),
                      faf95 = c(NA, NA, 2e-3, NA))
  expect_equal(.af_referencia(x), c(1e-5, 3e-4, 2e-3, NA))
  expect_equal(.af_referencia(x[, "af_gnomad"]), x$af_gnomad)
})

test_that("iei_filtrar_frecuencia usa la frecuencia del grupo mas alto", {
  gen_ad <- stats::na.omit(iuis_panel$gen[iuis_panel$herencia %in% "AD"])[1]
  x <- tibble::tibble(gen = gen_ad, af_gnomad = c(1e-6, 1e-6), af_grpmax = c(NA, 5e-3))
  expect_equal(nrow(iei_filtrar_frecuencia(x)), 1)
})

test_that("se puede anotar con la salida de un VEP local", {
  local_mocked_bindings(.vep_consultar = vep_falso)
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  f <- iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba())
  entrada <- withr::local_tempfile(fileext = ".txt")
  expect_message(iei_vep_entrada(f, entrada), "escritas")
  lineas <- readLines(entrada)
  expect_true(all(grepl("^X [0-9]+ \\. [ACGT]+ [ACGT]+$", lineas)))
  # Un VEP local escribe un objeto JSON por linea, con la entrada sin los puntos finales.
  json <- withr::local_tempfile(fileext = ".json")
  objetos <- lapply(vep_sintetico(), function(o) { o$input <- sub(" \\. \\. \\.$", "", o$input); o })
  writeLines(vapply(objetos, jsonlite::toJSON, "", auto_unbox = TRUE, null = "null"), json)
  local <- suppressMessages(suppressWarnings(iei_anotar(f, cache = FALSE, vep_json = json)))
  remoto <- suppressMessages(suppressWarnings(iei_anotar(f, cache = FALSE)))
  comunes <- c("gen", "consecuencia", "hgvs_p", "af_gnomad", "clinvar")
  expect_equal(as.data.frame(local[comunes]), as.data.frame(remoto[comunes]), ignore_attr = TRUE)
  expect_equal(attr(local, "vep_release"), "local")
  expect_error(iei_anotar(f, cache = FALSE, vep_json = "no_existe.json"), "No encuentro")
})

lista_btk <- function() {
  data.frame(variante = c("NM_000061.3(BTK):c.1630A>G (p.Arg544Gly)",
                          "NM_000061.3(BTK):c.1631G>C (p.Arg544Pro)",
                          "NM_000061.3(BTK):c.1625T>C (p.Leu542Pro)",
                          "NM_000061.3(BTK):c.1636A>G (p.Met546Val)",
                          "NM_000061.3(BTK):c.1639C>T (p.Arg547Trp)",
                          "NM_000061.3(BTK):c.1581C>T (p.Gln527Ter)"),
             estrellas = c(2L, 1L, 1L, 1L, 1L, 0L), stringsAsFactors = FALSE)
}

test_that("los criterios de ClinVar del gen marcan PS1, PM5 y PM1", {
  expect_equal(.partir_p("ENSP1:p.Arg544Gly"), data.frame(ref = "Arg", pos = 544L, alt = "Gly"))
  expect_true(is.na(.partir_p("p.?")$pos))
  expect_equal(.c_de_titulo("NM_1(BTK):c.1630A>G (p.Arg544Gly)"), "c.1630A>G")
  expect_equal(.p_de_titulo("NM_1(BTK):c.1630A>G (p.Arg544Gly)"), "p.Arg544Gly")
  x <- tibble::tibble(
    gen = c("BTK", "BTK", "BTK", "OTRO"),
    hgvs_c = c("ENST1:c.1630A>G", "ENST1:c.1632G>C", "ENST1:c.898C>T", "ENST2:c.1A>G"),
    hgvs_p = c("ENSP1:p.Arg544Gly", "ENSP1:p.Arg544Gly", "ENSP1:p.Arg300Trp", "ENSP2:p.Met1Val"),
    consecuencia = "missense_variant",
    clinvar = c("pathogenic", NA, NA, "pathogenic")
  )
  cr <- .criterios_clinvar(x, list(BTK = lista_btk()))
  expect_equal(cr$estrellas, c(2L, NA, NA, NA))
  expect_equal(cr$ps1, c(FALSE, TRUE, FALSE, FALSE))
  expect_equal(cr$pm5, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(cr$pm1, c(TRUE, TRUE, FALSE, FALSE))
})

test_that("la anotacion de ClinVar usa la cache por gen y deja la fecha", {
  withr::local_options(ieiprio.sin_red = FALSE, ieiprio.silenciar_privacidad = TRUE)
  llamadas <- 0
  local_mocked_bindings(.clinvar_gen = function(g) { llamadas <<- llamadas + 1; list(variantes = lista_btk()) })
  x <- tibble::tibble(gen = "BTK", hgvs_c = "ENST1:c.1630A>G", hgvs_p = "ENSP1:p.Arg544Gly",
                      consecuencia = "missense_variant", clinvar = "pathogenic")
  dir <- withr::local_tempdir()
  a <- .anotar_clinvar(x, dir)
  b <- .anotar_clinvar(x, dir)
  expect_equal(llamadas, 1)
  expect_equal(a$clinvar_estrellas, 2L)
  expect_true(a$pm5)
  expect_equal(attr(a, "ieiprio_clinvar_fecha"), format(Sys.Date()))
  # Sin red no consulta y deja los criterios en FALSE.
  withr::local_options(ieiprio.sin_red = TRUE)
  s <- suppressMessages(.anotar_clinvar(x, FALSE))
  expect_false(s$pm5)
})

test_that("iei_gnomad agrega homocigotos, FAF95 y restriccion del gen", {
  withr::local_options(ieiprio.sin_red = FALSE, ieiprio.silenciar_privacidad = TRUE)
  local_mocked_bindings(.gnomad_consulta = function(query, variables) {
    if (!is.null(variables$id)) {
      if (variables$id == "X-100-C-T") return(list(data = list(variant = NULL)))
      return(list(data = list(variant = list(
        exome = list(ac = 10, an = 1e5, homozygote_count = 1, hemizygote_count = 3,
                     faf95 = list(popmax = 2e-4, popmax_population = "afr")),
        genome = list(ac = 2, an = 1e4, homozygote_count = 1, hemizygote_count = 0,
                      faf95 = list(popmax = 5e-4, popmax_population = "nfe"))))))
    }
    list(data = list(gene = list(gnomad_constraint = list(oe_lof_upper = 0.25, mis_z = 3.5, pli = 1))))
  })
  x <- tibble::tibble(gen = "BTK", chr = "X", pos = c(100L, 200L), ref = "C", alt = "T")
  attr(x, "build") <- "GRCh38"
  g <- suppressMessages(iei_gnomad(x, cache = FALSE))
  expect_equal(g$gnomad_hom, c(0L, 2L))
  expect_equal(g$gnomad_hemi, c(0L, 3L))
  expect_equal(g$faf95, c(NA, 5e-4))
  expect_equal(g$grupo_faf95, c(NA, "nfe"))
  expect_equal(g$loeuf, c(0.25, 0.25))
  expect_equal(g$mis_z, c(3.5, 3.5))
  expect_equal(.gnomad_dataset("GRCh37"), "gnomad_r2_1")
  expect_warning(suppressMessages(iei_gnomad(x, max_variantes = 1, cache = FALSE)), "primeras 1")
})

test_that("gnomAD sin red no falla y deja NA", {
  withr::local_options(ieiprio.sin_red = TRUE)
  x <- tibble::tibble(gen = "BTK", chr = "X", pos = 100L, ref = "C", alt = "T")
  g <- suppressMessages(iei_gnomad(x, cache = FALSE))
  expect_true(is.na(g$gnomad_hom))
  expect_true(is.na(g$loeuf))
})
