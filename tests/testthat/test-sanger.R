# Escenario comun: referencia sintetica en X:5000000 (fuera de la region
# pseudoautosomica) y lecturas que cubren las bases 1001 a 1500.
INICIO <- 5000000L
REGION <- sprintf("X:%d-%d", INICIO, INICIO + 2999L)

ref_falsa <- function(ref) {
  function(chr, inicio, fin, build) list(seq = substr(ref, inicio - INICIO + 1, fin - INICIO + 1))
}

preparar_lecturas <- function(dir) {
  ref <- referencia_prueba()
  b <- strsplit(substr(ref, 1001, 1500), "")[[1]]
  alt_hom <- setdiff(c("A", "C", "G", "T"), b[100])[1]
  b[100] <- alt_hom
  segunda <- setdiff(c("A", "C", "G", "T"), b[250])[1]
  b <- b[-(350:351)]
  s <- paste(b, collapse = "")
  cal <- c(rep(5L, 15), rep(40L, length(b) - 30), rep(5L, 15))
  escribir_ab1(file.path(dir, "F.ab1"), s, het = list("250" = segunda), calidad = cal)
  n <- nchar(s)
  het_r <- list()
  het_r[[as.character(n - 250 + 1)]] <- c(A = "T", C = "G", G = "C", T = "A")[[segunda]]
  escribir_ab1(file.path(dir, "R.ab1"), revcomp_prueba(s), het = het_r, calidad = rev(cal))
  list(ref = ref, s = s, alt_hom = alt_hom, segunda = segunda)
}

test_that("iei_leer_ab1 lee bases, calidad, picos y trazas", {
  d <- withr::local_tempdir()
  escribir_ab1(file.path(d, "a.ab1"), "ACGTTGCA", calidad = c(10L, 20L, 30L, 40L, 40L, 30L, 20L, 10L),
               muestra = "P01")
  a <- iei_leer_ab1(file.path(d, "a.ab1"))
  expect_s3_class(a, "iei_ab1")
  expect_equal(a$bases, "ACGTTGCA")
  expect_equal(a$calidad, c(10L, 20L, 30L, 40L, 40L, 30L, 20L, 10L))
  expect_equal(a$picos, 12L * (1:8) + 1L)
  expect_named(a$trazas, c("A", "C", "G", "T"))
  # En el pico de la primera base (A) domina la traza A.
  expect_equal(names(which.max(unlist(a$trazas[a$picos[1], ]))), "A")
  expect_equal(a$muestra, "P01")
  expect_match(paste(cli::cli_fmt(print(a)), collapse = " "), "8 bases")
})

test_that("iei_leer_ab1 rechaza archivos que no son ABIF", {
  f <- withr::local_tempfile(fileext = ".ab1")
  writeLines("no soy un cromatograma", f)
  expect_error(iei_leer_ab1(f), "ABIF")
  expect_error(iei_leer_ab1("no_existe.ab1"), "No encuentro")
})

test_that("los picos dobles se llaman heterocigotos y el ruido no", {
  d <- withr::local_tempdir()
  s <- referencia_prueba(200, semilla = 3)
  segunda <- setdiff(c("A", "C", "G", "T"), substr(s, 100, 100))[1]
  escribir_ab1(file.path(d, "h.ab1"), s, het = list("100" = segunda))
  a <- iei_leer_ab1(file.path(d, "h.ab1"))
  h <- .llamar_heterocigotos(a, 0.33)
  b <- strsplit(h$bases, "")[[1]]
  expect_equal(which(b != strsplit(s, "")[[1]]), 100L)
  expect_equal(b[100], .iupac_de(substr(s, 100, 100), segunda))
  expect_gt(h$razon[100], 0.7)
  # Un umbral mas alto que la senal secundaria no llama nada.
  expect_equal(.llamar_heterocigotos(a, 0.9)$bases, s)
})

test_that("el recorte quita los extremos de baja calidad", {
  q <- c(rep(5L, 10), rep(40L, 50), rep(3L, 10))
  expect_equal(range(.recortar(q)), c(11L, 60L))
  expect_equal(.recortar(rep(NA_integer_, 5)), 1:5)
  expect_length(.recortar(rep(2L, 20)), 0)
})

test_that("iei_sanger encuentra SNV homocigota, heterocigota y delecion en ambas hebras", {
  d <- withr::local_tempdir()
  p <- preparar_lecturas(d)
  local_mocked_bindings(.referencia_consultar = ref_falsa(p$ref))
  v <- suppressMessages(iei_sanger(file.path(d, c("F.ab1", "R.ab1")), region = REGION, cache = FALSE))
  expect_equal(nrow(v), 3)
  expect_equal(v$pos, INICIO + 1000L + c(99L, 249L, 348L))
  expect_equal(v$genotipo, c("hom_alt", "het", "hom_alt"))
  expect_equal(v$alt[1], p$alt_hom)
  expect_equal(v$alt[2], p$segunda)
  expect_equal(nchar(v$ref[3]), 3)
  expect_equal(v$lecturas, c(2L, 2L, 2L))
  expect_equal(v$ab[1], 1)
  expect_equal(v$ab[2], 800 / 1800, tolerance = 0.01)
  expect_equal(unique(v$fuente), "sanger")
  info <- attr(v, "ieiprio_sanger")
  expect_equal(unname(vapply(info, `[[`, "", "hebra")), c("+", "-"))
  expect_true(all(vapply(info, `[[`, 0, "identidad") > 98))
  expect_equal(attr(v, "build"), "GRCh38")
})

test_that("en un varon la homocigota del X fuera de la PAR pasa a hemicigota", {
  d <- withr::local_tempdir()
  p <- preparar_lecturas(d)
  local_mocked_bindings(.referencia_consultar = ref_falsa(p$ref))
  v <- suppressMessages(iei_sanger(file.path(d, "F.ab1"), region = REGION, sexo = "M", cache = FALSE))
  expect_equal(v$genotipo, c("hemi", "het", "hemi"))
  expect_equal(v$lecturas, c(1L, 1L, 1L))
})

test_that("acepta FASTA y secuencia escrita, en cualquier hebra", {
  d <- withr::local_tempdir()
  p <- preparar_lecturas(d)
  local_mocked_bindings(.referencia_consultar = ref_falsa(p$ref))
  fa <- file.path(d, "lectura.fasta")
  writeLines(c(">lectura", substring(p$s, seq(1, nchar(p$s), 60), pmin(seq(60, nchar(p$s) + 59, 60), nchar(p$s)))), fa)
  v1 <- suppressMessages(iei_sanger(fa, region = REGION, cache = FALSE))
  v2 <- suppressMessages(iei_sanger(revcomp_prueba(p$s), region = REGION, cache = FALSE))
  expect_equal(v1$pos, v2$pos)
  expect_equal(v1$genotipo, c("hom_alt", "hom_alt"))
  expect_true(all(is.na(v1$ab)))
  expect_equal(attr(v2, "ieiprio_sanger")[[1]]$hebra, "-")
  expect_error(iei_sanger("ACGT", region = REGION, cache = FALSE), "no es un archivo")
})

test_that("las deleciones en repeticiones se normalizan a la izquierda", {
  ref <- list(inicio = 100L, secuencia = "GCATTTTTGCA")
  vars <- data.frame(pos = 105L, ref = "TT", alt = "T")
  n <- .normalizar_indels(vars, ref)
  expect_equal(n$pos, 102L)
  expect_equal(n$ref, "AT")
  expect_equal(n$alt, "A")
  ins <- .normalizar_indels(data.frame(pos = 106L, ref = "T", alt = "TT"), ref)
  expect_equal(c(ins$pos, ins$ref, ins$alt), c("102", "A", "AT"))
})

test_that("avisa de una posible indel heterocigota", {
  d <- withr::local_tempdir()
  ref <- referencia_prueba()
  s <- substr(ref, 1001, 1400)
  b <- strsplit(s, "")[[1]]
  het <- list()
  for (j in 200:(length(b) - 1)) if (b[j + 1] != b[j]) het[[as.character(j)]] <- b[j + 1]
  escribir_ab1(file.path(d, "i.ab1"), s, het = het)
  local_mocked_bindings(.referencia_consultar = ref_falsa(ref))
  # El segundo pico es siempre la base siguiente: una delecion de una base
  # heterocigota, que ahora se llama por deconvolucion.
  expect_message(v <- iei_sanger(file.path(d, "i.ab1"), region = REGION, cache = FALSE),
                 "indel heterocigota")
  pos <- attr(v, "ieiprio_sanger")[[1]]$indel_heterocigota
  expect_true(abs(pos - (INICIO + 1000L + 199L)) <= 3)
  del <- v[nchar(v$ref) == 2 & nchar(v$alt) == 1, ]
  expect_equal(nrow(del), 1)
  expect_equal(del$genotipo, "het")
  expect_true(abs(del$pos - (INICIO + 1000L + 198L)) <= 3)
  # Sin deconvolucion queda el aviso.
  expect_warning(suppressMessages(iei_sanger(file.path(d, "i.ab1"), region = REGION, cache = FALSE,
                                             indels_het = FALSE)),
                 "picos dobles sostenidos")
})

test_that("una lectura que no es del gen se descarta con aviso", {
  local_mocked_bindings(.referencia_consultar = ref_falsa(referencia_prueba()))
  ajena <- referencia_prueba(300, semilla = 99)
  expect_error(
    expect_warning(iei_sanger(ajena, region = REGION, cache = FALSE), "no se pudo ubicar"),
    "Ninguna lectura"
  )
  expect_error(iei_sanger(ajena), "gen")
})

test_that("con gen usa las coordenadas del panel y guarda la referencia en cache", {
  d <- withr::local_tempdir()
  btk <- ieiprio::genes_coordenadas
  btk <- btk[btk$gen == "BTK" & btk$build == "GRCh38", ]
  pedidos <- list()
  local_mocked_bindings(.referencia_consultar = function(chr, inicio, fin, build) {
    pedidos[[length(pedidos) + 1]] <<- c(chr, inicio, fin)
    list(seq = referencia_prueba(fin - inicio + 1))
  })
  lectura <- substr(referencia_prueba(btk$fin - btk$inicio + 2001), 5001, 5400)
  v1 <- suppressMessages(iei_sanger(lectura, gen = "btk", cache = d))
  v2 <- suppressMessages(iei_sanger(lectura, gen = "BTK", cache = d))
  expect_length(pedidos, 1)
  expect_equal(pedidos[[1]], c("X", btk$inicio - 1000, btk$fin + 1000))
  expect_equal(nrow(v1), 0)
  expect_equal(attr(v2, "ieiprio_sanger")[[1]]$region,
               sprintf("X:%d-%d", btk$inicio - 1000 + 5000, btk$inicio - 1000 + 5399))
})

test_that("iei_cromatograma dibuja las lecturas que cubren la variante", {
  skip_if_not_installed("ggplot2")
  d <- withr::local_tempdir()
  p <- preparar_lecturas(d)
  local_mocked_bindings(.referencia_consultar = ref_falsa(p$ref))
  v <- suppressMessages(iei_sanger(file.path(d, c("F.ab1", "R.ab1")), region = REGION, cache = FALSE))
  g <- iei_cromatograma(v, variante = 2)
  expect_s3_class(g, "ggplot")
  expect_match(g$labels$title, "het")
  expect_setequal(unique(g$data$lectura), c("F.ab1 (+)", "R.ab1 (-)"))
  expect_no_error(ggplot2::ggplot_build(iei_cromatograma(v, variante = 3)))
  expect_equal(unique(iei_cromatograma(v, posicion = v$pos[1], lectura = "R.ab1")$data$lectura), "R.ab1 (-)")
  expect_error(iei_cromatograma(v, variante = 9), "solo tiene")
  expect_error(iei_cromatograma(tibble::tibble(pos = 1)), "iei_sanger")
})

test_that("Sanger en vivo: secuencia de BTK con un cambio introducido", {
  saltar_si_no_en_vivo()
  btk <- ieiprio::genes_coordenadas
  btk <- btk[btk$gen == "BTK" & btk$build == "GRCh38", ]
  ini <- btk$inicio + 10000L
  r <- en_vivo(.referencia_consultar("X", ini, ini + 399L, "GRCh38"))
  b <- strsplit(toupper(r$seq), "")[[1]]
  b[200] <- setdiff(c("A", "C", "G", "T"), b[200])[1]
  lectura <- revcomp_prueba(paste(b, collapse = ""))
  v <- en_vivo(suppressMessages(iei_sanger(lectura, gen = "BTK", sexo = "M", cache = withr::local_tempdir())))
  expect_equal(v$pos, ini + 199L)
  expect_equal(v$genotipo, "hemi")
  a <- en_vivo(suppressMessages(iei_anotar(v, cache = FALSE)))
  expect_equal(a$gen, "BTK")
})

test_that("la salida de iei_sanger recorre anotar, priorizar y reporte", {
  d <- withr::local_tempdir()
  p <- preparar_lecturas(d)
  local_mocked_bindings(.referencia_consultar = ref_falsa(p$ref),
                        .vep_consultar = function(entrada, build) list(),
                        .vep_release = function(build) "999")
  v <- suppressMessages(iei_sanger(file.path(d, c("F.ab1", "R.ab1")), gen = "BTK",
                                   region = REGION, sexo = "M", cache = FALSE))
  expect_equal(unique(v$gen), "BTK")
  a <- suppressMessages(iei_anotar(v, cache = FALSE))
  r <- suppressMessages(iei_priorizar(a))
  expect_s3_class(r, "iei_resultado")
  expect_equal(nrow(r), 3)
  expect_false(is.null(attr(r, "ieiprio_sanger")))
  expect_s3_class(iei_cromatograma(r, variante = 1), "ggplot")
  html <- withr::local_tempfile(fileext = ".html")
  expect_no_error(suppressMessages(iei_reporte(r, archivo = html)))
  texto <- paste(readLines(html, warn = FALSE), collapse = " ")
  expect_match(texto, "Lecturas Sanger")
  expect_match(texto, "R.ab1: hebra -")
})

test_that("un archivo que no existe da un mensaje claro", {
  expect_error(iei_sanger("ruta/F.ab1", region = REGION, cache = FALSE), "No encuentro el archivo")
  expect_error(.sanger_leer("ACGTNNXZ", 0.33), "no son bases")
  expect_error(.sanger_leer("ACGTACGT", 0.33), "al menos 30")
})

ref_btk <- function() system.file("extdata", "btk_referencia.fa", package = "ieiprio")

test_that("iei_ejemplo_sanger genera lecturas que reproducen la variante de BTK sin red", {
  local_mocked_bindings(.referencia_consultar = function(...) stop("no deberia consultar Ensembl"))
  ab1 <- suppressMessages(iei_ejemplo_sanger(dir = withr::local_tempdir()))
  expect_length(ab1, 2)
  expect_true(all(file.exists(ab1)))
  v <- suppressMessages(iei_sanger(ab1, gen = "BTK", sexo = "M", referencia = ref_btk()))
  expect_s3_class(v, "iei_sanger")
  expect_equal(nrow(v), 1)
  expect_equal(c(v$chr, v$pos, v$ref, v$alt, v$genotipo), c("X", "101353990", "T", "C", "hemi"))
  expect_equal(v$lecturas, 2L)

  ab1_het <- suppressMessages(iei_ejemplo_sanger("het", dir = withr::local_tempdir()))
  h <- suppressMessages(iei_sanger(ab1_het, gen = "BTK", sexo = "F", referencia = ref_btk()))
  expect_equal(c(h$pos, h$genotipo), c("101353990", "het"))
  expect_equal(h$lecturas, 2L)
  expect_s3_class(iei_cromatograma(h), "ggplot")
})

test_that("la referencia propia acepta un FASTA con region o una lista", {
  r <- .referencia_propia(ref_btk(), gen = "BTK")
  expect_equal(c(r$chr, r$inicio, r$fin), c("X", "101353540", "101354439"))
  expect_equal(substr(r$secuencia, 101353990 - r$inicio + 1, 101353990 - r$inicio + 1), "T")
  l <- .referencia_propia(list(secuencia = "acgt acgt", inicio = 10), gen = "BTK")
  expect_equal(c(l$chr, l$fin, l$secuencia), c("X", "17", "ACGTACGT"))
  sin_region <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">BTK", "ACGT"), sin_region)
  expect_error(.referencia_propia(sin_region), "regi")
  expect_error(.referencia_propia(list(secuencia = "ACGT", inicio = 1)), "cromosoma")
  expect_error(iei_sanger("ACGT"), "referencia")
})

test_that("una secuencia escrita se analiza con la referencia del paquete", {
  r <- .referencia_propia(ref_btk(), gen = "BTK")
  b <- strsplit(substr(r$secuencia, 301, 600), "")[[1]]
  b[101353990 - r$inicio + 1 - 300] <- "C"
  v <- suppressMessages(iei_sanger(paste(b, collapse = ""), gen = "BTK", sexo = "M", referencia = ref_btk()))
  expect_equal(c(v$pos, v$ref, v$alt, v$genotipo), c("101353990", "T", "C", "hemi"))
})

test_that("el cromatograma clasico usa theme_minimal y el estilo iei el tema del paquete", {
  skip_if_not_installed("ggplot2")
  v <- suppressMessages(iei_sanger(suppressMessages(iei_ejemplo_sanger("het", dir = withr::local_tempdir())),
                                   gen = "BTK", referencia = ref_btk()))
  clasico <- iei_cromatograma(v)
  iei <- iei_cromatograma(v, estilo = "iei")
  expect_null(clasico$theme$plot.title$colour)
  expect_equal(iei$theme$plot.title$colour, .marca)
  expect_error(iei_cromatograma(v, estilo = "otro"))
})

test_that("iei_dotplot marca la diagonal, la variante y los desplazamientos por indel", {
  skip_if_not_installed("ggplot2")
  v <- suppressMessages(iei_sanger(suppressMessages(iei_ejemplo_sanger(dir = withr::local_tempdir())),
                                   gen = "BTK", sexo = "M", referencia = ref_btk()))
  g <- iei_dotplot(v)
  expect_s3_class(g, "ggplot")
  p <- attr(g, "puntos")
  expect_gt(nrow(p), 300)
  # La mayoria de los puntos cae en la diagonal principal.
  desfase <- p$j - p$i
  expect_gt(max(table(desfase)) / nrow(p), 0.5)
  expect_true(any(vapply(g$layers, function(l) inherits(l$geom, "GeomLabel"), logical(1))))
  expect_equal(g$theme$aspect.ratio, 3 / 4)
  expect_s3_class(plot(v, tipo = "dotplot"), "ggplot")
  expect_s3_class(plot(v), "ggplot")

  a <- "ACGTTGCAACGTAGGCTTACGATCGATCGGATCCATGCATTGACCA"
  b <- paste0(substr(a, 1, 20), "GGG", substring(a, 21))
  d <- attr(iei_dotplot(a, b, ventana = 6, umbral = 6), "puntos")
  expect_setequal(unique(d$j - d$i)[unique(d$j - d$i) %in% c(0, 3)], c(0, 3))
  propia <- attr(iei_dotplot(a, ventana = 6, umbral = 6), "puntos")
  expect_true(all(c(0) %in% (propia$j - propia$i)))
  inversa <- iei_dotplot(a, .revcomp(a), ventana = 6, umbral = 6)
  expect_match(inversa$labels$subtitle, "complemento reverso")
  expect_error(iei_dotplot(a, ventana = 5, umbral = 6), "umbral")
  expect_error(iei_dotplot("XYZ"), "v\u00e1lida")
  expect_error(iei_dotplot(tibble::tibble(a = 1)), "iei_sanger")
})

test_that("iei_dotplot ubica una secuencia en la referencia del gen", {
  skip_if_not_installed("ggplot2")
  r <- .referencia_propia(ref_btk(), gen = "BTK")
  s <- .revcomp(substr(r$secuencia, 301, 600))
  g <- iei_dotplot(s, gen = "BTK", referencia = ref_btk())
  p <- attr(g, "puntos")
  expect_true(all(p$j >= r$inicio))
  expect_match(g$labels$y, "Cromosoma X")
})

test_that("el reporte de un resultado Sanger incluye los cromatogramas y omite el control del VCF", {
  skip_if_not_installed("ggplot2")
  v <- suppressMessages(iei_sanger(suppressMessages(iei_ejemplo_sanger("het", dir = withr::local_tempdir())),
                                   gen = "BTK", sexo = "F", referencia = ref_btk()))
  vep <- system.file("extdata", "btk_vep.json", package = "ieiprio")
  r <- suppressMessages(iei_priorizar(iei_anotar(v, vep_json = vep, cache = FALSE), fenotipo = "agammaglobulinemia"))
  expect_equal(r$hgvs_c, "ENST00000308731.8:c.1632-2A>G")
  html <- withr::local_tempfile(fileext = ".html")
  suppressMessages(iei_reporte(r, html, dominios = FALSE))
  h <- paste(readLines(html, warn = FALSE), collapse = " ")
  expect_match(h, "<h2>Cromatogramas</h2>")
  expect_match(h, "data:image/png;base64,")
  expect_false(grepl("Control de calidad de la muestra", h))
  expect_match(h, "Sanger")
})

test_that("ejemplo Sanger en vivo llega a ClinVar patogenica", {
  saltar_si_no_en_vivo()
  ab1 <- en_vivo(suppressMessages(iei_ejemplo_sanger(dir = withr::local_tempdir())))
  v <- en_vivo(suppressMessages(iei_sanger(ab1, gen = "BTK", sexo = "M", cache = withr::local_tempdir())))
  r <- en_vivo(suppressMessages(iei_priorizar(iei_anotar(v, cache = FALSE), fenotipo = "HP:0004313")))
  expect_equal(r$categoria, "alta")
  expect_match(r$clinvar, "athogenic")
})
