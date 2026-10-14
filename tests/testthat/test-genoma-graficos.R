b_genoma <- function() suppressMessages(iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente"), n = 12))

test_that("el cariograma marca el panel, una busqueda y un resultado", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(iei_cariograma(), "ggplot")
  g <- iei_cariograma(b_genoma())
  expect_s3_class(g, "ggplot")
  expect_s3_class(plot(b_genoma(), tipo = "cariotipo"), "ggplot")
  expect_s3_class(iei_cariograma(build = "GRCh37"), "ggplot")
})

test_that("los graficos circulares dibujan sin error", {
  skip_if_not_installed("circlize")
  f <- withr::local_tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(iei_circos())
  expect_null(iei_circos(b_genoma()))
  expect_null(plot(b_genoma(), tipo = "acordes"))
})

test_that("los datos del genoma son coherentes", {
  g <- .genoma("GRCh38")
  expect_equal(nrow(g), 24)
  expect_true(all(g$centromero < g$largo))
  p <- .genes_panel_posicion("GRCh38")
  expect_true(all(p$chr %in% .cromosomas))
  expect_true(all(grepl("^#", p$color)))
})

test_that("el oncoplot y los nuevos estilos devuelven graficos", {
  skip_if_not_installed("ggplot2")
  anot <- tibble::tibble(muestra = c("A", "A", "B"), gen = c("BTK", "RAG1", "BTK"), chr = c("X", "11", "X"),
                         pos = c(101354013L, 36573500L, 101354900L), ref = "C", alt = "T",
                         genotipo = c("hemi", "het", "hemi"), impacto = c("HIGH", "MODERATE", "HIGH"),
                         consecuencia = c("stop_gained", "missense_variant", "splice_acceptor_variant"),
                         af_gnomad = NA_real_, clinvar = c("pathogenic", NA, NA))
  attr(anot, "ieiprio_sexo") <- c(A = "M", B = "M")
  prior <- suppressMessages(iei_priorizar(anot))
  expect_true(inherits(plot(prior, tipo = "oncoplot"), c("ggplot", "patchwork")))
  expect_s3_class(plot(prior, tipo = "genes"), "ggplot")
  expect_s3_class(plot(prior, tipo = "cariotipo"), "ggplot")
  expect_equal(.clase_consecuencia(c("stop_gained", "frameshift_variant", "splice_donor_variant",
                                     "missense_variant", "synonymous_variant", NA)),
               c("sin sentido", "frameshift o indel", "splicing", "missense", "sinonima", "otra"))
})
