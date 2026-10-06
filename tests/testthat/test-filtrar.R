coords_prueba <- function() {
  tibble::tibble(
    gen = c("GENA", "BTK", "BTK"),
    build = c("GRCh38", "GRCh38", "GRCh37"),
    chr = c("1", "X", "X"),
    inicio = c(900L, 101349000L, 100604435L),
    fin = c(3100L, 101391000L, 100641212L),
    ensembl_id = NA_character_
  )
}
vcf_mini <- function() {
  iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
}
filtrar_mini <- function(...) {
  iei_filtrar(vcf_mini(), genes = c("GENA", "BTK"), coordenadas = coords_prueba(), ...)
}

test_that("aplica los tres pasos y registra la bitacora", {
  f <- filtrar_mini()
  b <- iei_bitacora(f)
  expect_equal(b$paso, c("Inicial", "Genotipo", "Panel", "Calidad"))
  expect_equal(b$filas, c(20L, 15L, 8L, 7L))
  expect_equal(b$descartadas, c(0L, 5L, 7L, 1L))
  expect_equal(nrow(f), 7)
})

test_that("agrega la columna gen", {
  f <- filtrar_mini()
  expect_equal(names(f)[1:2], c("muestra", "gen"))
  expect_setequal(unique(f$gen), c("GENA", "BTK"))
  expect_true(all(f$gen[f$chr == "X"] == "BTK"))
})

test_that("quita la variante de baja calidad", {
  f <- filtrar_mini()
  expect_false(3000 %in% f$pos)
  expect_true(3000 %in% filtrar_mini(min_dp = 5, solo_pass = FALSE)$pos)
})

test_that("el margen amplia la busqueda", {
  expect_false(4000 %in% filtrar_mini()$pos)
  expect_true(4000 %in% filtrar_mini(margen = 1000)$pos)
})

test_that("conserva hemicigotos de BTK y respeta la fraccion alelica", {
  f <- filtrar_mini()
  expect_equal(f$genotipo[f$muestra == "P01" & f$pos == 101360000], "hemi")
  expect_equal(f$genotipo[f$muestra == "P02" & f$pos == 101360000], "het")
  estricto <- filtrar_mini(min_ab = 0.49)
  expect_false(any(estricto$genotipo == "het" & !is.na(estricto$ab)))
})

test_that("usa las coordenadas del build del VCF", {
  v <- vcf_mini()
  attr(v, "build") <- "GRCh37"
  f <- iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba())
  expect_equal(nrow(f), 0)
  expect_equal(iei_bitacora(f)$filas[3], 0L)
})

test_that("avisa de genes sin coordenadas", {
  expect_message(
    iei_filtrar(vcf_mini(), genes = c("BTK", "SINCOORD"), coordenadas = coords_prueba()),
    "sin coordenadas"
  )
})

test_that("valida entradas", {
  expect_error(iei_filtrar(data.frame(a = 1)), "iei_leer_vcf")
  expect_error(filtrar_mini(min_ab = 0.6), "min_ab")
  expect_error(iei_bitacora(vcf_mini()), "bitacora")
})

test_that("funciona con los datos del paquete", {
  f <- iei_filtrar(vcf_mini())
  expect_true("gen" %in% names(f))
  expect_equal(iei_bitacora(f)$filas[1], 20L)
})
