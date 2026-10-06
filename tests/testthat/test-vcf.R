leer_mini <- function(...) {
  iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"), ...)
}
gen <- function(v, muestra, pos, alt) {
  v$genotipo[v$muestra == muestra & v$pos == pos & v$alt == alt]
}

test_that("lee el VCF con la estructura esperada", {
  v <- leer_mini()
  expect_s3_class(v, "tbl_df")
  expect_named(v, c("muestra", "chr", "pos", "id", "ref", "alt", "genotipo",
                    "qual", "dp", "gq", "ab", "filtro", "build"))
  # 10 alelos alternativos utiles x 2 muestras
  expect_equal(nrow(v), 20)
  expect_equal(attr(v, "build"), "GRCh38")
})

test_that("lee igual el archivo comprimido", {
  expect_equal(nrow(iei_leer_vcf(fixture("mini.vcf.gz"))), 20)
})

test_that("normaliza los cromosomas", {
  v <- leer_mini()
  expect_setequal(unique(v$chr), c("1", "X", "MT"))
})

test_that("separa multialelicos y recalcula el genotipo", {
  v <- leer_mini()
  expect_equal(gen(v, "P01", 2000, "T"), "het")
  expect_equal(gen(v, "P01", 2000, "G"), "het")
  expect_equal(gen(v, "P02", 2000, "T"), "het")
  expect_equal(gen(v, "P02", 2000, "G"), "hom_ref")
  expect_equal(v$ab[v$muestra == "P01" & v$pos == 2000 & v$alt == "G"], 9 / 17)
})

test_that("descarta el alelo de deleccion solapada (*)", {
  v <- leer_mini()
  expect_false(any(v$alt == "*"))
  expect_equal(gen(v, "P01", 4000, "C"), "het")
})

test_that("marca hemicigotos en X para varones fuera de la PAR", {
  v <- leer_mini()
  expect_equal(gen(v, "P01", 101360000, "T"), "hemi")
  expect_equal(gen(v, "P02", 101360000, "T"), "het")
  # PAR1: el varon conserva el genotipo diploide
  expect_equal(gen(v, "P01", 100000, "C"), "het")
  # Llamado haploide
  expect_equal(gen(v, "P01", 101370000, "G"), "hemi")
})

test_that("sin sexo no inventa hemicigotos", {
  v <- iei_leer_vcf(fixture("mini.vcf"))
  expect_equal(gen(v, "P01", 101360000, "T"), "hom_alt")
})

test_that("maneja genotipos faltantes, fases y campos ausentes", {
  v <- leer_mini()
  expect_true(is.na(gen(v, "P02", 3000, "A")))
  expect_true(is.na(v$gq[v$muestra == "P01" & v$pos == 3000]))
  expect_equal(gen(v, "P02", 5000, "A"), "hom_alt")
  expect_true(is.na(v$ab[v$muestra == "P01" & v$pos == 5000]))
  expect_equal(v$filtro[v$muestra == "P01" & v$pos == 3000], "LowQual")
  expect_equal(v$id[v$pos == 2000][1], "rs1")
  expect_true(is.na(v$id[v$pos == 1000][1]))
})

test_that("permite elegir muestras", {
  v <- iei_leer_vcf(fixture("mini.vcf"), muestras = "P02")
  expect_equal(unique(v$muestra), "P02")
  expect_error(iei_leer_vcf(fixture("mini.vcf"), muestras = "P99"), "P99")
})

test_that("detecta GRCh37 por la longitud del contig", {
  expect_warning(v <- iei_leer_vcf(fixture("b37.vcf")), "NA")
  expect_equal(attr(v, "build"), "GRCh37")
})

test_that("pide el build si no lo puede detectar", {
  expect_error(iei_leer_vcf(fixture("sin_build.vcf")), "build")
  expect_warning(v <- iei_leer_vcf(fixture("sin_build.vcf"), build = "GRCh38"))
  expect_equal(v$build[1], "GRCh38")
})

test_that("valida entradas", {
  expect_error(iei_leer_vcf("no_existe.vcf"), "No encuentro")
  expect_error(iei_leer_vcf(fixture("mini.vcf"), sexo = "X"), "solo acepta")
  expect_error(iei_leer_vcf(fixture("mini.vcf"), sexo = c("M", "F")), "nombres")
})

test_that("funciones internas", {
  expect_equal(.normalizar_chr(c("chr1", "CHRX", "chrM", "MT", "7")),
               c("1", "X", "MT", "MT", "7"))
  expect_equal(.clasificar_genotipo(c("0/1", "1/1", "0/0", "./.", "1", "0", "2|1"),
                                    c(1, 1, 1, 1, 1, 1, 2)),
               c("het", "hom_alt", "hom_ref", NA, "hemi", "hom_ref", "het"))
  expect_true(.en_par(50000, "GRCh38"))
  expect_false(.en_par(50000, "GRCh37"))
})
