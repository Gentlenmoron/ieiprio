# Instantaneas de resultados clave: si un cambio en el codigo o en los datos
# mueve el ranking o los criterios, la prueba lo muestra para revisarlo.

test_that("el ranking por fenotipo de un caso clasico no cambia sin querer", {
  b <- suppressMessages(iei_buscar_fenotipo(c("HP:0004432", "HP:0006532"), n = 10))
  expect_snapshot_value(data.frame(gen = b$gen, similitud = b$similitud), style = "json2")
})

test_that("los criterios ACMG de un conjunto fijo no cambian sin querer", {
  x <- tibble::tibble(
    consecuencia = c("stop_gained", "missense_variant", "synonymous_variant", "splice_donor_variant"),
    af_gnomad = c(NA, 3e-6, 0.03, NA), exon = c("3/10", "5/10", "6/10", NA),
    ps1 = c(FALSE, TRUE, FALSE, FALSE), spliceai = c(NA, NA, 0.01, 0.9))
  a <- .acmg_sugeridos(x, rep(FALSE, 4), list(compuesto = rep(FALSE, 4)), c(3, 1, 0, 2.7),
                       list(favor = c(FALSE, TRUE, FALSE, TRUE), contra = c(FALSE, FALSE, TRUE, FALSE)),
                       c("AR", "AD", "AR", "XL"))
  expect_snapshot_value(a, style = "json2")
})
