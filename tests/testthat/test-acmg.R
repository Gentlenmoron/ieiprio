acmg_x <- function(...) {
  base <- list(consecuencia = "missense_variant", af_gnomad = NA_real_)
  tibble::as_tibble(utils::modifyList(base, list(...)))
}

acmg <- function(x, gof = FALSE, p_fenotipo = 0, herencia = "AR", h = NULL, pr = NULL) {
  n <- nrow(x)
  .acmg_sugeridos(x, rep(gof, n), h %||% list(compuesto = rep(FALSE, n)), rep(p_fenotipo, n),
                  pr, rep(herencia, n))
}

test_that("los puntos y las clases siguen a Tavtigian 2020", {
  expect_equal(.puntos_criterio(c("PVS1", "PS1", "PM2", "PP3", "BA1", "BS2", "BP4")),
               c(8, 4, 2, 1, -8, -4, -1))
  expect_equal(.puntos_criterio(c("PVS1_Strong", "PVS1_Moderate", "PM2_Supporting", "PM3_Supporting",
                                  "PP3_Strong", "PM1_Supporting")),
               c(4, 2, 1, 1, 4, 1))
  expect_equal(.clase_acmg(c(10, 6, 9, 0, 5, -1, -6, -7), rep(FALSE, 8)),
               c("Patog\u00e9nica", "Probablemente patog\u00e9nica", "Probablemente patog\u00e9nica",
                 "Significado incierto", "Significado incierto", "Probablemente benigna",
                 "Probablemente benigna", "Benigna"))
  expect_equal(.clase_acmg(8, TRUE), "Benigna")
})

test_that("PVS1 se gradua segun el tipo de variante y el exon", {
  x <- acmg_x(consecuencia = c("stop_gained", "frameshift_variant", "splice_donor_variant",
                               "start_lost", "missense_variant"),
              exon = c("3/10", "10/10", NA, "1/10", "4/10"))
  a <- acmg(x)
  pv <- vapply(strsplit(a$criterios, ", "), `[`, "", 1)
  expect_equal(pv, c("PVS1", "PVS1_Moderate", "PVS1_Strong", "PVS1_Moderate", "PM2_Supporting"))
  # En genes de ganancia de funcion pura no se aplica.
  expect_false(any(grepl("PVS1", acmg(x, gof = TRUE)$criterios)))
})

test_that("PM2, BS1 y BA1 dependen de la herencia del gen", {
  x <- acmg_x(af_gnomad = c(5e-6, 5e-5, 5e-3, 0.03, 0.1))
  dom <- acmg(x, herencia = "AD")
  rec <- acmg(x, herencia = "AR")
  expect_equal(grepl("PM2", dom$criterios), c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(grepl("PM2", rec$criterios), c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(grepl("BS1", dom$criterios), c(FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(grepl("BS1", rec$criterios), c(FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_equal(dom$clase[5], "Benigna")
  # FAF95 tiene prioridad sobre la frecuencia global.
  f <- acmg(acmg_x(af_gnomad = 1e-6, faf95 = 0.002), herencia = "AD")
  expect_match(f$criterios, "BS1")
})

test_that("PS1, PM5 y PM1 vienen de ClinVar; PS1 anula PM5", {
  x <- acmg_x(consecuencia = "missense_variant", ps1 = c(TRUE, FALSE, FALSE),
              pm5 = c(TRUE, TRUE, FALSE), pm1 = c(FALSE, FALSE, TRUE))
  a <- acmg(x)
  expect_match(a$criterios[1], "PS1")
  expect_false(grepl("PM5", a$criterios[1]))
  expect_match(a$criterios[2], "PM5")
  expect_match(a$criterios[3], "PM1_Supporting")
  expect_equal(a$puntos, c(4 + 1, 2 + 1, 1 + 1))
})

test_that("BS2 usa homocigotos o hemicigotos sanos de gnomAD", {
  x <- acmg_x(af_gnomad = 1e-4, gnomad_hom = c(3L, 0L, 0L), gnomad_hemi = c(0L, 4L, 1L))
  expect_equal(grepl("BS2", acmg(x, herencia = "AR")$criterios), c(TRUE, FALSE, FALSE))
  expect_equal(grepl("BS2", acmg(x, herencia = "XL")$criterios), c(FALSE, TRUE, FALSE))
})

test_that("PP2, BP7, PM4, PP3, PP4 y BP4", {
  x <- acmg_x(consecuencia = c("missense_variant", "synonymous_variant", "inframe_deletion",
                               "synonymous_variant"),
              mis_z = 3.5, spliceai = c(NA, 0.02, NA, 0.5))
  pr <- list(favor = c(TRUE, FALSE, FALSE, FALSE), contra = c(FALSE, TRUE, FALSE, FALSE))
  h <- list(compuesto = rep(FALSE, 4))
  a <- acmg(x, pr = pr, h = h, p_fenotipo = 2.7)
  expect_false(any(grepl("PP4", acmg(x, pr = pr, h = h, p_fenotipo = 2.4)$criterios)))
  expect_match(a$criterios[1], "PP2")
  expect_match(a$criterios[1], "PP3")
  expect_match(a$criterios[1], "PP4")
  expect_match(a$criterios[2], "BP7")
  expect_match(a$criterios[2], "BP4")
  expect_match(a$criterios[3], "PM4")
  expect_false(grepl("BP7", a$criterios[4]))
})

test_that("PM3 completo solo con fase trans confirmada", {
  x <- acmg_x(consecuencia = rep("missense_variant", 3))
  h <- list(compuesto = c(TRUE, TRUE, FALSE), trans = c(TRUE, FALSE, FALSE))
  a <- acmg(x, h = h)
  expect_match(a$criterios[1], "PM3(,|$)")
  expect_match(a$criterios[2], "PM3_Supporting")
  expect_false(grepl("PM3", a$criterios[3]))
})

test_that("una nula con fenotipo y compuesto trans llega a patogenica", {
  x <- acmg_x(consecuencia = "stop_gained", exon = "2/8")
  a <- acmg(x, h = list(compuesto = TRUE, trans = TRUE), p_fenotipo = 3)
  expect_equal(a$puntos, 8 + 1 + 2 + 1)
  expect_equal(a$clase, "Patog\u00e9nica")
})

test_that("no se sugieren criterios que dependen de la familia", {
  x <- acmg_x(segrega = c(TRUE, FALSE, NA), de_novo_confirmado = c(FALSE, FALSE, TRUE),
              origen = c("materno", "paterno", "de novo"))
  a <- acmg(x, herencia = "AD", h = list(compuesto = rep(FALSE, 3)))
  expect_false(any(grepl("PS2|PM6|PP1|BS4", a$criterios)))
})
