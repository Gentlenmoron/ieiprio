test_that("variantes_practica trae las columnas y genes del panel", {
  expect_s3_class(variantes_practica, "tbl_df")
  expect_named(variantes_practica, c("gen", "variante", "proteina", "tipo", "genotipo", "herencia",
                                     "enfermedad", "enfermedad_id", "fenotipo", "que_ensena"))
  expect_true(all(variantes_practica$gen %in% iuis_panel$gen))
  expect_true(all(variantes_practica$genotipo %in% c("het", "hom_alt", "hemi")))
  expect_true(all(startsWith(variantes_practica$variante, paste0(variantes_practica$gen, " c."))))
})

test_that("los signos de cada variante estan anotados en su enfermedad", {
  for (i in seq_len(nrow(variantes_practica))) {
    anot <- hpo_enfermedades[hpo_enfermedades$enfermedad_id == variantes_practica$enfermedad_id[i] &
                               !hpo_enfermedades$excluido, ]
    signos <- strsplit(variantes_practica$fenotipo[i], ";")[[1]]
    expect_length(signos, 3)
    expect_true(all(signos %in% anot$hpo_id), info = variantes_practica$gen[i])
  }
})

test_that("las variantes del banco se reconocen como texto", {
  ext <- suppressMessages(iei_extraer_variantes(variantes_practica$variante))
  expect_equal(nrow(ext), nrow(variantes_practica))
  expect_equal(ext$entrada, sub(" ", ":", variantes_practica$variante, fixed = TRUE))
  expect_true(all(ext$tipo == "gen_y_cambio"))
})

test_that("el fenotipo de cada variante pone a su gen entre los primeros", {
  for (g in c("BTK", "STAT3", "AIRE")) {
    v <- variantes_practica[variantes_practica$gen == g, ]
    b <- suppressMessages(iei_buscar_fenotipo(strsplit(v$fenotipo, ";")[[1]], n = 25))
    expect_true(g %in% b$gen, info = g)
  }
})
