test_that(".panelapp_procesar arma el panel", {
  pa <- pa_mini()
  expect_equal(pa$version, "8.78")
  expect_equal(nrow(pa$panel), 9)
  expect_named(pa$panel, c("gen", "hgnc_id", "herencia", "herencia_panelapp",
                           "evidencia", "fenotipos"))
  btk <- pa$panel[pa$panel$gen == "BTK", ]
  expect_equal(btk$herencia, "XL")
  expect_equal(btk$evidencia, "verde")
  expect_equal(pa$panel$herencia[pa$panel$gen == "GENROJO"], NA_character_)
  expect_true(is.na(pa$panel$fenotipos[pa$panel$gen == "IKZF2"]))
})

test_that("usa gene_symbol si falta hgnc_symbol", {
  expect_true("IL2RG" %in% pa_mini()$panel$gen)
})

test_that("junta alias y simbolos anteriores", {
  al <- pa_mini()$alias
  expect_equal(al$gen[al$alias == "TACI"], "TNFRSF13B")
  expect_equal(al$gen[al$alias == "AGMX1"], "BTK")
  expect_false(any(duplicated(al)))
})

test_that("un panel vacio da error claro", {
  expect_error(.panelapp_procesar(list(genes = list())), "no trae genes")
})

test_that(".construir_panel agrega categorias y ganancia de funcion", {
  curado <- data.frame(gen = c("STAT1", "BTK", "BTK"),
                       categoria_iuis = c(6, 3, 1),
                       ganancia_funcion = c(TRUE, FALSE, FALSE))
  p <- .construir_panel(pa_mini(), curado)
  expect_equal(p$categorias_iuis[p$gen == "BTK"], "1;3")
  expect_true(p$ganancia_funcion[p$gen == "STAT1"])
  expect_true(is.na(p$categorias_iuis[p$gen == "RAG1"]))
})

test_that("la URL incluye la version", {
  expect_match(.panelapp_url_panel(398, "8.78"), "panels/398/\\?version=8.78$")
})

test_that("PanelApp responde (en vivo)", {
  # Solo corre si activas: Sys.setenv(IEIPRIO_EN_VIVO = "true")
  skip_if_not(Sys.getenv("IEIPRIO_EN_VIVO") == "true", "test en vivo desactivado")
  pa <- .panelapp_procesar(.panelapp_descargar(398, "8.78"))
  expect_equal(pa$version, "8.78")
  expect_true("BTK" %in% pa$panel$gen)
})
