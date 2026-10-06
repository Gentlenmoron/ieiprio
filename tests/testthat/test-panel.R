# Estos tests corren sobre los datos incluidos en el paquete.

test_that("el panel tiene la estructura esperada", {
  p <- iuis_panel
  expect_s3_class(p, "tbl_df")
  expect_named(p, c("gen", "hgnc_id", "herencia", "herencia_panelapp",
                    "evidencia", "fenotipos", "categorias_iuis",
                    "ganancia_funcion"))
  expect_false(any(duplicated(p$gen)))
  expect_true(all(p$herencia %in% c("AD", "AR", "AR/AD", "XL", "MT", NA)))
  expect_true(all(p$evidencia %in% c("verde", "ambar", "rojo", NA)))
})

test_that("genes clasicos estan en el panel", {
  expect_true(all(c("BTK", "CYBB", "RAG1", "IL2RG") %in% iuis_panel$gen))
  expect_equal(iuis_panel$herencia[iuis_panel$gen == "BTK"], "XL")
})

test_that("los alias y categorias apuntan a genes del panel", {
  expect_true(all(iuis_alias$gen %in% iuis_panel$gen))
  expect_true(all(iuis_categorias$gen %in% iuis_panel$gen))
  expect_true(all(hpo_genes$gen %in% iuis_panel$gen))
})

test_that("iei_panel filtra por evidencia, categoria y herencia", {
  expect_true(all(iei_panel()$evidencia == "verde"))
  expect_gte(nrow(iei_panel(evidencia = c("verde", "ambar", "rojo"))),
             nrow(iei_panel()))
  expect_true("BTK" %in% iei_panel(categoria = 3)$gen)
  xl <- iei_panel(herencia = "XL")
  expect_true(all(xl$herencia == "XL"))
  ar <- iei_panel(herencia = "AR", evidencia = c("verde", "ambar", "rojo"))
  expect_true(all(grepl("AR", ar$herencia)))
})

test_that("iei_panel da errores claros", {
  expect_error(iei_panel(categoria = 11), "entre 1 y 10")
  expect_error(iei_panel(herencia = "mitocondrial"), "no reconocido")
  expect_error(iei_panel(evidencia = "azul"), "solo acepta")
})

test_that("iei_buscar_gen encuentra por simbolo y por alias", {
  r <- iei_buscar_gen(c("TACI", "btk"))
  expect_equal(r$gen, c("TNFRSF13B", "BTK"))
  expect_equal(r$consulta, c("TACI", "btk"))
})

test_that("genes desconocidos generan aviso y se omiten", {
  expect_warning(r <- iei_buscar_gen(c("BTK", "GENFALSO")), "GENFALSO")
  expect_equal(r$gen, "BTK")
  expect_warning(vacio <- iei_buscar_gen("GENFALSO"))
  expect_equal(nrow(vacio), 0)
})

test_that("iei_fuentes registra versiones", {
  f <- iei_fuentes()
  expect_true(all(c("PanelApp", "HPO") %in% f$fuente))
  expect_false(any(is.na(f$md5)))
})
