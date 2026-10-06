test_that("el panel tiene la estructura esperada", {
  p <- iei_panel()
  expect_s3_class(p, "tbl_df")
  expect_named(p, c("gen", "categoria_iuis", "categoria_nombre",
                    "enfermedad", "herencia", "ganancia_funcion"))
  expect_true(all(p$categoria_iuis %in% 1:10))
  expect_true(all(p$herencia %in% c("AR", "AD", "XL", "XLR", "AR/AD")))
})

test_that("no hay filas duplicadas de gen con la misma herencia y enfermedad", {
  p <- iei_panel()
  expect_false(any(duplicated(p[, c("gen", "herencia", "enfermedad")])))
})

test_that("todos los alias apuntan a genes del panel", {
  expect_true(all(iuis_alias$gen %in% iuis_panel$gen))
})

test_that("iei_panel filtra por categoria", {
  p <- iei_panel(categoria = 3)
  expect_true(all(p$categoria_iuis == 3))
  expect_true("BTK" %in% p$gen)
})

test_that("iei_panel filtra por herencia y respeta AR/AD", {
  ar <- iei_panel(herencia = "AR")
  expect_true("TNFRSF13B" %in% ar$gen)
  expect_false("BTK" %in% ar$gen)
  expect_true(all(iei_panel(herencia = "xl")$herencia == "XL"))
})

test_that("iei_panel da errores claros", {
  expect_error(iei_panel(categoria = 11), "entre 1 y 10")
  expect_error(iei_panel(herencia = "mitocondrial"), "no reconocido")
})

test_that("iei_buscar_gen encuentra por simbolo y por alias", {
  r <- iei_buscar_gen("TACI")
  expect_equal(unique(r$gen), "TNFRSF13B")
  expect_equal(unique(r$consulta), "TACI")
  expect_equal(unique(iei_buscar_gen("btk")$gen), "BTK")
})

test_that("un gen con varias enfermedades devuelve varias filas", {
  r <- iei_buscar_gen("STAT1")
  expect_equal(nrow(r), 2)
  expect_setequal(r$herencia, c("AD", "AR"))
})

test_that("genes desconocidos generan aviso y se omiten", {
  expect_warning(r <- iei_buscar_gen(c("BTK", "GENFALSO")), "GENFALSO")
  expect_equal(unique(r$gen), "BTK")
  expect_warning(vacio <- iei_buscar_gen("GENFALSO"))
  expect_equal(nrow(vacio), 0)
  expect_true("consulta" %in% names(vacio))
})

test_that("iei_buscar_gen valida la entrada", {
  expect_error(iei_buscar_gen(123), "vector de texto")
})
