test_that(".obo_procesar lee terminos, padres, sinonimos exactos y alternos", {
  t <- .obo_procesar(fixture("mini.obo"))
  expect_equal(t$id, c("HP:0000001", "HP:0000002", "HP:0000003", "HP:0000004"))
  expect_equal(t$padres[[4]], "HP:0000003")
  expect_equal(t$sinonimos[[3]], "Infection")
  expect_equal(t$alternos[[4]], "HP:0000099")
})

test_that("los ancestros incluyen al termino y suben hasta la raiz", {
  t <- .obo_procesar(fixture("mini.obo"))
  a <- .ancestros_todos(t$id, t$padres)
  expect_setequal(t$id[a[[4]]], t$id)
  expect_equal(t$id[a[[1]]], "HP:0000001")
})

test_that("el contenido de informacion crece con la especificidad", {
  t <- .obo_procesar(fixture("mini.obo"))
  ic <- .ic_calcular(t, fixture("mini_g2p.txt"))
  expect_equal(ic[1], 0)                       # todos los genes tienen la raiz
  expect_true(all(diff(ic) >= 0))
  expect_equal(ic[4], -log(2 / 4), tolerance = 1e-4)
})

test_that("hpo_terminos cubre la ontologia completa con IC", {
  expect_gt(nrow(hpo_terminos), 15000)
  expect_true(all(hpo_genes$hpo_id %in% hpo_terminos$id))
  expect_equal(hpo_terminos$ic[hpo_terminos$id == "HP:0000001"], 0)
})

test_that("la similitud vale 1 para terminos iguales o mas especificos", {
  # Neumonia recurrente frente a un gen con neumonia recurrente
  expect_equal(.similitud_terminos("HP:0006532", "HP:0006532"), 1)
  # El caso es general (Neumonia) y el gen tiene algo mas especifico: tambien 1
  expect_equal(.similitud_terminos("HP:0002090", "HP:0006532"), 1)
  # Emparentados: entre 0 y 1
  s <- .similitud_terminos("HP:0006532", "HP:0002205")
  expect_gt(s, 0); expect_lt(s, 1)
  # Sin gen, 0
  expect_equal(.similitud_terminos("HP:0006532", character(0)), 0)
})

test_that("el texto se traduce a HPO en espanol y en ingles", {
  i <- .interpretar_fenotipo(c("HP:0004313", "neumon\u00eda recurrente", "Agammaglobulinemia",
                               "chronic diarrhea", "zzqx"))
  expect_equal(i$hpo_id[1:4], c("HP:0004313", "HP:0006532", "HP:0004432", "HP:0002028"))
  expect_equal(i$como, c("codigo", "diccionario", "diccionario", "nombre exacto", "no encontrado"))
  expect_error(.fenotipo_codigos("zzqx"), "interpretar")
  expect_message(.fenotipo_codigos("diarrea"), "HP:0002014")
  expect_null(.fenotipo_codigos(NULL))
})

test_that("los codigos antiguos pasan a su termino vigente", {
  t <- hpo_terminos
  k <- which(lengths(t$alternos) > 0)[1]
  viejo <- t$alternos[[k]][1]
  expect_equal(.hpo_vigente(viejo), t$id[k])
  expect_true(is.na(.hpo_vigente("HP:9999999")))
})

test_that("iei_buscar_hpo encuentra terminos en espanol e ingles", {
  n <- iei_buscar_hpo("neumonia")
  expect_equal(n$hpo_id[1], "HP:0006532")       # del diccionario, primero
  expect_true("HP:0002090" %in% n$hpo_id)
  expect_named(n, c("hpo_id", "termino", "genes_panel", "ic"))
  expect_true(any(grepl("andidiasis", iei_buscar_hpo("candidiasis")$termino)))
  expect_message(vacio <- iei_buscar_hpo("zzqxw"), "No encontr.")
  expect_equal(nrow(vacio), 0)
  expect_error(iei_buscar_hpo(""), "palabra")
})
