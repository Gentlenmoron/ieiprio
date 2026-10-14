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

test_that(".construir_panel agrega categorias y mecanismo", {
  curado <- data.frame(gen = c("STAT1", "BTK"), categoria_iuis = c(6, 1))
  mec <- data.frame(gen = c("STAT1", "PIK3CD"), mecanismo = c("LOF y GOF", "LOF y GOF"))
  p <- .construir_panel(pa_mini(), curado, mec)
  expect_equal(p$categorias_iuis[p$gen == "BTK"], "1;3")        # fenotipo + curado
  expect_equal(p$categorias_iuis[p$gen == "RAG1"], "1")         # SCID en el fenotipo
  expect_equal(p$categorias_iuis[p$gen == "TNFRSF13B"], "3")    # CVID
  expect_true(is.na(p$categorias_iuis[p$gen == "GENROJO"]))
  expect_true(p$ganancia_funcion[p$gen == "STAT1"])
  expect_equal(p$mecanismo[p$gen == "PIK3CD"], "LOF y GOF")
  expect_false(p$ganancia_funcion[p$gen == "BTK"])
  # Formato anterior del curado, con columna ganancia_funcion.
  viejo <- .construir_panel(pa_mini(), data.frame(gen = "STAT1", categoria_iuis = 6,
                                                  ganancia_funcion = TRUE))
  expect_equal(viejo$mecanismo[viejo$gen == "STAT1"], "GOF")
})

test_that(".categorias_iuis_derivar registra la fuente de cada categoria", {
  g <- c("A", "B", "C", "D")
  f <- c("Predominantly antibody deficiencies; Agammaglobulinemia",
         "Chronic granulomatous disease", "Fenotipo sin relacion", NA)
  cats <- .categorias_iuis_derivar(g, f, data.frame(gen = "C", categoria_iuis = 8))
  expect_equal(cats$categoria_iuis[cats$gen == "A"], 3L)
  expect_equal(cats$fuente[cats$gen == "A"], "PanelApp, tabla IUIS")
  expect_equal(cats$fuente[cats$gen == "B"], "PanelApp, fenotipo")
  expect_equal(cats$fuente[cats$gen == "C"], "curado")
  expect_false("D" %in% cats$gen)
  expect_true(all(iuis_categorias$gen %in% iuis_panel$gen))
  expect_gt(length(unique(iuis_categorias$gen)), 400)
})

test_that("la URL incluye la version", {
  expect_match(.panelapp_url_panel(398, "8.78"), "panels/398/\\?version=8.78$")
})

test_that("PanelApp responde (en vivo)", {
  saltar_si_no_en_vivo()
  pa <- .panelapp_procesar(en_vivo(.panelapp_descargar(398, "8.78")))
  expect_equal(pa$version, "8.78")
  expect_true("BTK" %in% pa$panel$gen)
  btk <- pa$coordenadas[pa$coordenadas$gen == "BTK", ]
  expect_setequal(btk$build, c("GRCh37", "GRCh38"))
  expect_true(all(btk$chr == "X"))
})

test_that("extrae coordenadas por build y descarta ubicaciones malformadas", {
  co <- pa_mini()$coordenadas
  expect_named(co, c("gen", "build", "chr", "inicio", "fin", "ensembl_id"))
  btk <- co[co$gen == "BTK", ]
  expect_setequal(btk$build, c("GRCh37", "GRCh38"))
  expect_equal(btk$chr, c("X", "X"))
  expect_equal(btk$inicio[btk$build == "GRCh38"], 101349338L)
  expect_false("GENROJO" %in% co$gen)
  expect_false("STAT1" %in% co$gen)
})
