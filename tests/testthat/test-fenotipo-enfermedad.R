test_that("la similitud vectorizada coincide con la similitud por pares", {
  o <- .ontologia()
  caso <- c("HP:0006532", "HP:0004432", "HP:0002719")
  otros <- c("HP:0002205", "HP:0004313", "HP:0001903", "HP:0006532", "HP:0000964")
  for (t in caso) {
    v <- .similitud_vector(o$idx[[t]])[o$idx[otros]]
    p <- vapply(otros, function(g) .similitud_terminos(t, g), numeric(1))
    expect_equal(unname(v), unname(p), tolerance = 1e-9)
  }
  expect_equal(sum(.similitud_vector(NA_integer_)), 0)
})

test_that("la similitud por enfermedad nombra la enfermedad mas parecida", {
  r <- .fenotipo_por_gen(c("HP:0004432", "HP:0006532"), c("BTK", "STAT1", "GEN_INVENTADO"), "enfermedad")
  expect_equal(r$gen, c("BTK", "STAT1", "GEN_INVENTADO"))
  expect_match(r$enfermedad[1], "agammaglobulinemia", ignore.case = TRUE)
  expect_gt(r$similitud[1], r$similitud[2])
  expect_equal(r$similitud[3], 0)
  expect_true(is.na(r$enfermedad[3]))
  expect_equal(dim(attr(r, "matriz")), c(2L, 3L))
})

test_that("un signo ausente cambia la enfermedad elegida y baja el puntaje", {
  caso <- c("HP:0001873", "HP:0002719")        # trombocitopenia, infecciones recurrentes
  con <- .fenotipo_por_gen(caso, "WAS", "enfermedad")
  sin_eccema <- .fenotipo_por_gen(caso, "WAS", "enfermedad", ausentes = "HP:0000964")
  expect_lte(sin_eccema$similitud, con$similitud)
  b <- iei_buscar_fenotipo(caso, ausentes = "HP:0000964", n = 400)
  expect_equal(attr(b, "ausentes"), "HP:0000964")
  expect_message(utils::capture.output(print(b)), "Ausentes")
})

test_that("un signo que la enfermedad excluye resta", {
  # La osteopetrosis intermedia (CLCN7) tiene la trombocitopenia como signo excluido.
  r <- .similitud_enfermedades("HP:0001873")
  k <- which(r$info$gen == "CLCN7" & r$info$enfermedad_id == "ORPHA:210110")
  expect_equal(r$matriz[1, k], -1)
  expect_equal(r$puntaje[k], 0)
  e <- hpo_enfermedades[hpo_enfermedades$enfermedad_id == "ORPHA:210110" & hpo_enfermedades$gen == "CLCN7", ]
  expect_match(.mejor_coincidencia("HP:0001873", e$hpo_id[!e$excluido], e$hpo_id[e$excluido]), "excluye")
})

test_that("los tres metodos dan valores entre 0 y 1", {
  for (m in c("enfermedad", "semantica", "exacta")) {
    r <- .fenotipo_por_gen(c("HP:0004432", "HP:0006532"), c("BTK", "IL2RG", "CYBB"), m)
    expect_true(all(r$similitud >= 0 & r$similitud <= 1))
  }
  expect_equal(nrow(.fenotipo_por_gen(character(0), "BTK")), 1)
})

test_that("el valor p empirico esta entre 0 y 1", {
  withr::local_seed(1)
  b <- iei_buscar_fenotipo(c("HP:0004432", "HP:0006532"), n = 5, p_valor = TRUE, permutaciones = 20)
  expect_true(all(b$p_valor > 0 & b$p_valor <= 1))
  expect_lt(b$p_valor[1], 0.2)
})

test_that("la validacion simulada compara los metodos", {
  v <- suppressMessages(iei_validar_fenotipo(n = 8, semilla = 3))
  expect_s3_class(v, "iei_validacion")
  expect_setequal(unique(v$metodo), c("enfermedad", "semantica", "exacta"))
  expect_true(all(v$puesto >= 1))
  r <- .resumen_validacion(v)
  expect_equal(nrow(r), 3)
  expect_true(all(r$top10 >= r$top5 & r$top5 >= r$primero))
  expect_message(utils::capture.output(print(v)), "Validacion")
  # Misma semilla, mismo resultado.
  expect_equal(v$puesto, suppressMessages(iei_validar_fenotipo(n = 8, semilla = 3))$puesto)
})

test_that("iei_priorizar descuenta los signos ausentes", {
  x <- resultado_mini()
  x <- x[x$gen == "BTK", ][1, ]
  a <- suppressMessages(iei_priorizar(x, fenotipo = c("HP:0004432", "HP:0006532")))
  b <- suppressMessages(iei_priorizar(x, fenotipo = c("HP:0004432", "HP:0006532"),
                                      fenotipo_ausente = "HP:0001903"))
  expect_match(a$razones, "agammaglobulinemia", ignore.case = TRUE)
  expect_lte(b$p_fenotipo, a$p_fenotipo)
  expect_equal(attr(b, "ieiprio_fenotipo_ausente"), "HP:0001903")
})
