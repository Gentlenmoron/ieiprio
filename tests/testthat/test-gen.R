test_that("la ficha local se arma sin conexion", {
  g <- iei_gen("TACI", clinvar = FALSE, proteina = FALSE)
  expect_s3_class(g, "iei_gen")
  expect_equal(g$gen, "TNFRSF13B")
  expect_true("TACI" %in% g$alias)
  expect_equal(g$categorias$categoria_iuis, 3L)
  expect_gt(nrow(g$hpo), 10)
  expect_true(all(diff(g$hpo$ic) <= 0))         # ordenados de mas a menos especifico
  expect_null(g$clinvar)
  out <- paste(cli::cli_fmt(print(g)), collapse = " ")
  expect_match(out, "Categor.a IUIS: 3")
  s <- iei_gen("STAT3", clinvar = FALSE, proteina = FALSE)
  expect_equal(s$mecanismo, "LOF y GOF")
  expect_error(iei_gen(c("BTK", "RAG1")), "un solo gen")
  expect_warning(iei_gen("NOEXISTE1", clinvar = FALSE, proteina = FALSE), "no est. en el panel")
})

test_that("con ClinVar y Ensembl la ficha cuenta y ubica las variantes", {
  withr::local_options(ieiprio.sin_red = FALSE)
  local_mocked_bindings(.clinvar_peticion = clinvar_falso,
                        .proteina_canonica = function(gen, build) "ENSP00000308176",
                        .proteina_consultar = proteina_falsa)
  g <- iei_gen("BTK", cache = FALSE)
  expect_equal(g$clinvar$conteo$variantes, c(120L, 40L, 300L, 90L, 60L))
  v <- g$clinvar$variantes
  expect_equal(nrow(v), 4)                     # la de otro gen se descarta
  expect_equal(v$clinvar_id[1], 103L)          # panel de expertos primero
  expect_equal(v$estrellas[v$clinvar_id == 103], 3L)
  expect_equal(v$pos_proteina[v$clinvar_id == 101], 544L)
  expect_equal(v$tipo[v$clinvar_id == 104], "sin sentido")
  expect_equal(v$tipo[v$clinvar_id == 102], "splicing")
  expect_equal(g$proteina$largo, 659)
  out <- paste(cli::cli_fmt(print(g)), collapse = " ")
  expect_match(out, "patogenica 120")
  expect_match(out, "SH3 1")
})

test_that("si el servicio falla, avisa y la ficha sale igual; la cache evita repetir", {
  withr::local_options(ieiprio.sin_red = FALSE)
  local_mocked_bindings(.clinvar_peticion = function(...) stop("sin red"),
                        .proteina_canonica = function(...) stop("sin red"))
  expect_warning(g <- iei_gen("BTK", cache = FALSE, proteina = FALSE), "en l")
  expect_null(g$clinvar)
  expect_equal(g$herencia, "XL")
  d <- withr::local_tempdir()
  llamadas <- 0
  local_mocked_bindings(.clinvar_peticion = function(...) { llamadas <<- llamadas + 1; clinvar_falso(...) })
  iei_gen("BTK", proteina = FALSE, cache = d)
  n1 <- llamadas
  iei_gen("BTK", proteina = FALSE, cache = d)
  expect_equal(llamadas, n1)
})

test_that("ClinVar en vivo para BTK", {
  saltar_si_no_en_vivo()
  g <- iei_gen("BTK", cache = FALSE)
  skip_if(is.null(g$clinvar), "ClinVar no respondio")
  expect_gt(g$clinvar$conteo$variantes[1], 50)
  expect_gt(nrow(g$clinvar$variantes), 50)
  expect_gt(mean(!is.na(g$clinvar$variantes$pos_proteina)), 0.5)
  expect_equal(g$proteina$largo, 659)
})
