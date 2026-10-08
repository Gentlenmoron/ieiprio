test_that("iei_analizar corre el flujo completo y genera reportes", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999",
                        .proteina_consultar = function(id, build) NULL)
  withr::local_options(list(rlib_message_verbosity = "quiet"))
  dir <- withr::local_tempdir()
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  p <- iei_analizar(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"),
                    fenotipo = "HP:0004313", carpeta = dir, responsable = "Prueba")
  expect_true(all(c("puntaje", "categoria", "acmg") %in% names(p)))
  expect_setequal(basename(attr(p, "reportes")), c("P01.html", "P02.html"))
  expect_true(all(file.exists(attr(p, "reportes"))))
  expect_equal(attr(p, "ieiprio_fenotipo"), "HP:0004313")
})

test_that("iei_analizar puede omitir reportes y pasar argumentos al filtro", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999")
  withr::local_options(list(rlib_message_verbosity = "quiet"))
  withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
  p <- iei_analizar(fixture("mini.vcf"), sexo = "M", carpeta = NULL, min_dp = 50)
  expect_null(attr(p, "reportes"))
  expect_equal(attr(p, "ieiprio_parametros")$min_dp, 50)
})
