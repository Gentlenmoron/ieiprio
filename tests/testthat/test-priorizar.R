panel_prueba <- tibble::tibble(
  gen = c("GXL", "GAR", "GAD", "GGOF", "GSIN"),
  herencia = c("XL", "AR", "AD", "AD", NA),
  ganancia_funcion = c(FALSE, FALSE, FALSE, TRUE, FALSE)
)
hpo_prueba <- tibble::tibble(gen = c("GAD", "GAD", "GXL"),
                             hpo_id = c("HP:0002719", "HP:0004313", "HP:0002719"),
                             hpo_nombre = "x")
cat_prueba <- tibble::tibble(gen = "GAD", categoria_iuis = 3L, categoria_nombre = "x")

variantes <- function() {
  x <- tibble::tribble(
    ~muestra, ~gen,   ~pos, ~genotipo, ~impacto,   ~consecuencia,       ~af_gnomad, ~clinvar,
    "P01",    "GXL",  1,    "hemi",    "HIGH",     "stop_gained",       NA,         "pathogenic",
    "P02",    "GXL",  1,    "het",     "HIGH",     "stop_gained",       NA,         "pathogenic",
    "P01",    "GAR",  2,    "het",     "HIGH",     "frameshift_variant", 1e-5,      "likely_pathogenic",
    "P03",    "GAR",  3,    "het",     "MODERATE", "missense_variant",  NA,         "uncertain_significance",
    "P03",    "GAR",  4,    "het",     "HIGH",     "stop_gained",       NA,         NA,
    "P01",    "GAD",  5,    "het",     "MODERATE", "missense_variant",  NA,         NA,
    "P01",    "GGOF", 6,    "het",     "MODERATE", "missense_variant",  NA,         NA,
    "P01",    "GAD",  7,    "het",     "LOW",      "synonymous_variant", 0.001,     "benign;likely_benign",
    "P04",    "GAD",  8,    "het",     "MODERATE", "missense_variant",  2e-5,      "conflicting_interpretations_of_pathogenicity",
    "P04",    "GSIN", 9,    "hom_alt", "MODERATE", "missense_variant",  NA,         NA
  )
  attr(x, "ieiprio_bitacora") <- tibble::tibble(paso = "Inicial", filas = 10L, descartadas = 0L)
  attr(x, "build") <- "GRCh38"
  x
}
priorizar <- function(...) {
  iei_priorizar(variantes(), panel = panel_prueba, hpo = hpo_prueba,
                categorias = cat_prueba, ...)
}
fila <- function(p, muestra, pos) p[p$muestra == muestra & p$pos == pos, ]

test_that("hemicigota patogenica con impacto alto queda en categoria alta", {
  r <- fila(priorizar(), "P01", 1)
  expect_equal(r$puntaje, 10)
  expect_equal(r$categoria, "alta")
  expect_match(r$razones, "ClinVar patog\u00e9nica \\(\\+4\\)")
  expect_match(r$razones, "hemicigota en gen ligado al X \\(\\+2\\)")
})

test_that("heterocigota en gen XL no suma herencia y queda anotada", {
  r <- fila(priorizar(), "P02", 1)
  expect_equal(r$puntaje, 8)
  expect_equal(r$categoria, "media")
  expect_match(r$nota, "portadora")
})

test_that("heterocigoto unico en gen AR resta y se marca como portador", {
  r <- fila(priorizar(), "P01", 2)
  expect_equal(r$puntaje, 4 + 3 - 2)
  expect_match(r$razones, "heterocigoto \u00fanico en gen recesivo \\(-2\\)")
  expect_match(r$nota, "Portador")
})

test_that("dos heterocigotos en gen AR se marcan como posible compuesto", {
  p <- priorizar()
  expect_equal(fila(p, "P03", 3)$puntaje, 1 + 1 + 2)
  expect_equal(fila(p, "P03", 4)$puntaje, 3 + 1 + 2)
  expect_true(all(grepl("compuesto", p$nota[p$muestra == "P03"])))
})

test_that("missense en gen de ganancia de funcion cuenta como impacto alto", {
  r <- fila(priorizar(), "P01", 6)
  expect_equal(r$puntaje, 3 + 1 + 2)
  expect_match(r$razones, "ganancia de funci.n")
})

test_that("descarta benignas y lo registra en la bitacora", {
  p <- priorizar()
  expect_false(7 %in% p$pos)
  b <- iei_bitacora(p)
  expect_equal(b$paso[nrow(b)], "ClinVar benigna")
  expect_equal(b$descartadas[nrow(b)], 1L)
})

test_that("ClinVar en conflicto suma poco", {
  r <- fila(priorizar(), "P04", 8)
  expect_equal(r$puntaje, 1 + 1 + 2)
  expect_match(r$razones, "conflicto")
})

test_that("gen sin herencia conocida no suma herencia y lo avisa", {
  r <- fila(priorizar(), "P04", 9)
  expect_equal(r$puntaje, 1 + 1)
  expect_match(r$nota, "sin modo de herencia")
})

test_that("el fenotipo HPO suma de forma proporcional", {
  sin <- fila(priorizar(), "P01", 5)
  expect_equal(sin$puntaje, 4)
  expect_equal(sin$categoria, "baja")
  todo <- fila(priorizar(fenotipo = c("HP:0002719", "HP:0004313"), similitud = "exacta"), "P01", 5)
  expect_equal(todo$puntaje, 7)
  expect_match(todo$razones, "fenotipo 2 de 2 t.rminos HPO")
  mitad <- fila(priorizar(fenotipo = c("HP:0002719", "HP:0001903"), similitud = "exacta"), "P01", 5)
  expect_equal(mitad$puntaje, 5.5)
})

test_that("un portador en gen AR nunca queda en prioridad alta", {
  cats <- tibble::tibble(gen = "GAR", categoria_iuis = 1L, categoria_nombre = "x")
  p <- iei_priorizar(variantes(), categoria = 1, pesos = list(portador_ar = 0),
                     panel = panel_prueba, hpo = hpo_prueba, categorias = cats)
  r <- fila(p, "P01", 2)
  expect_gte(r$puntaje, 9)
  expect_equal(r$categoria, "media")
  expect_match(r$razones, "limitada a media")
})

test_that("la categoria IUIS sospechada suma el maximo", {
  r <- fila(priorizar(categoria = 3), "P01", 5)
  expect_equal(r$puntaje, 7)
  expect_match(r$razones, "categor.a IUIS")
})

test_that("los pesos se pueden ajustar", {
  r <- fila(priorizar(pesos = list(clinvar_patogenica = 6)), "P01", 1)
  expect_equal(r$puntaje, 12)
  r2 <- fila(priorizar(pesos = list(corte_alta = 11)), "P01", 1)
  expect_equal(r2$categoria, "media")
})

test_that("ordena por muestra y puntaje y conserva atributos", {
  p <- priorizar()
  expect_equal(p$muestra, sort(p$muestra))
  p01 <- p[p$muestra == "P01", ]
  expect_equal(p01$puntaje, sort(p01$puntaje, decreasing = TRUE))
  expect_equal(attr(p, "build"), "GRCh38")
})

test_that("valida entradas y maneja tablas vacias", {
  expect_error(iei_priorizar(data.frame(a = 1)), "faltan columnas")
  expect_error(priorizar(fenotipo = "zzqx palabra inexistente"), "interpretar")
  vacio <- iei_priorizar(variantes()[0, ], panel = panel_prueba)
  expect_equal(nrow(vacio), 0)
  expect_true("puntaje" %in% names(vacio))
})

test_that(".clasificar_clinvar entiende los formatos de VEP", {
  expect_equal(
    .clasificar_clinvar(c("pathogenic", "likely_pathogenic;pathogenic", "benign",
                          "likely_benign;benign", "pathogenic;benign",
                          "conflicting_interpretations_of_pathogenicity",
                          "uncertain_significance", NA, "pathogenic/likely_pathogenic")),
    c("patogenica", "patogenica", "benigna", "benigna", "conflicto",
      "conflicto", "otra", NA, "patogenica")
  )
})

test_that("flujo completo con datos del paquete", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999")
  v <- iei_leer_vcf(fixture("mini.vcf"), sexo = c(P01 = "M", P02 = "F"))
  f <- iei_filtrar(v, genes = "BTK", coordenadas = coords_prueba())
  a <- suppressMessages(iei_anotar(f, cache = FALSE))
  p <- iei_priorizar(a)
  r <- p[p$muestra == "P01" & p$pos == 101360000, ]
  # patogenica (+4), missense (+1), predictores (+1), hemicigota (+2); en gnomAD (0)
  expect_equal(r$puntaje, 8)
  expect_equal(r$categoria, "media")
  expect_match(r$razones, "REVEL 0.82")
  expect_match(r$razones, "hemicigota")
  expect_match(p$nota[p$muestra == "P02"], "portadora")
})

test_that("la similitud semantica premia terminos emparentados", {
  hpo <- tibble::tibble(gen = c("GAD", "GXL"),
                        hpo_id = c("HP:0002205", "HP:0001875"),   # infecciones respiratorias recurrentes; neutropenia
                        hpo_nombre = "x")
  p <- iei_priorizar(variantes(), fenotipo = "HP:0006532",      # neumonia recurrente
                     panel = panel_prueba, hpo = hpo, categorias = cat_prueba)
  gad <- fila(p, "P01", 5); gxl <- fila(p, "P01", 1)
  expect_gt(gad$p_fenotipo, 0)
  expect_lt(gad$p_fenotipo, 3)
  expect_match(gad$razones, "similitud")
  expect_gt(gad$p_fenotipo, gxl$p_fenotipo)
  ex <- iei_priorizar(variantes(), fenotipo = "HP:0006532", similitud = "exacta",
                      panel = panel_prueba, hpo = hpo, categorias = cat_prueba)
  expect_equal(fila(ex, "P01", 5)$p_fenotipo, 0)
  # Un termino identico vale el maximo tambien con similitud semantica.
  ig <- iei_priorizar(variantes(), fenotipo = "HP:0002205",
                      panel = panel_prueba, hpo = hpo, categorias = cat_prueba)
  expect_equal(fila(ig, "P01", 5)$p_fenotipo, 3)
})

test_that("el fenotipo se puede escribir como texto", {
  hpo <- tibble::tibble(gen = "GAD", hpo_id = "HP:0004432", hpo_nombre = "x")
  expect_message(
    p <- iei_priorizar(variantes(), fenotipo = "Agammaglobulinemia",
                       panel = panel_prueba, hpo = hpo, categorias = cat_prueba),
    "HP:0004432")
  expect_equal(fila(p, "P01", 5)$p_fenotipo, 3)
  expect_equal(attr(p, "ieiprio_fenotipo"), "HP:0004432")
})

test_that("PVS1 solo se omite en genes de ganancia de funcion pura", {
  panel <- panel_prueba
  panel$mecanismo <- c(NA, NA, "LOF y GOF", "GOF", NA)
  v <- variantes()
  v$consecuencia[v$gen == "GGOF"] <- "stop_gained"; v$impacto[v$gen == "GGOF"] <- "HIGH"
  v$consecuencia[v$pos == 5] <- "stop_gained"; v$impacto[v$pos == 5] <- "HIGH"
  p <- iei_priorizar(v, panel = panel, hpo = hpo_prueba, categorias = cat_prueba)
  expect_false(grepl("PVS1", fila(p, "P01", 6)$acmg))
  expect_match(fila(p, "P01", 5)$acmg, "PVS1")
})

test_that("la fase del VCF distingue compuestos en trans de pares en cis", {
  x <- tibble::tibble(muestra = "P01", gen = "RAG1", gt = c("0|1", "1|0", "0|1", "0|1"),
                      ps = c("100", "100", "200", "200"))
  x$gen[3:4] <- "RAG2"
  f <- .fase_compuestos(x, rep(TRUE, 4))
  expect_equal(f, c("trans", "trans", "cis", "cis"))
  sin_fase <- x
  sin_fase$gt <- "0/1"
  expect_equal(.fase_compuestos(sin_fase, rep(TRUE, 4)), rep("desconocida", 4))
})
