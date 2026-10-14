recoder_falso <- function(entradas, build) {
  r <- jsonlite::read_json(fixture("recoder_sintetico.json"))
  Filter(function(el) {
    b <- if (!is.null(el$input)) el else el[[1]]
    b$input %in% entradas
  }, r)
}

test_that("extrae HGVS, rsID y gen con cambio, con genotipo", {
  v <- iei_extraer_variantes(c(
    "Se identifico la variante NM_000061.3(BTK):c.1632-2A>G en hemicigosis.",
    "Tambien rs121908711 y STAT1 p.(Arg274Gln) en heterocigosis.",
    "Genomica NC_000023.11:g.101353990T>C",
    "Sin variantes aqui. El paciente tiene IgG baja."
  ))
  expect_equal(v$entrada, c("NM_000061.3:c.1632-2A>G", "rs121908711", "STAT1:p.Arg274Gln",
                            "NC_000023.11:g.101353990T>C"))
  expect_equal(v$tipo, c("hgvs_transcrito", "rsid", "gen_y_cambio", "hgvs_genomico"))
  expect_equal(v$genotipo, c("hemi", "het", "het", NA))
  expect_equal(v$linea, c(1L, 2L, 2L, 3L))
})

test_that("no duplica la variante escrita con el gen entre parentesis", {
  v <- iei_extraer_variantes("NM_000061.3(BTK):c.1632-2A>G")
  expect_equal(nrow(v), 1)
})

test_that("normaliza guiones tipograficos y signos raros", {
  v <- iei_extraer_variantes("BTK c.1632\u22122A&gt;G homocigoto")
  expect_equal(v$entrada, "BTK:c.1632-2A>G")
  expect_equal(v$genotipo, "hom_alt")
})

test_that("texto sin variantes devuelve tabla vacia", {
  v <- iei_extraer_variantes(c("nada", NA, ""))
  expect_equal(nrow(v), 0)
  expect_named(v, c("entrada", "tipo", "genotipo", "linea"))
})

test_that("lee Word, PDF y texto plano", {
  w <- iei_leer_texto(fixture("informe.docx"))
  expect_true(any(grepl("c.1632-2A>G", w, fixed = TRUE)))
  vw <- iei_extraer_variantes(w)
  expect_equal(vw$entrada, c("NM_000061.3:c.1632-2A>G", "RAG1:p.Arg396Cys"))
  expect_equal(vw$genotipo, c("hemi", "het"))

  skip_if_not_installed("pdftools")
  vp <- iei_extraer_variantes(iei_leer_texto(fixture("informe.pdf")))
  expect_equal(vp$entrada, c("BTK:c.1632-2A>G", "rs999"))

  tmp <- withr::local_tempfile(fileext = ".txt")
  writeLines("rs999", tmp)
  expect_equal(iei_leer_texto(tmp), "rs999")
  expect_error(iei_leer_texto("no_existe.pdf"), "No encuentro")
})

test_that(".recoder_aplanar entiende la respuesta de Ensembl", {
  r <- .recoder_aplanar(jsonlite::read_json(fixture("recoder_sintetico.json")))
  expect_equal(nrow(r), 4)
  btk <- r[r$entrada == "NM_000061.3:c.1632-2A>G", ]
  expect_equal(nrow(btk), 1)
  expect_equal(btk$chr, "X"); expect_equal(btk$pos, 101353990L)
  expect_false(btk$ambigua)
  expect_false(any(grepl("LRG", r$chr)))
  expect_true(all(r$ambigua[r$entrada == "RAG1:p.Arg396Cys"]))
  expect_equal(nrow(.recoder_aplanar(list(error = "x"))), 0)
})

test_that("iei_variantes consulta cada variante por separado", {
  llamadas <- list()
  local_mocked_bindings(.recoder_consultar = function(entradas, build) {
    llamadas[[length(llamadas) + 1]] <<- entradas
    recoder_falso(entradas, build)
  })
  suppressMessages(iei_variantes(c("NM_000061.3:c.1632-2A>G", "rs999")))
  expect_equal(lengths(llamadas), c(1L, 1L))
})

test_that("iei_variantes arma una tabla compatible con el flujo", {
  local_mocked_bindings(.recoder_consultar = recoder_falso)
  co <- tibble::tibble(gen = c("BTK", "RAG1"), build = "GRCh38", chr = c("X", "11"),
                       inicio = c(101349338L, 36510709L), fin = c(101390796L, 36593156L),
                       ensembl_id = NA_character_)
  local_mocked_bindings(.genes_de_variantes = function(tabla, build, margen = 20) {
    tabla$gen <- vapply(seq_len(nrow(tabla)), function(i) {
      k <- which(co$chr == tabla$chr[i] & tabla$pos[i] >= co$inicio & tabla$pos[i] <= co$fin)
      if (length(k) == 0) NA_character_ else co$gen[k[1]]
    }, character(1))
    tabla[, c("muestra", "gen", setdiff(names(tabla), c("muestra", "gen")))]
  })
  expect_warning(
    v <- suppressMessages(iei_variantes(c("NM_000061.3:c.1632-2A>G homocigoto",
                                           "RAG1 p.Arg396Cys", "BTK c.999999X>Y"), sexo = "M")),
    "no reconocio")
  expect_equal(v$gen[v$entrada == "NM_000061.3:c.1632-2A>G"], "BTK")
  expect_equal(v$genotipo[v$entrada == "NM_000061.3:c.1632-2A>G"], "hemi")
  expect_equal(sum(v$entrada == "RAG1:p.Arg396Cys"), 2)
  expect_true(all(c("muestra", "chr", "pos", "ref", "alt", "genotipo", "build") %in% names(v)))
  expect_equal(attr(v, "build"), "GRCh38")
  expect_equal(iei_bitacora(v)$filas, c(3L, 3L))

  g <- suppressMessages(iei_variantes("rs999", genotipo = "het"))
  expect_equal(g$genotipo, "het")
  expect_true(is.na(g$gen))
})

test_that("variantes sin genotipo se priorizan con una nota", {
  x <- tibble::tibble(muestra = "consulta", gen = "BTK", genotipo = NA_character_,
                      impacto = "HIGH", consecuencia = "splice_acceptor_variant",
                      af_gnomad = NA_real_, clinvar = "pathogenic")
  p <- iei_priorizar(x)
  expect_equal(p$p_herencia, 0)
  expect_match(p$nota, "Genotipo no indicado")
})

test_that("iei_variantes valida la entrada", {
  expect_error(iei_variantes(123), "texto")
  expect_error(iei_variantes("aqui no hay nada"), "No encontre")
})

test_that("Variant Recoder en vivo", {
  saltar_si_no_en_vivo()
  btk <- en_vivo(.recoder_aplanar(.recoder_consultar("NM_000061.3:c.1632-2A>G", "GRCh38")))
  expect_equal(btk$chr, "X")
  expect_equal(btk$pos, 101353990L)
  expect_false(btk$ambigua)
  rs <- en_vivo(.recoder_aplanar(.recoder_consultar("rs699", "GRCh38")))
  expect_equal(rs$chr, "1")
  expect_equal(rs$pos, 230710048L)
})

test_that("al anotar se resuelven las posiciones ambiguas por la nomenclatura", {
  x <- tibble::tibble(
    muestra = "consulta", gen = "BTK", chr = "X", pos = c(1L, 2L, 3L), ref = "A", alt = "G",
    genotipo = "hemi", entrada = c(rep("BTK:c.1632-2A>G", 3)), ambigua = TRUE)
  attr(x, "build") <- "GRCh38"
  attr(x, "ieiprio_bitacora") <- tibble::tibble(paso = "Ubicadas", filas = 3L, descartadas = 0L)
  x$hgvs_c <- c("ENST1:c.1632-2A>G", "ENST1:c.1630A>G", "ENST1:c.1729A>G")
  x$hgvs_p <- c(NA, "ENSP1:p.Arg544Gly", "ENSP1:p.Lys577Glu")
  r <- .resolver_ambiguas(x)
  expect_equal(r$pos, 1L)
  expect_false(r$ambigua)
  b <- iei_bitacora(r)
  expect_equal(b$paso[nrow(b)], "Posiciones ambiguas resueltas")
  expect_equal(b$descartadas[nrow(b)], 2L)

  # Por proteina
  y <- x; y$entrada <- "BTK:p.Lys577Glu"
  expect_equal(.resolver_ambiguas(y)$pos, 3L)

  # Si nada coincide, se conservan todas
  z <- x; z$entrada <- "BTK:c.999A>G"
  expect_equal(nrow(.resolver_ambiguas(z)), 3)

  # Tablas sin estas columnas no se tocan
  expect_identical(.resolver_ambiguas(x[, c("muestra", "gen")]), x[, c("muestra", "gen")])
})

test_that("la ficha resume los fenotipos repetidos de PanelApp", {
  f <- "Agammaglobulinemia, X-linked; Agammaglobulinemia, X-linked 1, 300755; agammaglobulinemia, x-linked; CVID; Otro; Mas"
  expect_equal(.resumir_fenotipos(f), "Agammaglobulinemia, X-linked; Agammaglobulinemia, X-linked 1; CVID (y 2 mas)")
  expect_equal(.resumir_fenotipos(NA), "-")
})
