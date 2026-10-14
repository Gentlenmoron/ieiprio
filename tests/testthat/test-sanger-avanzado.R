test_that("un SCF version 3 se escribe y se vuelve a leer igual", {
  a <- cromatograma_sintetico(ref_sintetica(80))
  f <- withr::local_tempfile(fileext = ".scf")
  .scf_escribir(a, f)
  b <- iei_leer_scf(f)
  expect_s3_class(b, "iei_ab1")
  expect_equal(b$bases, a$bases)
  expect_equal(b$picos, a$picos)
  expect_equal(b$calidad, a$calidad)
  expect_equal(as.data.frame(b$trazas), as.data.frame(a$trazas), ignore_attr = TRUE)
  expect_error(iei_leer_scf(fixture("mini.vcf")), "no es un archivo SCF")
  expect_error(iei_leer_scf("no_existe.scf"), "No encuentro")
})

test_that("tambien se lee un SCF version 2", {
  a <- cromatograma_sintetico(ref_sintetica(40))
  f <- withr::local_tempfile(fileext = ".scf")
  n <- nrow(a$trazas); nb <- nchar(a$bases)
  u32 <- function(v) writeBin(as.integer(v), raw(), size = 4, endian = "big")
  m <- as.matrix(a$trazas)
  muestras <- writeBin(as.integer(t(m)), raw(), size = 2, endian = "big")
  vb <- strsplit(a$bases, "")[[1]]
  bases <- do.call(c, lapply(seq_len(nb), function(i) {
    c(u32(a$picos[i] - 1L), as.raw(ifelse(c("A", "C", "G", "T") == vb[i], 30L, 0L)), charToRaw(vb[i]), raw(3))
  }))
  off_b <- 128L + length(muestras)
  cab <- c(charToRaw(".scf"), u32(n), u32(128L), u32(nb), u32(0), u32(0), u32(off_b), u32(0),
           u32(off_b + length(bases)), charToRaw("2.00"), u32(2), u32(0), u32(0), u32(0), raw(72))
  writeBin(c(cab, muestras, bases), f)
  b <- iei_leer_scf(f)
  expect_equal(b$bases, a$bases)
  expect_equal(b$calidad, rep(30L, nb))
  expect_equal(b$trazas$G, a$trazas$G)
})

test_that("los cebadores se recortan al inicio y al final", {
  sec <- ref_sintetica(200, 3)
  cf <- "ACGTTGCAAGGCTTAC"; cr <- "TTGACCGTAGGATCCA"
  lectura <- paste0("NNN", cf, sec, .revcomp(cr), "GGG")
  k <- .recortar_cebadores(lectura, c(cf, cr))
  expect_equal(substr(lectura, min(k), max(k)), sec)
  # Tolera dos diferencias.
  mal <- paste0(substr(cf, 1, 3), "T", substr(cf, 5, 9), "A", substr(cf, 11, 16))
  k2 <- .recortar_cebadores(paste0(mal, sec), cf)
  expect_equal(min(k2), nchar(cf) + 1)
  expect_equal(.recortar_cebadores(sec, NULL), seq_len(nchar(sec)))
})

alelos_indel <- function(ref, desde, hasta, tipo = "del", ins = "GA", corte = 200) {
  a1 <- substr(ref, desde, hasta)
  a2 <- if (tipo == "del") paste0(substr(ref, desde, corte - 1), substr(ref, corte + 3, hasta + 3)) else
    paste0(substr(ref, desde, corte - 1), ins, substr(ref, corte, hasta))
  list(a1 = a1, a2 = a2)
}

esperada <- function(ref, pos, r, a) {
  .normalizar_indels(data.frame(pos = pos, ref = r, alt = a), list(inicio = 1001L, secuencia = ref))
}

sanger_sintetico <- function(lecturas, ref, ...) {
  local_mocked_bindings(.referencia_consultar = function(chr, inicio, fin, build) list(seq = ref),
                        .env = parent.frame())
  archivos <- vapply(seq_along(lecturas), function(i) {
    f <- withr::local_tempfile(fileext = ".scf", .local_envir = parent.frame(3))
    .scf_escribir(lecturas[[i]], f)
    f
  }, "")
  iei_sanger(archivos, region = "X:1001-1400", cache = FALSE, ...)
}

test_that("una delecion heterocigota se llama por deconvolucion en ambas hebras", {
  ref <- ref_sintetica()
  al <- alelos_indel(ref, 51, 330, "del")
  esp <- esperada(ref, 1001L + 198L, substr(ref, 199, 202), substr(ref, 199, 199))
  directa <- cromatograma_sintetico(al$a1, al$a2, archivo = "F")
  v <- suppressWarnings(suppressMessages(sanger_sintetico(list(directa), ref)))
  expect_equal(nrow(v), 1)
  expect_equal(c(v$pos, v$ref, v$alt, v$genotipo), c(esp$pos, esp$ref, esp$alt, "het"))
  expect_match(attr(v, "ieiprio_sanger")[[1]]$indel_llamada, as.character(esp$pos))
  # En la hebra reversa los dos alelos empiezan juntos en el extremo derecho.
  a2r <- paste0(substr(ref, 48, 199), substr(ref, 203, 330))
  reversa <- cromatograma_sintetico(.revcomp(al$a1), .revcomp(a2r), archivo = "R")
  r <- suppressWarnings(suppressMessages(sanger_sintetico(list(reversa), ref)))
  expect_equal(c(r$pos, r$ref, r$alt), c(esp$pos, esp$ref, esp$alt))
})

test_that("una insercion heterocigota se llama con su secuencia", {
  ref <- ref_sintetica()
  al <- alelos_indel(ref, 51, 330, "ins", ins = "GA")
  esp <- esperada(ref, 1001L + 198L, substr(ref, 199, 199), paste0(substr(ref, 199, 199), "GA"))
  d <- cromatograma_sintetico(al$a1, al$a2)
  v <- suppressWarnings(suppressMessages(sanger_sintetico(list(d), ref)))
  expect_equal(c(v$pos, v$ref, v$alt, v$genotipo), c(esp$pos, esp$ref, esp$alt, "het"))
})

test_that("sin deconvolucion queda el aviso para revision manual", {
  ref <- ref_sintetica()
  al <- alelos_indel(ref, 51, 330, "del")
  d <- cromatograma_sintetico(al$a1, al$a2)
  expect_warning(suppressMessages(sanger_sintetico(list(d), ref, indels_het = FALSE)), "picos dobles")
})

test_that("el consenso junta las lecturas y respeta los heterocigotos", {
  ref <- ref_sintetica()
  a1 <- substr(ref, 51, 330)
  a2 <- a1; substr(a2, 100, 100) <- if (substr(a1, 100, 100) == "A") "G" else "A"
  f <- cromatograma_sintetico(a1, a2, alto2 = 0.8, archivo = "F")
  r <- cromatograma_sintetico(.revcomp(a1), .revcomp(a2), alto2 = 0.8, archivo = "R")
  v <- suppressMessages(sanger_sintetico(list(f, r), ref))
  expect_equal(v$pos, 1001L + 149L)
  expect_equal(v$genotipo, "het")
  expect_equal(v$lecturas, 2L)
  arch <- withr::local_tempfile(fileext = ".fasta")
  cons <- suppressMessages(iei_consenso(v, arch))
  expect_true(all(cons$lecturas == 2))
  expect_false(any(cons$discordante))
  s <- attr(cons, "secuencia")
  k <- which(cons$pos == v$pos)
  expect_equal(substr(s, k, k), .iupac_de(v$ref, v$alt))
  substr(s, k, k) <- cons$ref[k]
  expect_equal(s, paste(cons$ref, collapse = ""))
  lineas <- readLines(arch)
  expect_match(lineas[1], "^>sanger X:")
  expect_true(all(nchar(lineas[-1]) <= 60))
  expect_error(iei_consenso(tibble::tibble(a = 1)), "iei_sanger")
})
