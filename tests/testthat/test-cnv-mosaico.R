vcf_cnv <- function() {
  f <- withr::local_tempfile(fileext = ".vcf", .local_envir = parent.frame())
  writeLines(c(
    "##fileformat=VCFv4.2", "##reference=GRCh38",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Tipo\">",
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Fin\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copias\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "P01", "P02"), collapse = "\t"),
    "chrX\t101360000\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=101370000\tGT:CN\t1:0\t0/0:2",
    "chr22\t19000000\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=21500000\tGT:CN\t0/0:2\t0/1:1",
    "chr5\t1000000\t.\tN\t<DUP>\t50\tPASS\tEND=1005000\tGT:CN\t0/1:3\t0/0:2",
    "chrX\t101365000\t.\tC\tT\t60\tPASS\t.\tGT:AD:DP:GQ\t0/1:10,10:20:99\t0/0:20,0:20:99"
  ), f)
  suppressWarnings(suppressMessages(iei_leer_vcf(f, sexo = c(P01 = "M", P02 = "F"))))
}

test_that("los registros simbolicos se separan como CNV y no como variantes", {
  v <- vcf_cnv()
  expect_equal(nrow(v[v$genotipo != "hom_ref", ]), 1)
  cnv <- attr(v, "ieiprio_cnv")
  expect_equal(nrow(cnv), 3)
  expect_setequal(cnv$tipo, c("DEL", "DUP"))
})

test_that("iei_cnv cruza con el panel y avisa la region 22q11.2", {
  v <- vcf_cnv()
  c <- iei_cnv(v)
  expect_equal(c$genes[c$muestra == "P01" & c$chr == "X"], "BTK")
  expect_match(c$genes[c$chr == "22"], "TBX1")
  expect_match(c$nota[c$chr == "22"], "22q11.2")
  expect_true(is.na(c$genes[c$chr == "5"]))
  expect_equal(c$tamano[c$chr == "X"], 10001L)
  expect_equal(nrow(iei_cnv(v, genes = "BTK")[!is.na(iei_cnv(v, genes = "BTK")$genes), ]), 1)
  expect_message(vacio <- iei_cnv(tibble::tibble(a = 1)), "no trae CNV")
  expect_equal(nrow(vacio), 0)
})

test_that("el reporte muestra los CNV del paciente", {
  local_mocked_bindings(.vep_consultar = vep_falso, .vep_release = function(build) "999")
  v <- vcf_cnv()
  f <- suppressWarnings(suppressMessages(iei_filtrar(v, genes = "BTK", min_gq = 0, min_dp = 0)))
  a <- suppressWarnings(suppressMessages(iei_anotar(f, cache = FALSE)))
  p <- suppressMessages(iei_priorizar(a, fenotipo = "agammaglobulinemia"))
  arch <- withr::local_tempfile(fileext = ".html")
  suppressMessages(iei_reporte(p, arch, paciente = "P01"))
  h <- paste(readLines(arch), collapse = "\n")
  expect_match(h, "CNV y variantes estructurales")
  expect_match(h, "Control de calidad de la muestra")
  expect_match(h, "Entrada \\(VCF\\)")
  expect_match(h, "<b>BTK</b>")
  expect_false(grepl("familiar", h))
  expect_match(h, "puntos\\)")
  expect_match(h, "@page")
})

test_that("en genes de mosaico una fraccion alelica baja se conserva", {
  co <- tibble::tibble(gen = c("FAS", "BTK"), build = "GRCh38", chr = c("10", "X"),
                       inicio = c(1000L, 5000L), fin = c(2000L, 6000L), ensembl_id = NA)
  x <- tibble::tibble(muestra = "P01", chr = c("10", "X"), pos = c(1500L, 5500L), id = NA, ref = "C",
                      alt = "T", genotipo = "het", qual = 99, dp = 100, gq = 99, ab = 0.08,
                      filtro = "PASS", build = "GRCh38")
  attr(x, "build") <- "GRCh38"
  f <- suppressMessages(iei_filtrar(x, coordenadas = co))
  expect_equal(f$gen, "FAS")
  expect_true(f$mosaico_posible)
  f2 <- suppressMessages(iei_filtrar(x, coordenadas = co, min_ab_mosaico = 0.1))
  expect_equal(nrow(f2), 0)
})

test_that("las notas avisan los genes dificiles y los mosaicos", {
  x <- tibble::tibble(muestra = "P01", gen = c("NCF1", "FAS"), chr = c("7", "10"), pos = 1:2,
                      ref = "C", alt = "T", genotipo = "het", impacto = "HIGH",
                      consecuencia = "stop_gained", af_gnomad = NA_real_, clinvar = NA_character_,
                      ab = c(0.5, 0.08), mosaico_posible = c(FALSE, TRUE))
  p <- suppressMessages(iei_priorizar(x))
  expect_match(p$nota[p$gen == "NCF1"], "pseudog")
  expect_match(p$nota[p$gen == "FAS"], "Posible mosaico")
})
