# Genera inst/extdata/ejemplo.vcf.gz, un paciente ficticio para la vineta,
# los ejemplos y los tests.
# Correr desde la raiz del paquete:  source("data-raw/ejemplo_vcf.R")
#
# El VCF mezcla:
#   - 5 variantes REALES de ClinVar en genes del panel (coordenadas y alelos
#     correctos en GRCh38, para que Ensembl VEP las anote).
#   - Unas 4000 variantes de fondo SINTETICAS, ubicadas fuera de los genes del
#     panel para que nunca lleguen a VEP. Solo sirven para el control de calidad
#     (histogramas, Ti/Tv, verificacion de sexo).
#
# La primera vez consulta ClinVar y guarda la seleccion en
# data-raw/ejemplo_clinvar.csv. Despues usa ese archivo, asi el VCF se
# reconstruye siempre igual. Para elegir variantes nuevas, borra el CSV.

devtools::load_all()
set.seed(20261007)

diseno <- data.frame(
  rol           = c("causal", "portador", "vus", "benigna", "benigna"),
  gen           = c("BTK", "RAG1", "STAT1", "BTK", "CYBB"),
  clasificacion = c("Pathogenic", "Pathogenic", "Uncertain significance", "Benign", "Benign"),
  genotipo      = c("1/1", "0/1", "0/1", "1/1", "1/1"),  # paciente varon
  stringsAsFactors = FALSE
)

ruta_csv <- "data-raw/ejemplo_clinvar.csv"
if (file.exists(ruta_csv)) {
  clinvar <- utils::read.csv(ruta_csv, stringsAsFactors = FALSE)
  message("Usando la seleccion fijada en ", ruta_csv)
} else {
  message("Consultando ClinVar...")
  usados <- integer(0)
  filas <- list()
  for (i in seq_len(nrow(diseno))) {
    v <- .clinvar_elegir(diseno$gen[i], diseno$clasificacion[i], excluir = usados)
    usados <- c(usados, v$clinvar_id)
    filas[[i]] <- cbind(diseno[i, c("rol", "genotipo")], as.data.frame(v))
    message("  ", diseno$rol[i], ": ", v$gen, " ClinVar ", v$clinvar_id, " ", v$nombre)
  }
  clinvar <- do.call(rbind, filas)
  clinvar$fecha_consulta <- format(Sys.Date())
  utils::write.csv(clinvar, ruta_csv, row.names = FALSE)
}

# ---- Variantes de fondo sinteticas fuera de los genes del panel -------------
largos <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
            159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
            114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
            58617616, 64444167, 46709983, 50818468)
names(largos) <- as.character(1:22)
largo_x <- 156040895
co <- genes_coordenadas[genes_coordenadas$build == "GRCh38", ]

fuera_de_genes <- function(chr, pos) {
  !vapply(seq_along(pos), function(i) {
    any(co$chr == chr[i] & pos[i] >= co$inicio - 1000 & pos[i] <= co$fin + 1000)
  }, logical(1))
}
sortear <- function(n, chr_posibles, pesos, minimo = 1, maximo = NULL) {
  chr <- sample(chr_posibles, n, replace = TRUE, prob = pesos)
  tope <- if (is.null(maximo)) largos[chr] else rep(maximo, n)
  pos <- as.integer(round(stats::runif(n, minimo, tope - 10)))
  ok <- fuera_de_genes(chr, pos)
  data.frame(chr = chr[ok], pos = pos[ok], stringsAsFactors = FALSE)
}
auto <- sortear(4000, names(largos), largos)
x <- sortear(140, "X", 1, minimo = 3e6, maximo = 1.5e8)

transicion <- c(A = "G", G = "A", C = "T", T = "C")
alelos <- function(n) {
  ref <- sample(c("A", "C", "G", "T"), n, replace = TRUE)
  alt <- ifelse(stats::runif(n) < 0.74, transicion[ref],
                vapply(ref, function(r) sample(setdiff(c("A", "C", "G", "T"), c(r, transicion[[r]])), 1), ""))
  data.frame(ref = ref, alt = unname(alt), stringsAsFactors = FALSE)
}
fondo <- rbind(
  cbind(auto, alelos(nrow(auto)), genotipo = ifelse(stats::runif(nrow(auto)) < 0.62, "0/1", "1/1"),
        origen = "sintetico"),
  cbind(x, alelos(nrow(x)), genotipo = "1/1", origen = "sintetico")
)
reales <- data.frame(chr = clinvar$chr, pos = clinvar$pos, ref = clinvar$ref, alt = clinvar$alt,
                     genotipo = clinvar$genotipo,
                     origen = paste0("clinvar;CLINVAR_ID=", clinvar$clinvar_id))
todas <- rbind(fondo, reales)

# Calidad simulada
n <- nrow(todas)
dp <- pmax(6L, as.integer(round(stats::rnorm(n, 45, 14))))
gq <- pmin(99L, pmax(8L, as.integer(round(stats::rnorm(n, 86, 14)))))
k <- ifelse(todas$genotipo == "0/1", as.integer(round(dp * stats::runif(n, 0.35, 0.65))), dp)
ad <- paste0(dp - k, ",", k)
qual <- as.integer(round(stats::runif(n, 40, 900)))
es_real <- startsWith(todas$origen, "clinvar")
dp[es_real] <- 52L; gq[es_real] <- 99L
ad[es_real] <- ifelse(todas$genotipo[es_real] == "0/1", "25,27", "0,52")

orden_chr <- match(todas$chr, c(1:22, "X"))
o <- order(orden_chr, todas$pos)
lineas <- sprintf("chr%s\t%d\t.\t%s\t%s\t%d\tPASS\tORIGEN=%s\tGT:AD:DP:GQ\t%s:%s:%d:%d",
                  todas$chr[o], todas$pos[o], todas$ref[o], todas$alt[o], qual[o],
                  todas$origen[o], todas$genotipo[o], ad[o], dp[o], gq[o])
encabezado <- c(
  "##fileformat=VCFv4.2",
  "##source=ieiprio data-raw/ejemplo_vcf.R",
  "##reference=GRCh38",
  sprintf("##contig=<ID=chr%s,length=%d>", c(names(largos), "X"), c(largos, largo_x)),
  '##INFO=<ID=ORIGEN,Number=1,Type=String,Description="clinvar o sintetico">',
  '##INFO=<ID=CLINVAR_ID,Number=1,Type=Integer,Description="Identificador de variacion en ClinVar">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
  '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
  '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEJEMPLO01"
)
con <- gzfile("inst/extdata/ejemplo.vcf.gz", "w")
writeLines(c(encabezado, lineas), con)
close(con)
message("Listo. inst/extdata/ejemplo.vcf.gz con ", length(lineas), " variantes (",
        sum(es_real), " de ClinVar).")
