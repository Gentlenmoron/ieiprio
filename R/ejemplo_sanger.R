# Escritor minimo del formato ABIF. Lo usan iei_ejemplo_sanger() y las
# pruebas. Cada base es un pico gaussiano; en las posiciones de `het` se
# agrega un segundo pico de la base indicada con 80% de la altura.
.escribir_ab1 <- function(ruta, secuencia, het = list(), calidad = NULL,
                          separacion = 12L, muestra = "prueba", alto = NULL) {
  bases <- strsplit(secuencia, "")[[1]]
  n <- length(bases)
  if (is.null(calidad)) calidad <- rep(40L, n)
  picos <- as.integer(separacion * seq_len(n))          # base 0 en el archivo
  largo <- max(picos) + separacion
  trazas <- matrix(0, largo, 4, dimnames = list(NULL, c("G", "A", "T", "C")))
  x <- seq_len(largo) - 1
  pico <- function(centro, alto) alto * exp(-((x - centro)^2) / (2 * 2^2))
  if (is.null(alto)) alto <- rep(1000, n)
  for (i in seq_len(n)) {
    trazas[, bases[i]] <- trazas[, bases[i]] + pico(picos[i], alto[i])
    segunda <- het[[as.character(i)]]
    if (!is.null(segunda)) trazas[, segunda] <- trazas[, segunda] + pico(picos[i], 0.8 * alto[i])
  }
  trazas <- round(trazas)

  i16 <- function(v) writeBin(as.integer(v), raw(), size = 2, endian = "big")
  i32 <- function(v) writeBin(as.integer(v), raw(), size = 4, endian = "big")
  entradas <- list(
    list("PBAS", 2L, 2L, 1L, charToRaw(secuencia)),
    list("PLOC", 2L, 4L, 2L, i16(picos)),
    list("PCON", 2L, 2L, 1L, as.raw(calidad)),
    list("DATA", 9L, 4L, 2L, i16(trazas[, "G"])),
    list("DATA", 10L, 4L, 2L, i16(trazas[, "A"])),
    list("DATA", 11L, 4L, 2L, i16(trazas[, "T"])),
    list("DATA", 12L, 4L, 2L, i16(trazas[, "C"])),
    list("FWO_", 1L, 2L, 1L, charToRaw("GATC")),
    list("SMPL", 1L, 18L, 1L, c(as.raw(nchar(muestra)), charToRaw(muestra)))
  )
  datos <- raw(0)
  directorio <- raw(0)
  desplazamiento <- 128L
  for (e in entradas) {
    d <- e[[5]]
    n_el <- length(d) %/% e[[4]]
    if (length(d) <= 4) {
      campo <- c(d, raw(4 - length(d)))
    } else {
      campo <- i32(desplazamiento + length(datos))
      datos <- c(datos, d)
    }
    directorio <- c(directorio, charToRaw(e[[1]]), i32(e[[2]]), i16(e[[3]]), i16(e[[4]]),
                    i32(n_el), i32(length(d)), campo, i32(0))
  }
  inicio_dir <- 128L + length(datos)
  cabecera <- c(charToRaw("ABIF"), i16(101), charToRaw("tdir"), i32(1), i16(1023), i16(28),
                i32(length(entradas)), i32(length(directorio)), i32(inicio_dir), i32(0))
  cabecera <- c(cabecera, raw(128 - length(cabecera)))
  writeBin(c(cabecera, datos, directorio), ruta)
  invisible(ruta)
}


#' Cromatogramas Sanger de ejemplo
#'
#' Genera dos cromatogramas `.ab1` (directo y reverso) de un amplicon de 500
#' pb de BTK que contiene la variante patogenica real NM_000061.3:c.1632-2A>G
#' (ClinVar 265455), la misma del VCF de ejemplo. Sirven para probar
#' [iei_sanger()] y [iei_cromatograma()] sin tener archivos propios. La
#' secuencia sale de la referencia GRCh38 de BTK que trae el paquete
#' (`inst/extdata/btk_referencia.fa`, bajada de Ensembl), asi que no necesita
#' conexion.
#'
#' Los picos son simulados: el cromatograma es limpio salvo en los extremos,
#' donde la calidad cae como en una lectura real.
#'
#' @param genotipo `"hemi"` para un varon afectado o `"het"` para una
#'   portadora, que se ve como pico doble.
#' @param dir Carpeta donde escribir los archivos.
#' @return Las rutas de los dos archivos, de forma invisible.
#' @export
#' @examples
#' ab1 <- iei_ejemplo_sanger()
#' ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
#' v <- iei_sanger(ab1, gen = "BTK", sexo = "M", referencia = ref)
#' v
iei_ejemplo_sanger <- function(genotipo = c("hemi", "het"), dir = tempfile("sanger_")) {
  genotipo <- match.arg(genotipo)
  pos <- 101353990L                      # X:101353990 T>C en GRCh38
  ref <- .referencia_propia(system.file("extdata", "btk_referencia.fa", package = "ieiprio"), gen = "BTK")
  amplicon <- strsplit(substr(ref$secuencia, pos - 250L - ref$inicio + 1L, pos + 249L - ref$inicio + 1L), "")[[1]]
  if (amplicon[251] != "T") cli::cli_abort("La referencia de BTK no coincide con la esperada en X:{pos}.")
  het <- list()
  if (genotipo == "hemi") amplicon[251] <- "C" else het[["251"]] <- "C"
  n <- length(amplicon)
  # Calidad baja en los extremos y algo de variacion en la altura de los picos.
  calidad <- as.integer(pmin(45, pmax(4, round(45 * pmin(seq_len(n) / 40, rev(seq_len(n)) / 50)))))
  alto <- round(1000 + 300 * sin(seq_len(n) * 1.7))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  directa <- file.path(dir, "BTK_c1632-2_F.ab1")
  reversa <- file.path(dir, "BTK_c1632-2_R.ab1")
  .escribir_ab1(directa, paste(amplicon, collapse = ""), het = het, calidad = calidad,
                muestra = "BTK_F", alto = alto)
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  het_r <- list()
  if (genotipo == "het") het_r[[as.character(n - 251L + 1L)]] <- "G"
  .escribir_ab1(reversa, paste(rev(comp[amplicon]), collapse = ""), het = het_r,
                calidad = rev(calidad), muestra = "BTK_R", alto = rev(alto))
  cli::cli_inform("Cromatogramas de ejemplo en {.path {dir}}.")
  invisible(c(directa, reversa))
}
