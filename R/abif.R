# Lector del formato ABIF (.ab1 / .abi) de los secuenciadores Sanger de
# Applied Biosystems, escrito en R base. Especificacion: "Applied Biosystems
# Genetic Analysis Data File Format" (ABIF). Todos los enteros son big-endian.

.abif_leer_entradas <- function(bytes) {
  if (length(bytes) < 34 || rawToChar(bytes[1:4]) != "ABIF") {
    cli::cli_abort("El archivo no es un cromatograma ABIF (.ab1).")
  }
  int32 <- function(pos) readBin(bytes[pos:(pos + 3)], "integer", size = 4, endian = "big")
  int16 <- function(pos) readBin(bytes[pos:(pos + 1)], "integer", size = 2, endian = "big")
  # La entrada raiz del directorio esta en el byte 7 (base 1).
  n_entradas <- int32(7 + 12)
  inicio_dir <- int32(7 + 20)
  entradas <- lapply(seq_len(n_entradas), function(k) {
    p <- inicio_dir + (k - 1) * 28 + 1
    list(
      nombre = rawToChar(bytes[p:(p + 3)]),
      numero = int32(p + 4),
      tipo = int16(p + 8),
      tam_elemento = int16(p + 10),
      n_elementos = int32(p + 12),
      tam_datos = int32(p + 16),
      # Datos de 4 bytes o menos van dentro del propio campo de posicion.
      inicio_datos = if (int32(p + 16) <= 4) p + 20 else int32(p + 20) + 1
    )
  })
  stats::setNames(entradas, vapply(entradas, function(e) paste0(e$nombre, e$numero), ""))
}

.abif_bytes <- function(bytes, e) {
  if (is.null(e) || e$tam_datos == 0) return(NULL)
  bytes[e$inicio_datos:(e$inicio_datos + e$tam_datos - 1)]
}

.abif_valor <- function(bytes, e) {
  datos <- .abif_bytes(bytes, e)
  if (is.null(datos)) return(NULL)
  switch(as.character(e$tipo),
    "1" = as.integer(datos),                                             # byte sin signo
    "2" = rawToChar(datos[datos != as.raw(0)]),                          # char
    "3" = readBin(datos, "integer", n = e$n_elementos, size = 2, signed = FALSE, endian = "big"),
    "4" = readBin(datos, "integer", n = e$n_elementos, size = 2, endian = "big"),
    "5" = readBin(datos, "integer", n = e$n_elementos, size = 4, endian = "big"),
    "7" = readBin(datos, "double", n = e$n_elementos, size = 4, endian = "big"),
    "8" = readBin(datos, "double", n = e$n_elementos, size = 8, endian = "big"),
    "18" = rawToChar(datos[-1][datos[-1] != as.raw(0)]),                 # pString
    "19" = rawToChar(datos[datos != as.raw(0)]),                         # cString
    NULL
  )
}

#' Leer un cromatograma Sanger (.ab1)
#'
#' Lee un archivo ABIF de un secuenciador Applied Biosystems y devuelve las
#' bases llamadas, su calidad, la posicion de cada pico y las cuatro trazas.
#' Esta escrito en R base, sin dependencias de Bioconductor.
#'
#' @param archivo Ruta al archivo `.ab1` o `.abi`.
#' @return Una lista de clase `iei_ab1` con `bases` (texto), `calidad`
#'   (Phred por base, `NA` si el archivo no la trae), `picos` (posicion de
#'   cada base en la traza), `trazas` (tibble con las columnas `A`, `C`, `G`,
#'   `T`), `muestra` y `archivo`.
#' @export
#' @examples
#' \dontrun{
#' ab1 <- iei_leer_ab1("BTK_exon15_F.ab1")
#' ab1$bases
#' }
iei_leer_ab1 <- function(archivo) {
  if (!file.exists(archivo)) cli::cli_abort("No encuentro {.path {archivo}}.")
  bytes <- readBin(archivo, "raw", n = file.size(archivo))
  e <- .abif_leer_entradas(bytes)
  val <- function(clave) .abif_valor(bytes, e[[clave]])

  bases <- val("PBAS2") %||% val("PBAS1")
  picos <- val("PLOC2") %||% val("PLOC1")
  if (is.null(bases) || is.null(picos)) {
    cli::cli_abort("{.path {basename(archivo)}} no trae bases llamadas (PBAS) ni posiciones de picos (PLOC).")
  }
  orden <- strsplit(val("FWO_1") %||% "GATC", "")[[1]]
  canales <- lapply(9:12, function(k) val(paste0("DATA", k)))
  if (any(vapply(canales, is.null, logical(1)))) {
    cli::cli_abort("{.path {basename(archivo)}} no trae las trazas analizadas (DATA9 a DATA12).")
  }
  trazas <- stats::setNames(canales, orden)[c("A", "C", "G", "T")]
  # La calidad se guarda como "char", pero cada byte es un valor Phred.
  calidad <- as.integer(.abif_bytes(bytes, e[["PCON2"]]) %||% .abif_bytes(bytes, e[["PCON1"]]))
  if (length(calidad) != nchar(bases)) calidad <- rep(NA_integer_, nchar(bases))

  structure(list(
    bases = bases, calidad = as.integer(calidad), picos = as.integer(picos) + 1L,
    trazas = tibble::as_tibble(trazas),
    muestra = val("SMPL1") %||% NA_character_,
    archivo = basename(archivo)
  ), class = "iei_ab1")
}

#' @rdname iei_leer_ab1
#' @param x Un objeto `iei_ab1`.
#' @param ... No se usa.
#' @export
print.iei_ab1 <- function(x, ...) {
  q <- x$calidad
  cli::cli_text("{.strong Cromatograma Sanger} {.file {x$archivo}}")
  cli::cli_text("{nchar(x$bases)} bases {cli::symbol$bullet} ",
                if (all(is.na(q))) "sin calidad" else
                  "calidad mediana {stats::median(q, na.rm = TRUE)}, {round(100 * mean(q >= 20, na.rm = TRUE))}% con Phred 20 o mas",
                " {cli::symbol$bullet} {nrow(x$trazas)} puntos de traza")
  cli::cli_text("{.emph {substr(x$bases, 1, 60)}...}")
  invisible(x)
}
