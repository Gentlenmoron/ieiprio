# Sanger avanzado: lectura de archivos SCF, recorte de cebadores, llamado de
# indels heterocigotas por deconvolucion de picos dobles y secuencia
# consenso.

# ---- SCF -------------------------------------------------------------------

#' Leer un cromatograma SCF
#'
#' Lee un archivo SCF (versiones 2 y 3), el formato abierto que exportan
#' muchos secuenciadores antiguos y programas como Staden, Chromas o 4Peaks.
#' Devuelve el mismo objeto que [iei_leer_ab1()], asi que sirve para
#' [iei_sanger()] y para graficar con `plot()`.
#'
#' @param archivo Ruta al archivo `.scf`.
#' @return Un objeto `iei_ab1` con `bases`, `calidad`, `picos`, `trazas`,
#'   `muestra` y `archivo`.
#' @export
#' @examples
#' \dontrun{
#' scf <- iei_leer_scf("BTK_exon15_F.scf")
#' plot(scf)
#' }
iei_leer_scf <- function(archivo) {
  if (!file.exists(archivo)) cli::cli_abort("No encuentro {.path {archivo}}.")
  b <- readBin(archivo, "raw", n = file.size(archivo))
  if (length(b) < 128 || !identical(rawToChar(b[1:4]), ".scf")) {
    cli::cli_abort("{.path {basename(archivo)}} no es un archivo SCF.")
  }
  u32 <- function(off) readBin(b[(off + 1):(off + 4)], "integer", size = 4, endian = "big")
  n_puntos <- u32(4); off_puntos <- u32(8); n_bases <- u32(12); off_bases <- u32(24)
  version <- suppressWarnings(as.numeric(rawToChar(b[37:40])))
  tam <- u32(40)
  if (!tam %in% c(1L, 2L)) tam <- 2L
  leer_u <- function(off, n, size) {
    if (n == 0) return(integer(0))
    v <- b[(off + 1):(off + n * size)]
    if (size == 1) as.integer(v) else readBin(v, "integer", size = 2, n = n, signed = FALSE, endian = "big")
  }
  modulo <- if (tam == 1) 256 else 65536
  if (!is.na(version) && version >= 3) {
    canales <- lapply(0:3, function(k) {
      v <- leer_u(off_puntos + k * n_puntos * tam, n_puntos, tam)
      # Doble codificacion delta: se deshace con dos sumas acumuladas.
      v <- cumsum(v) %% modulo
      cumsum(v) %% modulo
    })
    picos <- readBin(b[(off_bases + 1):(off_bases + 4 * n_bases)], "integer", size = 4, n = n_bases, endian = "big")
    p <- off_bases + 4 * n_bases
    prob <- lapply(0:3, function(k) as.integer(b[(p + k * n_bases + 1):(p + (k + 1) * n_bases)]))
    bases <- rawToChar(b[(p + 4 * n_bases + 1):(p + 5 * n_bases)])
  } else {
    todo <- leer_u(off_puntos, 4 * n_puntos, tam)
    canales <- lapply(1:4, function(k) todo[seq(k, length(todo), by = 4)])
    reg <- matrix(b[(off_bases + 1):(off_bases + 12 * n_bases)], nrow = 12)
    picos <- vapply(seq_len(n_bases), function(i) readBin(reg[1:4, i], "integer", size = 4, endian = "big"), integer(1))
    prob <- lapply(5:8, function(k) as.integer(reg[k, ]))
    bases <- rawToChar(reg[9, ])
  }
  trazas <- tibble::tibble(A = canales[[1]], C = canales[[2]], G = canales[[3]], T = canales[[4]])
  vb <- strsplit(toupper(bases), "")[[1]]
  canal <- match(vb, c("A", "C", "G", "T"))
  calidad <- vapply(seq_along(vb), function(i) if (is.na(canal[i])) 0L else prob[[canal[i]]][i], integer(1))
  structure(list(bases = toupper(bases), calidad = calidad, picos = as.integer(picos) + 1L,
                 trazas = trazas, muestra = NA_character_, archivo = basename(archivo)),
            class = "iei_ab1")
}

# Escribe un SCF version 3 (con muestras de 2 bytes). Se usa en las pruebas y
# para exportar cromatogramas sinteticos.
.scf_escribir <- function(a, archivo) {
  n <- nrow(a$trazas); nb <- nchar(a$bases)
  u32 <- function(v) writeBin(as.integer(v), raw(), size = 4, endian = "big")
  canales <- lapply(c("A", "C", "G", "T"), function(k) {
    v <- as.integer(round(a$trazas[[k]])) %% 65536
    v <- c(v[1], diff(v)) %% 65536
    v <- c(v[1], diff(v)) %% 65536
    writeBin(as.integer(v), raw(), size = 2, endian = "big")
  })
  muestras <- do.call(c, canales)
  vb <- strsplit(a$bases, "")[[1]]
  cal <- ifelse(is.na(a$calidad), 0L, as.integer(a$calidad))
  prob <- lapply(c("A", "C", "G", "T"), function(k) as.raw(ifelse(vb == k, pmin(cal, 255L), 0L)))
  bases <- c(u32(a$picos - 1L), do.call(c, prob), charToRaw(a$bases), raw(3 * nb))
  off_m <- 128L; off_b <- off_m + length(muestras); off_c <- off_b + length(bases)
  cab <- c(charToRaw(".scf"), u32(n), u32(off_m), u32(nb), u32(0), u32(0), u32(off_b),
           u32(0), u32(off_c), charToRaw("3.00"), u32(2), u32(0), u32(0), u32(off_c), raw(72))
  writeBin(c(cab, muestras, bases), archivo)
  invisible(archivo)
}

# ---- Cebadores -------------------------------------------------------------

# Indices de la lectura que quedan despues de quitar los cebadores: un
# cebador cerca del inicio se quita junto con lo anterior, y el complemento
# reverso de un cebador cerca del final se quita junto con lo posterior. Las
# bases del cebador copian al oligonucleotido, no al paciente, asi que una
# variante ahi no es real.
.recortar_cebadores <- function(bases, cebadores, max_dif = 2L, zona = 40L) {
  n <- nchar(bases)
  ini <- 1L; fin <- n
  if (length(cebadores) == 0) return(seq_len(n))
  for (p in toupper(gsub("[^A-Za-z]", "", cebadores))) {
    if (nchar(p) < 10) next
    for (patron in c(p, .revcomp(p))) {
      m <- utils::aregexec(patron, bases, fixed = TRUE,
                           max.distance = list(all = max_dif, insertions = 0, deletions = 0,
                                               substitutions = max_dif))[[1]]
      if (m[1] < 1) next
      desde <- m[1]; hasta <- m[1] + attr(m, "match.length")[1] - 1L
      if (desde <= zona) ini <- max(ini, hasta + 1L)
      else if (hasta >= n - zona) fin <- min(fin, desde - 1L)
    }
  }
  if (fin < ini) integer(0) else ini:fin
}

# ---- Indels heterocigotas --------------------------------------------------

# Despues de una insercion o delecion heterocigota la lectura es la suma de
# dos secuencias desfasadas. Se prueba cada desfase posible (delecion o
# insercion de 1 a `max_largo` bases) y se elige el que explica los picos
# dobles: en cada posicion los dos picos deben ser la base de referencia y
# la base de referencia corrida. Es la misma idea de Indigo, TIDE o Poly Peak
# Parser.
.deconvolucionar_indel <- function(lect, al, hebra, ref, inicio_idx, max_largo = 30L,
                                   ventana = 60L, min_puntaje = 0.8) {
  if (is.null(lect$trazas) || is.null(lect$picos) || is.na(inicio_idx)) return(NULL)
  s <- .senal_picos(lect)
  ok <- !is.na(al$pos) & !is.na(al$iq_original) & al$iq_original < inicio_idx & al$q == al$r
  if (!any(ok)) return(NULL)
  k <- which(ok)[which.max(al$iq_original[ok])]
  i0 <- al$iq_original[k]; p0 <- al$pos[k]
  base_gen <- function(p) {
    j <- p - ref$inicio + 1L
    if (j < 1 || j > nchar(ref$secuencia)) "" else substr(ref$secuencia, j, j)
  }
  comp <- function(x) if (x == "") "" else .complemento[[x]]
  M <- ventana + max_largo + 1L
  R <- if (hebra == "+") vapply(seq_len(M), function(m) base_gen(p0 + m), "") else
    vapply(seq_len(M), function(m) comp(base_gen(p0 - m)), "")
  idx <- i0 + seq_len(min(ventana, nrow(s) - i0))
  if (length(idx) < 20) return(NULL)
  conj <- lapply(idx, function(i) {
    v <- s[i, ]
    o <- order(v, decreasing = TRUE)
    if (v[o[2]] >= 0.2 * v[o[1]]) colnames(s)[o[1:2]] else colnames(s)[o[1]]
  })
  L <- length(idx)
  explica <- function(m, a, b) {
    cc <- conj[[m]]
    if (a == "" || b == "") return(NA)
    if (length(cc) == 2) a != b && setequal(cc, c(a, b)) else a == b && a == cc
  }
  puntaje <- function(d, tipo) {
    m <- seq_len(L)
    if (tipo == "ins") m <- m[m > d]
    b <- if (tipo == "del") R[m + d] else R[m - d]
    mean(mapply(explica, m, R[m], b), na.rm = TRUE)
  }
  cand <- data.frame(tipo = rep(c("del", "ins"), each = max_largo), d = rep(seq_len(max_largo), 2),
                     stringsAsFactors = FALSE)
  cand$puntaje <- mapply(puntaje, cand$d, cand$tipo)
  cand <- cand[order(-cand$puntaje, cand$d), ]
  mejor <- cand[1, ]
  if (is.na(mejor$puntaje) || mejor$puntaje < min_puntaje) return(NULL)
  d <- mejor$d
  secuencia <- function(desde, hasta) paste(vapply(desde:hasta, base_gen, ""), collapse = "")
  if (mejor$tipo == "del") {
    pos <- if (hebra == "+") p0 else p0 - d - 1L
    r <- secuencia(pos, pos + d)
    a <- base_gen(pos)
  } else {
    ins <- paste(vapply(seq_len(d), function(m) {
      o <- setdiff(conj[[m]], R[m])
      if (length(o) > 0) o[1] else R[m]
    }, ""), collapse = "")
    pos <- if (hebra == "+") p0 else p0 - 1L
    r <- base_gen(pos)
    a <- paste0(r, if (hebra == "+") ins else .revcomp(ins))
  }
  norm <- .normalizar_indels(data.frame(pos = pos, ref = r, alt = a, stringsAsFactors = FALSE), ref)
  list(pos = as.integer(norm$pos), ref = norm$ref, alt = norm$alt, puntaje = round(mejor$puntaje, 3),
       iq = i0, tipo = mejor$tipo, largo = d)
}

# ---- Consenso --------------------------------------------------------------

#' Secuencia consenso de las lecturas Sanger
#'
#' Junta las lecturas alineadas por [iei_sanger()] en una sola secuencia en
#' coordenadas del genoma. En cada posicion gana la base con mas calidad
#' sumada entre las lecturas; si dos lecturas buenas no coinciden, la
#' posicion se marca como discordante y se escribe con el codigo IUPAC de las
#' dos bases. Las inserciones respecto de la referencia no entran en el
#' consenso; estan en la tabla de variantes.
#'
#' @param x Resultado de [iei_sanger()] (o de [iei_anotar()] sobre el).
#' @param archivo Opcional, ruta de un FASTA donde guardar el consenso.
#' @param min_calidad Calidad Phred minima para que una lectura cuente como
#'   discordante.
#' @return Un tibble con `pos`, `ref`, `consenso`, `lecturas` (cuantas cubren
#'   la posicion) y `discordante`. El atributo `secuencia` tiene el consenso
#'   como texto.
#' @export
#' @examples
#' \dontrun{
#' v <- iei_sanger(c("F.ab1", "R.ab1"), gen = "BTK")
#' cons <- iei_consenso(v, "BTK_consenso.fasta")
#' }
iei_consenso <- function(x, archivo = NULL, min_calidad = 20) {
  info <- attr(x, "ieiprio_sanger")
  if (is.null(info)) cli::cli_abort("{.arg x} no viene de {.fn iei_sanger}.")
  partes <- lapply(info, function(l) {
    al <- l$alineamiento
    cal <- if ("calidad" %in% names(al)) al$calidad else rep(NA_real_, nrow(al))
    ok <- !is.na(al$pos)
    data.frame(pos = al$pos[ok], ref = al$r[ok], base = al$q[ok], calidad = cal[ok],
               lectura = l$lectura, stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, partes)
  d$peso <- ifelse(is.na(d$calidad), 20, d$calidad)
  filas <- lapply(split(d, d$pos), function(g) {
    suma <- tapply(g$peso, g$base, sum)
    suma <- sort(suma, decreasing = TRUE)
    buenas <- unique(g$base[g$peso >= min_calidad & g$base != "N"])
    discord <- length(buenas) > 1
    b <- names(suma)[1]
    if (discord && all(buenas %in% c("A", "C", "G", "T")) && length(buenas) == 2) b <- .iupac_de(buenas[1], buenas[2])
    data.frame(pos = g$pos[1], ref = g$ref[1], consenso = b, lecturas = length(unique(g$lectura)),
               discordante = discord, stringsAsFactors = FALSE)
  })
  out <- tibble::as_tibble(do.call(rbind, filas))
  out <- out[order(out$pos), ]
  sec <- paste(out$consenso, collapse = "")
  attr(out, "secuencia") <- sec
  if (!is.null(archivo)) {
    chr <- x$chr[1] %||% ""
    cab <- sprintf(">%s %s:%d-%d consenso de %d lectura%s", x$muestra[1] %||% "sanger", chr,
                   min(out$pos), max(out$pos), length(info), if (length(info) == 1) "" else "s")
    cortes <- seq(1, nchar(sec), by = 60)
    writeLines(c(cab, substring(sec, cortes, pmin(cortes + 59, nchar(sec)))), archivo)
    cli::cli_inform("Consenso de {nrow(out)} bases guardado en {.path {archivo}}.")
  }
  if (any(out$discordante)) {
    cli::cli_warn("{sum(out$discordante)} posici{?\u00f3n/ones} con lecturas discordantes: {.val {utils::head(out$pos[out$discordante], 5)}}.")
  }
  out
}
