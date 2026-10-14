# Dot plot de secuencias por ventanas, en R base.

# Convierte cualquier entrada de secuencia en una cadena de bases.
.secuencia_de <- function(x, nombre_arg = "x") {
  if (inherits(x, c("iei_ab1", "iei_scf"))) return(list(bases = toupper(x$bases), nombre = x$archivo %||% "lectura"))
  if (!is.character(x) || length(x) != 1) {
    cli::cli_abort("{.arg {nombre_arg}} debe ser una secuencia, una ruta a FASTA, .ab1 o .scf, o un resultado de {.fn iei_sanger}.")
  }
  nombre <- "secuencia"
  if (file.exists(x)) {
    nombre <- basename(x)
    ext <- tolower(tools::file_ext(x))
    if (ext %in% c("ab1", "abi")) return(list(bases = toupper(iei_leer_ab1(x)$bases), nombre = basename(x)))
    if (ext == "scf") return(list(bases = toupper(iei_leer_scf(x)$bases), nombre = basename(x)))
    lineas <- readLines(x, warn = FALSE)
    cab <- which(startsWith(lineas, ">"))
    if (length(cab) > 1) lineas <- lineas[seq_len(cab[2] - 1)]
    x <- paste(lineas[!startsWith(lineas, ">")], collapse = "")
  }
  bases <- toupper(gsub("[^A-Za-z]", "", x))
  if (nchar(bases) < 10 || grepl("[^ACGTNRYSWKMBDHV]", bases)) {
    cli::cli_abort("{.arg {nombre_arg}} no es un archivo ni una secuencia de ADN v\u00e1lida.")
  }
  list(bases = bases, nombre = nombre)
}

# Puntos del dot plot: posiciones (i, j) donde la ventana de largo `ventana`
# que empieza en i de `a` y en j de `b` coincide en al menos `umbral` bases.
.puntos_dotplot <- function(a, b, ventana, umbral) {
  av <- strsplit(a, "")[[1]]; bv <- strsplit(b, "")[[1]]
  n <- length(av) - ventana + 1L; m <- length(bv) - ventana + 1L
  if (n < 1 || m < 1) cli::cli_abort("Las secuencias son m\u00e1s cortas que la {.arg ventana}.")
  iguales <- outer(av, bv, "==") & outer(av != "N", bv != "N")
  cuenta <- matrix(0L, n, m)
  for (k in seq_len(ventana) - 1L) cuenta <- cuenta + iguales[k + seq_len(n), k + seq_len(m), drop = FALSE]
  idx <- which(cuenta >= umbral, arr.ind = TRUE)
  data.frame(i = as.integer(idx[, 1]), j = as.integer(idx[, 2]))
}

#' Dot plot de dos secuencias
#'
#' Compara dos secuencias con el metodo clasico de ventanas. Se marca un punto
#' en (i, j) cuando la ventana de `ventana` bases que empieza en la posicion i
#' de `x` y la que empieza en la posicion j de `y` coinciden en al menos
#' `umbral` bases. Dos secuencias iguales forman una diagonal. Una sustitucion
#' apenas la interrumpe, una insercion o delecion la corta y la continua
#' desplazada, y una repeticion aparece como diagonales paralelas.
#'
#' Si `y` no se da, `x` se compara con su referencia: la del alineamiento si
#' `x` viene de [iei_sanger()], o la de `gen`, `region` o `referencia`. Sin
#' nada de eso, `x` se compara consigo misma, lo que muestra sus repeticiones.
#' La lectura se orienta sola si esta en la hebra reversa.
#'
#' @param x Secuencia escrita, ruta a un FASTA, `.ab1` o `.scf`, o un
#'   resultado de [iei_sanger()].
#' @param y Opcional, otra secuencia en cualquiera de esas formas.
#' @param gen,region,referencia,build,cache Para comparar con la referencia
#'   genomica, como en [iei_sanger()].
#' @param ventana Largo de la ventana en bases.
#' @param umbral Coincidencias minimas dentro de la ventana para marcar un
#'   punto.
#' @param variante Fila de `x` cuya posicion se marca, si `x` viene de
#'   [iei_sanger()]. `NULL` no marca nada.
#' @param lectura Lectura de `x` a usar si trae varias (nombre o numero).
#' @return Un grafico de ggplot2. Los puntos quedan en el atributo `puntos`.
#' @export
#' @examples
#' ab1 <- iei_ejemplo_sanger()
#' ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
#' v <- iei_sanger(ab1, gen = "BTK", sexo = "M", referencia = ref)
#' iei_dotplot(v)
#' iei_dotplot("ACGTTGCAACGTAGGCTTACGATCGATCGGATCCATGCA",
#'             "ACGTTGCAACGTAGGCTTACGGGGATCGATCGGATCCATGCA", ventana = 6, umbral = 5)
iei_dotplot <- function(x, y = NULL, gen = NULL, region = NULL, referencia = NULL,
                        ventana = 10, umbral = 8, variante = 1, lectura = 1,
                        build = c("GRCh38", "GRCh37"), cache = tools::R_user_dir("ieiprio", "cache")) {
  .necesita_ggplot()
  build <- match.arg(build)
  ventana <- as.integer(ventana); umbral <- as.integer(umbral)
  if (ventana < 1 || umbral < 1 || umbral > ventana) {
    cli::cli_abort("{.arg umbral} debe estar entre 1 y {.arg ventana}.")
  }
  marca <- NULL
  eje_y <- "Posici\u00f3n en y"
  desde_y <- 1L
  info <- if (is.data.frame(x)) attr(x, "ieiprio_sanger") else NULL
  if (is.data.frame(x) && is.null(info)) cli::cli_abort("{.arg x} no viene de {.fn iei_sanger}.")

  if (!is.null(info)) {
    l <- info[[lectura]]
    if (is.null(l)) cli::cli_abort("No encuentro la lectura {.val {lectura}}.")
    al <- l$alineamiento
    a <- paste(al$q[!is.na(al$iq)], collapse = "")
    nombre_x <- sprintf("%s (%s)", l$lectura, l$hebra)
    if (is.null(y)) {
      con_ref <- !is.na(al$jr)
      b <- paste(al$r[con_ref], collapse = "")
      desde_y <- min(al$pos, na.rm = TRUE)
      gen_x <- if ("gen" %in% names(x)) stats::na.omit(x$gen)[1] else NA
      nombre_y <- if (is.na(gen_x)) "la referencia" else sprintf("la referencia de %s", gen_x)
      eje_y <- sprintf("Cromosoma %s", sub(":.*$", "", l$region))
      if (!is.null(variante) && nrow(x) >= variante) {
        p <- x$pos[variante]
        fila <- which(al$pos == p & !is.na(al$iq))
        if (length(fila) > 0) {
          marca <- data.frame(i = sum(!is.na(al$iq[seq_len(fila[1])])), j = p,
                              etiqueta = sprintf("%s>%s", x$ref[variante], x$alt[variante]))
        }
      }
    }
  } else {
    sx <- .secuencia_de(x, "x")
    a <- sx$bases
    nombre_x <- sx$nombre
  }

  if (!is.null(y)) {
    sy <- .secuencia_de(y, "y")
    b <- sy$bases
    nombre_y <- sy$nombre
  } else if (is.null(info) && (!is.null(gen) || !is.null(region) || !is.null(referencia))) {
    if (!is.null(gen)) {
      g <- .normalizar_simbolo(gen, ieiprio::iuis_panel$gen, ieiprio::iuis_alias)
      if (!is.na(g)) gen <- g
    }
    ref <- if (is.null(referencia)) .referencia(gen, region, build, cache) else
      .referencia_propia(referencia, gen, region, build)
    ubic <- .ubicar(a, ref$secuencia)
    if (is.null(ubic)) cli::cli_abort("No pude ubicar {.arg x} en la referencia de {ref$nombre}.")
    if (ubic$hebra == "-") a <- .revcomp(a)
    b <- substr(ref$secuencia, ubic$inicio, ubic$fin)
    desde_y <- ref$inicio + ubic$inicio - 1L
    nombre_y <- sprintf("la referencia de %s", ref$nombre)
    eje_y <- sprintf("Cromosoma %s", ref$chr)
  } else if (is.null(info)) {
    b <- a
    nombre_y <- nombre_x
  }

  puntos <- .puntos_dotplot(a, b, ventana, umbral)
  invertida <- FALSE
  # Si y esta en la otra hebra, la comparacion directa casi no da puntos.
  if (!is.null(y)) {
    inversa <- .puntos_dotplot(a, .revcomp(b), ventana, umbral)
    if (nrow(inversa) > 2 * nrow(puntos)) {
      puntos <- inversa
      invertida <- TRUE
    }
  }
  puntos$j <- puntos$j + desde_y - 1L

  subtitulo <- sprintf("Ventana de %d bases con al menos %d coincidencias%s.", ventana, umbral,
                       if (invertida) ", y en complemento reverso" else "")
  g <- ggplot2::ggplot(puntos, ggplot2::aes(x = .data$i, y = .data$j)) +
    ggplot2::geom_point(size = 0.35, colour = .marca, shape = 15) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE),
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::labs(title = sprintf("%s frente a %s", nombre_x, nombre_y), subtitle = subtitulo,
                  x = sprintf("Posici\u00f3n en %s", nombre_x), y = eje_y) +
    .tema_iei(base_size = 11) +
    ggplot2::theme(aspect.ratio = 3 / 4, panel.border = ggplot2::element_rect(fill = NA, colour = .tinta[["eje"]]))
  if (!is.null(marca)) {
    g <- g +
      ggplot2::geom_vline(xintercept = marca$i, colour = "#e34948", linewidth = 0.4, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = marca$j, colour = "#e34948", linewidth = 0.4, linetype = "dashed") +
      ggplot2::geom_point(data = marca, colour = "#e34948", size = 3, shape = 21, stroke = 1.1) +
      ggplot2::geom_label(data = marca, ggplot2::aes(label = .data$etiqueta), colour = "#a32626",
                          hjust = -0.15, vjust = 1.4, size = 3.3, fontface = "bold", fill = "white",
                          linewidth = 0)
  }
  attr(g, "puntos") <- puntos
  g
}

#' Graficos de un resultado de iei_sanger
#'
#' `plot()` sobre la salida de [iei_sanger()] dibuja el cromatograma de la
#' variante ([iei_cromatograma()]) o el dot plot de la lectura frente a la
#' referencia ([iei_dotplot()]).
#'
#' @param x Resultado de [iei_sanger()].
#' @param tipo `"cromatograma"` o `"dotplot"`.
#' @param ... Otros argumentos para [iei_cromatograma()] o [iei_dotplot()].
#' @return Un grafico de ggplot2.
#' @export
#' @examples
#' ab1 <- iei_ejemplo_sanger()
#' ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
#' v <- iei_sanger(ab1, gen = "BTK", sexo = "M", referencia = ref)
#' plot(v, tipo = "dotplot")
plot.iei_sanger <- function(x, tipo = c("cromatograma", "dotplot"), ...) {
  tipo <- match.arg(tipo)
  g <- switch(tipo, cromatograma = iei_cromatograma(x, ...), dotplot = iei_dotplot(x, ...))
  g
}
