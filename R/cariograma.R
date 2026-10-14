# Cariograma con ggplot2, al estilo de karyoploteR.

#' Cariograma del panel
#'
#' Dibuja los 24 cromosomas uno al lado del otro con los genes del panel
#' marcados a la izquierda en el color de su categoria IUIS. Con un resultado
#' de [iei_priorizar()] marca las variantes a la derecha, del color de su
#' prioridad, y con una busqueda de [iei_buscar_fenotipo()] marca los genes
#' candidatos, del color de su herencia. Los mejores se rotulan.
#'
#' @inheritParams iei_circos
#' @return Un objeto ggplot.
#' @seealso [iei_circos()]
#' @export
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   iei_cariograma()
#'   iei_cariograma(iei_buscar_fenotipo("candidiasis mucocutanea cronica", n = 10))
#' }
iei_cariograma <- function(x = NULL, build = c("GRCh38", "GRCh37"), etiquetas = 12, titulo = NULL) {
  .necesita_ggplot()
  build <- match.arg(build)
  if (inherits(x, "iei_resultado")) build <- attr(x, "build") %||% build
  gen <- .genoma(build)
  gen$x <- as.integer(gen$chr)
  panel <- .genes_panel_posicion(build)
  panel$x <- as.integer(.cromosomes_factor(panel$chr))
  panel$categoria_nombre <- factor(.nombres_iuis[panel$categoria], levels = .nombres_iuis)
  marcas <- .marcas_genoma(x, build)
  mb <- 1e6
  brazo <- rbind(
    data.frame(x = gen$x, y0 = 0, y1 = gen$centromero - 1.6e6),
    data.frame(x = gen$x, y0 = gen$centromero + 1.6e6, y1 = gen$largo)
  )

  g <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = brazo, ggplot2::aes(x = .data$x, xend = .data$x, y = .data$y0 / mb, yend = .data$y1 / mb),
                          colour = "#dcdde3", linewidth = 4.2, lineend = "round") +
    ggplot2::geom_segment(data = brazo, ggplot2::aes(x = .data$x, xend = .data$x, y = .data$y0 / mb, yend = .data$y1 / mb),
                          colour = "#eef0f4", linewidth = 2.6, lineend = "round") +
    ggplot2::geom_point(data = gen, ggplot2::aes(x = .data$x, y = .data$centromero / mb), colour = "#a32626",
                        size = 1.1) +
    ggplot2::geom_segment(data = panel, ggplot2::aes(x = .data$x - 0.36, xend = .data$x - 0.14,
                                                     y = .data$medio / mb, yend = .data$medio / mb,
                                                     colour = .data$categoria_nombre),
                          linewidth = 0.55, na.rm = TRUE)
  escala_col <- ggplot2::scale_colour_manual(values = stats::setNames(unname(.colores_iuis), .nombres_iuis),
                                             na.value = "#c9c8c2", name = "Categor\u00eda IUIS",
                                             breaks = intersect(.nombres_iuis, as.character(panel$categoria_nombre)),
                                             guide = ggplot2::guide_legend(override.aes = list(linewidth = 2),
                                                                           nrow = 4))
  if (!is.null(marcas)) {
    marcas$x <- as.integer(.cromosomes_factor(marcas$chr))
    top <- marcas[!duplicated(marcas$gen), ]
    top <- utils::head(top, etiquetas)
    if (inherits(x, "iei_resultado")) {
      marcas$grupo <- factor(marcas$grupo, levels = .niveles_categoria)
      relleno <- ggplot2::scale_fill_manual(values = .colores_categoria, name = "Prioridad", drop = FALSE)
    } else {
      relleno <- ggplot2::scale_fill_manual(values = .colores_herencia, name = "Herencia", na.value = "#a3a29b")
    }
    g <- g +
      ggplot2::geom_point(data = marcas, ggplot2::aes(x = .data$x + 0.25, y = .data$pos / mb, fill = .data$grupo),
                          shape = 21, colour = "white", size = 2.6, stroke = 0.5) +
      ggplot2::geom_text(data = top, ggplot2::aes(x = .data$x + 0.42, y = .data$pos / mb, label = .data$gen),
                         hjust = 0, size = 2.7, fontface = "bold", colour = .tinta[["primaria"]]) +
      relleno
  }
  if (is.null(titulo)) {
    titulo <- if (inherits(x, "iei_resultado")) "Variantes del caso en los cromosomas"
              else if (inherits(x, "iei_busqueda")) "Genes candidatos en los cromosomas"
              else "Genes del panel en los cromosomas"
  }
  g + escala_col +
    ggplot2::scale_x_continuous(breaks = gen$x, labels = as.character(gen$chr), expand = ggplot2::expansion(add = 0.7)) +
    ggplot2::scale_y_reverse(breaks = seq(0, 250, 50), labels = function(v) paste(v, "Mb")) +
    ggplot2::labs(x = NULL, y = NULL, title = titulo,
                  subtitle = sprintf("%d genes del panel, coloreados por categor\u00eda IUIS. Genoma %s.",
                                     length(unique(panel$gen)), build)) +
    .tema_iei(leyenda = "bottom") +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), legend.box = "vertical",
                   legend.box.just = "left",
                   axis.text.x = ggplot2::element_text(face = "bold", colour = .marca, size = 10))
}
