# Graficos circulares con circlize.

.necesita_circlize <- function() rlang::check_installed("circlize", reason = "para los gr\u00e1ficos circulares.")

.alfa <- function(col, a) grDevices::adjustcolor(col, alpha.f = a)

#' Mapa circular del panel en el genoma
#'
#' Dibuja los cromosomas en circulo con los genes del panel coloreados por
#' categoria IUIS y, si se le pasa un resultado o una busqueda, marca encima
#' las variantes o los genes candidatos. Esta inspirado en los graficos de
#' circlize y sirve para ver de un vistazo donde estan los genes del panel,
#' en que categoria cae cada candidato y como se reparten las variantes de un
#' caso o de una cohorte.
#'
#' * Sin `x`, muestra el panel completo con la densidad de genes por region.
#' * Con un resultado de [iei_priorizar()], cada variante es un punto cuya
#'   distancia al centro es su puntaje y cuyo color es su prioridad. Se
#'   rotulan los genes de las mejores variantes.
#' * Con una busqueda de [iei_buscar_fenotipo()], cada gen candidato es una
#'   barra tan alta como su similitud, y los candidatos que comparten
#'   categoria IUIS se unen por dentro.
#'
#' @param x `NULL`, un resultado de [iei_priorizar()] o una busqueda de
#'   [iei_buscar_fenotipo()].
#' @param build Version del genoma.
#' @param etiquetas Cuantos genes rotular.
#' @param titulo Titulo del grafico; `NULL` usa uno por defecto.
#' @return No devuelve nada; dibuja en el dispositivo grafico activo.
#' @seealso [iei_cariograma()] para la misma informacion en cromosomas
#'   lineales.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("circlize", quietly = TRUE)) {
#'   iei_circos()
#'   b <- iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente"), n = 15)
#'   iei_circos(b)
#' }
#' }
iei_circos <- function(x = NULL, build = c("GRCh38", "GRCh37"), etiquetas = 15, titulo = NULL) {
  .necesita_circlize()
  build <- match.arg(build)
  if (inherits(x, "iei_resultado")) build <- attr(x, "build") %||% build
  gen <- .genoma(build)
  panel <- .genes_panel_posicion(build)
  marcas <- .marcas_genoma(x, build)

  lado <- .lienzo_ancho()
  op <- graphics::par(mar = if (lado) c(1, 1, 3, 17) else c(5.5, 1, 3, 1), bg = "white", xpd = NA)
  on.exit({ circlize::circos.clear(); graphics::par(op) }, add = TRUE)
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = c(rep(1.1, 23), 7),
                       cell.padding = c(0, 0, 0, 0), track.margin = c(0.004, 0.004),
                       points.overflow.warning = FALSE)
  cyto <- data.frame(chr = paste0("chr", gen$chr), inicio = 0, fin = gen$largo)
  circlize::circos.genomicInitialize(cyto, plotType = NULL)

  # Rotulos por fuera.
  if (!is.null(marcas) && etiquetas > 0) {
    top <- marcas[!duplicated(marcas$gen), ]
    top <- utils::head(top, etiquetas)
    bed <- data.frame(chr = paste0("chr", top$chr), inicio = top$pos, fin = top$pos + 1, gen = top$gen)
    col_et <- if (inherits(x, "iei_resultado")) .colores_categoria[top$grupo] else unname(.colores_iuis[panel$categoria[match(top$gen, panel$gen)]])
    col_et[is.na(col_et)] <- .tinta[["secundaria"]]
    circlize::circos.genomicLabels(bed, labels.column = 4, side = "outside", cex = 0.62, col = col_et,
                                   line_col = .alfa(col_et, 0.8), connection_height = circlize::mm_h(4),
                                   labels_height = circlize::cm_h(1.2), padding = 0.4)
  }

  # Cromosomas.
  circlize::circos.track(ylim = c(0, 1), track.height = 0.055, bg.border = NA, panel.fun = function(xx, yy) {
    chr <- sub("^chr", "", circlize::CELL_META$sector.index)
    i <- match(chr, .cromosomas)
    fondo <- if (i %% 2 == 1) "#2b3a55" else "#4b5d7c"
    cen <- gen$centromero[i]
    circlize::circos.rect(0, 0.15, cen - 1.5e6, 0.85, col = fondo, border = NA)
    circlize::circos.rect(cen + 1.5e6, 0.15, gen$largo[i], 0.85, col = fondo, border = NA)
    circlize::circos.rect(cen - 1.5e6, 0.38, cen + 1.5e6, 0.62, col = "#a32626", border = NA)
    circlize::circos.text(circlize::CELL_META$xcenter, 0.5, chr, cex = 0.55, col = "white", font = 2,
                          facing = "inside", niceFacing = TRUE)
  })

  # Genes del panel, una raya por gen con el color de su categoria.
  circlize::circos.track(ylim = c(0, 1), track.height = 0.07, bg.border = NA, bg.col = "#f6f6f4",
                         panel.fun = function(xx, yy) {
    chr <- sub("^chr", "", circlize::CELL_META$sector.index)
    g <- panel[panel$chr == chr, ]
    if (nrow(g)) circlize::circos.segments(g$medio, 0.08, g$medio, 0.92, col = g$color, lwd = 1.3)
  })

  if (inherits(x, "iei_resultado")) {
    ymax <- max(10, marcas$valor, na.rm = TRUE)
    ymin <- min(0, marcas$valor, na.rm = TRUE)
    circlize::circos.track(ylim = c(ymin, ymax), track.height = 0.28, bg.border = NA, bg.col = "#fbfbfa",
                           panel.fun = function(xx, yy) {
      for (h in c(5, 9)) circlize::circos.lines(c(0, circlize::CELL_META$xlim[2]), c(h, h),
                                                col = if (h == 9) .alfa("#b42318", 0.35) else .alfa("#b54708", 0.25),
                                                lty = 2, lwd = 0.7)
      chr <- sub("^chr", "", circlize::CELL_META$sector.index)
      m <- marcas[as.character(marcas$chr) == chr, ]
      if (nrow(m)) {
        col <- .colores_categoria[m$grupo]; col[is.na(col)] <- "#98a2b3"
        circlize::circos.segments(m$pos, ymin, m$pos, m$valor, col = .alfa(col, 0.35), lwd = 0.8)
        circlize::circos.points(m$pos, m$valor, pch = 21, bg = .alfa(col, 0.9), col = "white",
                                cex = ifelse(m$grupo == "alta", 1.25, 0.8))
      }
    })
  } else if (inherits(x, "iei_busqueda")) {
    circlize::circos.track(ylim = c(0, 1), track.height = 0.22, bg.border = NA, bg.col = "#fbfbfa",
                           panel.fun = function(xx, yy) {
      chr <- sub("^chr", "", circlize::CELL_META$sector.index)
      m <- marcas[as.character(marcas$chr) == chr, ]
      if (nrow(m)) {
        col <- .colores_herencia[m$grupo]; col[is.na(col)] <- "#a3a29b"
        ancho <- 2.2e6
        circlize::circos.rect(m$pos - ancho, 0, m$pos + ancho, m$valor, col = .alfa(col, 0.85), border = NA)
        circlize::circos.points(m$pos, m$valor, pch = 21, bg = col, col = "white", cex = 1.1)
      }
    })
    # Lazos entre candidatos de la misma categoria IUIS.
    top <- utils::head(marcas, 30)
    top$categoria <- panel$categoria[match(top$gen, panel$gen)]
    for (cat_i in unique(stats::na.omit(top$categoria))) {
      m <- top[top$categoria %in% cat_i, ]
      if (nrow(m) < 2) next
      pares <- utils::combn(nrow(m), 2)
      for (j in seq_len(ncol(pares))) {
        a <- m[pares[1, j], ]; b <- m[pares[2, j], ]
        circlize::circos.link(paste0("chr", a$chr), a$pos, paste0("chr", b$chr), b$pos,
                              col = .alfa(.colores_iuis[[cat_i]], 0.35), lwd = 1.2, h.ratio = 0.6)
      }
    }
  } else {
    # Densidad de genes del panel cada 5 Mb.
    circlize::circos.track(ylim = c(0, 1), track.height = 0.18, bg.border = NA, panel.fun = function(xx, yy) {
      chr <- sub("^chr", "", circlize::CELL_META$sector.index)
      i <- match(chr, .cromosomas)
      cortes <- seq(0, gen$largo[i] + 10e6, by = 10e6)
      g <- panel[panel$chr == chr, ]
      n <- if (nrow(g)) as.numeric(table(cut(g$medio, cortes))) else rep(0, length(cortes) - 1)
      maximo <- 16
      for (k in seq_along(n)) if (n[k] > 0) {
        circlize::circos.rect(cortes[k], 0, cortes[k + 1], min(1, n[k] / maximo),
                              col = .alfa(.marca_2, 0.35 + 0.6 * min(1, n[k] / maximo)), border = NA)
      }
    })
  }

  if (is.null(titulo)) {
    titulo <- if (inherits(x, "iei_resultado")) "Variantes del caso sobre el genoma"
              else if (inherits(x, "iei_busqueda")) "Genes candidatos por fenotipo"
              else sprintf("Panel de errores innatos de la inmunidad (%d genes)", length(unique(panel$gen)))
  }
  graphics::title(titulo, col.main = .marca, cex.main = 1.15, font.main = 2, line = 1.2)
  cats <- names(.nombres_iuis)
  usr <- graphics::par("usr")
  segunda <- if (inherits(x, "iei_resultado")) {
    list(legend = c("alta", "media", "baja"), pt.bg = .colores_categoria, pch = 21, col = "white",
         pt.cex = 1.4, title = "Prioridad")
  } else if (inherits(x, "iei_busqueda")) {
    h <- intersect(c("AR", "AD", "XL", "AR/AD", "MT", "sin dato"), unique(marcas$grupo))
    list(legend = h, fill = .colores_herencia[h], border = NA, title = "Herencia")
  }
  if (lado) {
    # Lienzo apaisado: las leyendas van en una columna a la derecha del circulo.
    x0 <- usr[2] + 0.04 * diff(usr[1:2])
    y0 <- usr[4]
    if (!is.null(segunda)) {
      caja <- do.call(graphics::legend, c(list(x = x0, y = y0, xjust = 0, yjust = 1, bty = "n", cex = 0.8,
                                               title.font = 2, title.col = .marca, title.adj = 0), segunda))
      y0 <- caja$rect$top - caja$rect$h - 0.04 * diff(usr[3:4])
    }
    graphics::legend(x = x0, y = y0, xjust = 0, yjust = 1, legend = .nombres_iuis, fill = .colores_iuis[cats],
                     border = NA, bty = "n", cex = 0.72, ncol = 1, title = "Categor\u00eda IUIS del gen",
                     title.font = 2, title.adj = 0, text.col = .tinta[["secundaria"]], title.col = .marca)
  } else {
    graphics::legend(x = mean(usr[1:2]), y = usr[3] - 0.03 * diff(usr[3:4]), xjust = 0.5, yjust = 1,
                     legend = .nombres_iuis, fill = .colores_iuis[cats], border = NA, bty = "n",
                     cex = 0.7, ncol = 5, title = "Categor\u00eda IUIS de cada gen del panel", title.font = 2,
                     text.col = .tinta[["secundaria"]], title.col = .marca)
    if (!is.null(segunda)) {
      do.call(graphics::legend, c(list(x = "topright", bty = "n", cex = 0.75, title.font = 2,
                                       title.col = .marca), segunda))
    }
  }
  invisible(NULL)
}

# TRUE si el dispositivo es apaisado (por ejemplo 4:3), para dibujar las
# leyendas al costado del circulo en lugar de abajo.
.lienzo_ancho <- function() {
  d <- graphics::par("din")
  isTRUE(d[1] / d[2] > 1.15)
}

# Diagrama de cuerdas entre los terminos del caso y los genes candidatos.
.grafico_acordes <- function(object, n = 12) {
  .necesita_circlize()
  m <- attr(object, "matriz")
  if (is.null(m)) cli::cli_abort("La b\u00fasqueda no guard\u00f3 la matriz de similitud.")
  d <- utils::head(tibble::as_tibble(unclass(object)), n)
  m <- m[, intersect(d$gen, colnames(m)), drop = FALSE]
  m[m < 0.25] <- 0
  m <- m[rowSums(m) > 0, colSums(m) > 0, drop = FALSE]
  if (length(m) == 0) cli::cli_abort("Ning\u00fan t\u00e9rmino del caso se parece lo suficiente a los genes candidatos.")
  terminos <- .acortar(.hpo_nombre(rownames(m)), 30)
  terminos <- make.unique(ifelse(is.na(terminos), rownames(m), terminos))
  rownames(m) <- terminos
  herencia <- d$herencia[match(colnames(m), d$gen)]
  herencia[is.na(herencia)] <- "sin dato"
  col_gen <- .colores_herencia[herencia]; col_gen[is.na(col_gen)] <- "#a3a29b"
  col_term <- grDevices::colorRampPalette(c("#0b2545", "#1d6fa5", "#6da7ec"))(nrow(m))
  grid_col <- c(stats::setNames(col_term, rownames(m)), stats::setNames(unname(col_gen), colnames(m)))

  lado <- .lienzo_ancho()
  op <- graphics::par(mar = if (lado) c(1, 1, 3, 12) else c(1, 1, 3, 1), bg = "white", xpd = NA)
  on.exit({ circlize::circos.clear(); graphics::par(op) }, add = TRUE)
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.after = c(rep(2, nrow(m) - 1), 10, rep(2, ncol(m) - 1), 10))
  circlize::chordDiagram(m, grid.col = grid_col, transparency = 0.35, annotationTrack = "grid",
                         annotationTrackHeight = circlize::mm_h(3), link.sort = TRUE, link.decreasing = TRUE,
                         preAllocateTracks = list(track.height = circlize::mm_h(22)))
  circlize::circos.track(track.index = 1, bg.border = NA, panel.fun = function(xx, yy) {
    sec <- circlize::CELL_META$sector.index
    es_gen <- sec %in% colnames(m)
    circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], sec,
                          facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                          cex = if (es_gen) 0.75 else 0.6, font = if (es_gen) 2 else 1,
                          col = if (es_gen) .tinta[["primaria"]] else .tinta[["secundaria"]])
  })
  graphics::title("Qu\u00e9 signos del caso explican cada gen", col.main = .marca, cex.main = 1.1, font.main = 2,
                  line = 1.2)
  h <- intersect(c("AR", "AD", "XL", "AR/AD", "MT", "sin dato"), unique(herencia))
  usr <- graphics::par("usr")
  graphics::legend(x = if (lado) usr[2] + 0.06 * diff(usr[1:2]) else "topright", y = if (lado) usr[4],
                   xjust = 0, yjust = 1, legend = h, fill = .colores_herencia[h], border = NA, bty = "n",
                   cex = 0.8, title = "Herencia del gen", title.font = 2, title.adj = 0, title.col = .marca)
  invisible(NULL)
}
