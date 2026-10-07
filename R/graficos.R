# Graficos SVG para el reporte. Sin dependencias: cada funcion devuelve texto.

.colores_componente <- c(
  p_clinvar = "#7a5af8", p_impacto = "#e04f16", p_frecuencia = "#0ba5ec",
  p_herencia = "#12b76a", p_fenotipo = "#ee46bc"
)
.nombres_componente <- c(
  p_clinvar = "ClinVar", p_impacto = "Impacto", p_frecuencia = "Frecuencia",
  p_herencia = "Herencia", p_fenotipo = "Fenotipo"
)
.colores_dominio <- c("#d1e9ff", "#d1fadf", "#fef0c7", "#fce7f6", "#ebe9fe", "#ccfbef")

.num <- function(x) formatC(x, format = "f", digits = 1)

.svg_histograma <- function(h, titulo, etiqueta_x, color = "#475467", ancho = 320, alto = 170) {
  n <- h$n
  if (sum(n) == 0) {
    return(sprintf('<figure class="graf"><figcaption>%s</figcaption><p class="sub">Sin datos</p></figure>',
                   titulo))
  }
  m_izq <- 34; m_inf <- 34; m_sup <- 10; m_der <- 8
  w <- ancho - m_izq - m_der; hh <- alto - m_inf - m_sup
  bw <- w / length(n)
  ymax <- max(n)
  barras <- paste(sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="%s"/>',
                          .num(m_izq + (seq_along(n) - 1) * bw + 1),
                          .num(m_sup + hh - n / ymax * hh), .num(max(bw - 2, 1)),
                          .num(n / ymax * hh), color), collapse = "")
  cortes <- h$cortes
  lab <- ifelse(is.infinite(cortes), "", formatC(cortes, format = "g"))
  paso <- max(1, floor(length(cortes) / 5))
  idx <- seq(1, length(cortes), by = paso)
  marcas_x <- paste(sprintf('<text x="%s" y="%s" text-anchor="middle">%s</text>',
                            .num(m_izq + (idx - 1) * bw), alto - m_inf + 14, lab[idx]),
                    collapse = "")
  sprintf(paste0('<figure class="graf"><figcaption>%s</figcaption>',
                 '<svg viewBox="0 0 %d %d" role="img" aria-label="%s" class="svg">',
                 '<line x1="%d" y1="%s" x2="%d" y2="%s" class="eje"/>%s%s',
                 '<text x="%d" y="%d" class="eti">%s</text>',
                 '<text x="%d" y="%d" text-anchor="end">%s</text>',
                 '<text x="%d" y="%d" text-anchor="middle" class="eti">%s</text></svg></figure>'),
          titulo, ancho, alto, titulo,
          m_izq, .num(m_sup + hh), ancho - m_der, .num(m_sup + hh),
          barras, marcas_x, 2, m_sup + 10, "", m_izq - 4, m_sup + 10, ymax,
          m_izq + w %/% 2, alto - 4, etiqueta_x)
}

.svg_puntaje <- function(d, ancho = 760) {
  comp <- names(.colores_componente)
  if (nrow(d) == 0 || !all(comp %in% names(d))) return("")
  d <- utils::head(d, 15)
  fila_h <- 28; m_izq <- 190; m_der <- 40; m_sup <- 34
  pos_max <- max(rowSums(pmax(as.matrix(d[, comp]), 0)), 1)
  neg_min <- min(rowSums(pmin(as.matrix(d[, comp]), 0)), 0)
  rango <- pos_max - neg_min
  w <- ancho - m_izq - m_der
  x0 <- m_izq + (-neg_min) / rango * w
  escala <- w / rango
  alto <- m_sup + nrow(d) * fila_h + 16

  filas <- vapply(seq_len(nrow(d)), function(i) {
    y <- m_sup + (i - 1) * fila_h
    etiqueta <- paste0(.esc(d$gen[i]), " ", .esc(sub("^[^:]+:", "", .celda(d, "hgvs_p")[i])))
    acum_pos <- 0; acum_neg <- 0
    segs <- character(0)
    for (cn in comp) {
      v <- d[[cn]][i]
      if (is.na(v) || v == 0) next
      if (v > 0) { xa <- x0 + acum_pos * escala; acum_pos <- acum_pos + v }
      else { acum_neg <- acum_neg + v; xa <- x0 + acum_neg * escala }
      segs <- c(segs, sprintf('<rect x="%s" y="%s" width="%s" height="%d" fill="%s"><title>%s %+g</title></rect>',
                              .num(xa), .num(y + 6), .num(abs(v) * escala), fila_h - 12,
                              .colores_componente[[cn]], .nombres_componente[[cn]], v))
    }
    sprintf(paste0('<text x="%d" y="%s" text-anchor="end" class="gen">%s</text>%s',
                   '<text x="%s" y="%s" class="total">%s</text>'),
            m_izq - 8, .num(y + fila_h / 2 + 4), etiqueta, paste(segs, collapse = ""),
            .num(x0 + acum_pos * escala + 6), .num(y + fila_h / 2 + 4), .esc(d$puntaje[i]))
  }, character(1))

  leyenda <- paste(vapply(seq_along(comp), function(k) {
    sprintf('<rect x="%d" y="8" width="12" height="12" rx="2" fill="%s"/><text x="%d" y="18">%s</text>',
            m_izq + (k - 1) * 110, .colores_componente[[k]], m_izq + (k - 1) * 110 + 16,
            .nombres_componente[[k]])
  }, character(1)), collapse = "")

  sprintf(paste0('<figure class="graf ancho"><svg viewBox="0 0 %d %d" role="img" ',
                 'aria-label="Desglose del puntaje por variante" class="svg">%s',
                 '<line x1="%s" y1="%d" x2="%s" y2="%d" class="cero"/>%s</svg>',
                 '<figcaption>Cada barra suma los puntos de cada componente. A la izquierda de la linea, ',
                 'los puntos que restan.</figcaption></figure>'),
          ancho, alto, leyenda, .num(x0), m_sup - 4, .num(x0), alto - 10,
          paste(filas, collapse = ""))
}

.svg_embudo <- function(b, ancho = 760) {
  if (is.null(b) || nrow(b) == 0) return("")
  fila_h <- 34; alto <- nrow(b) * fila_h + 10
  w_max <- ancho - 360
  # Escala logaritmica: los pasos suelen ir de miles a unas pocas variantes.
  base <- log10(max(b$filas[1], 1) + 1)
  filas <- vapply(seq_len(nrow(b)), function(i) {
    w <- max(log10(b$filas[i] + 1) / max(base, 1) * w_max, 2)
    x <- 200 + (w_max - w) / 2
    y <- 4 + (i - 1) * fila_h
    sprintf(paste0('<text x="190" y="%d" text-anchor="end" class="gen">%s</text>',
                   '<rect x="%s" y="%d" width="%s" height="%d" rx="4" fill="#475467" opacity="%s"/>',
                   '<text x="%d" y="%d" class="total">%s%s</text>'),
            y + 20, .esc(b$paso[i]), .num(x), y + 4, .num(w), fila_h - 10,
            .num(1 - 0.08 * (i - 1)), 200 + w_max + 12, y + 20,
            format(b$filas[i], big.mark = " "),
            if (i > 1 && b$descartadas[i] > 0) sprintf(' <tspan class="sub">(-%s)</tspan>',
                                                     format(b$descartadas[i], big.mark = " ")) else "")
  }, character(1))
  sprintf(paste0('<figure class="graf ancho"><svg viewBox="0 0 %d %d" role="img" aria-label="Embudo del filtrado" class="svg">%s</svg>',
                 '<figcaption>Ancho de las barras en escala logar&iacute;tmica. Entre par&eacute;ntesis, las filas descartadas en cada paso.</figcaption></figure>'),
          ancho, alto, paste(filas, collapse = ""))
}

.svg_proteina <- function(info, marcas, ancho = 760) {
  if (is.null(info)) return('<p class="sub">Sin informacion de la proteina.</p>')
  largo <- info$largo; dom <- info$dominios
  m_izq <- 20; m_der <- 20; w <- ancho - m_izq - m_der
  X <- function(p) m_izq + (p - 1) / max(largo - 1, 1) * w
  y_eje <- 96; alto <- 150
  doms <- if (nrow(dom) == 0) "" else paste(vapply(seq_len(nrow(dom)), function(i) {
    x1 <- X(dom$inicio[i]); x2 <- X(dom$fin[i])
    color <- .colores_dominio[(i - 1) %% length(.colores_dominio) + 1]
    nombre <- dom$nombre[i]
    max_c <- floor((x2 - x1) / 6.5)
    corto <- if (nchar(nombre) > max_c) if (max_c > 3) paste0(substr(nombre, 1, max_c - 1), ".") else "" else nombre
    sprintf(paste0('<rect x="%s" y="%d" width="%s" height="22" rx="4" fill="%s" stroke="#98a2b3" stroke-width=".5">',
                   '<title>%s (%g-%g)</title></rect><text x="%s" y="%d" text-anchor="middle" class="dom">%s</text>'),
            .num(x1), y_eje - 11, .num(max(x2 - x1, 2)), color, .esc(nombre), dom$inicio[i], dom$fin[i],
            .num((x1 + x2) / 2), y_eje + 4, .esc(corto))
  }, character(1)), collapse = "")
  marcas <- marcas[!is.na(marcas$pos), ]
  colores <- c(alta = "#b42318", media = "#b54708", baja = "#475467")
  pins <- if (nrow(marcas) == 0) "" else paste(vapply(seq_len(nrow(marcas)), function(i) {
    x <- X(marcas$pos[i]); y_top <- 30 + (i - 1) %% 2 * 16
    col <- colores[[marcas$categoria[i]]] %||% "#475467"
    sprintf(paste0('<line x1="%s" y1="%d" x2="%s" y2="%d" stroke="%s" stroke-width="1.5"/>',
                   '<circle cx="%s" cy="%d" r="6" fill="%s"/>',
                   '<text x="%s" y="%d" text-anchor="%s" class="pin">%s</text>'),
            .num(x), y_top, .num(x), y_eje - 11, col, .num(x), y_top, col,
            .num(x), y_top - 10, if (x < 80) "start" else if (x > ancho - 80) "end" else "middle",
            .esc(marcas$etiqueta[i]))
  }, character(1)), collapse = "")
  sprintf(paste0('<svg viewBox="0 0 %d %d" role="img" aria-label="Diagrama de la proteina" class="svg prot">',
                 '<rect x="%d" y="%d" width="%d" height="8" rx="4" fill="#e4e7ec"/>%s%s',
                 '<text x="%d" y="%d" class="eti">1</text>',
                 '<text x="%d" y="%d" text-anchor="end" class="eti">%d aa</text></svg>'),
          ancho, alto, m_izq, y_eje - 4, w, doms, pins,
          m_izq, y_eje + 32, ancho - m_der, y_eje + 32, largo)
}

.tabla_fenotipo <- function(d, fenotipo, hpo) {
  if (is.null(fenotipo) || nrow(d) == 0) return("")
  genes <- unique(d$gen)
  terminos <- split(hpo$hpo_id, hpo$gen)
  nombres <- hpo$hpo_nombre[match(fenotipo, hpo$hpo_id)]
  enc <- c("Gen", sprintf('<span title="%s">%s</span><br><span class="sub">%s</span>',
                          .esc(nombres), .esc(fenotipo), .esc(ifelse(is.na(nombres), "", nombres))),
           "Coincidencias")
  filas <- vapply(genes, function(g) {
    tiene <- fenotipo %in% terminos[[g]]
    sprintf('<tr><td><b>%s</b></td>%s<td class="num">%d de %d</td></tr>', .esc(g),
            paste(ifelse(tiene, '<td class="si">&#9679;</td>', '<td class="no">&middot;</td>'),
                  collapse = ""), sum(tiene), length(fenotipo))
  }, character(1))
  .tabla_html(enc, filas, clase = "matriz")
}
