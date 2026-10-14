# Objeto iei_resultado: la tabla priorizada con metodos propios.
# Sigue siendo un tibble, asi que funciona con dplyr, corchetes y ggplot2.

.columnas_resultado <- c("muestra", "gen", "puntaje", "categoria")

.como_resultado <- function(x) {
  if (!inherits(x, "iei_resultado")) class(x) <- c("iei_resultado", class(x))
  x
}

.es_resultado_completo <- function(x) all(.columnas_resultado %in% names(x))

.niveles_categoria <- c("alta", "media", "baja")

#' @export
print.iei_resultado <- function(x, n = 10, ...) {
  if (!.es_resultado_completo(x)) {
    class(x) <- setdiff(class(x), "iei_resultado")
    return(print(x, n = n, ...))
  }
  n_cat <- table(factor(x$categoria, levels = .niveles_categoria))
  cli::cli_text("{.strong Resultado de ieiprio} {cli::symbol$bullet} ",
                "{length(unique(x$muestra))} muestra{?s}, {nrow(x)} variante{?s}")
  cli::cli_text("Genoma {attr(x, 'build') %||% '?'} {cli::symbol$bullet} ",
                "Ensembl VEP {attr(x, 'vep_release') %||% '?'} {cli::symbol$bullet} ",
                "prioridad alta {n_cat[['alta']]}, media {n_cat[['media']]}, baja {n_cat[['baja']]}")
  cols <- intersect(c("muestra", "gen", "hgvs_c", "hgvs_p", "genotipo", "puntaje",
                      "categoria", "acmg"), names(x))
  vista <- tibble::as_tibble(as.data.frame(x)[, cols, drop = FALSE])
  if ("hgvs_c" %in% cols) vista$hgvs_c <- sub("^[^:]+:", "", vista$hgvs_c)
  if ("hgvs_p" %in% cols) vista$hgvs_p <- sub("^[^:]+:", "", vista$hgvs_p)
  print(utils::head(vista, n))
  if (nrow(x) > n) cli::cli_text("{.emph ... y {nrow(x) - n} m\u00e1s. Usa {.fn iei_tabla} para ver todo.}")
  invisible(x)
}

#' Resumen de un resultado
#'
#' @param object Resultado de [iei_priorizar()] o [iei_analizar()].
#' @param ... No se usa.
#' @return Un objeto `summary.iei_resultado`, una lista con las tablas
#'   `por_categoria`, `por_gen`, `por_muestra`, `acmg` y `bitacora`.
#' @export
summary.iei_resultado <- function(object, ...) {
  x <- as.data.frame(object)
  cat_f <- factor(x$categoria, levels = .niveles_categoria)
  por_gen <- if (nrow(x) == 0) tibble::tibble(gen = character(0), variantes = integer(0),
                                               mejor_puntaje = numeric(0), mejor_categoria = character(0)) else {
    g <- split(x, x$gen)
    tibble::tibble(
      gen = names(g),
      variantes = vapply(g, nrow, integer(1)),
      mejor_puntaje = vapply(g, function(d) max(d$puntaje), numeric(1)),
      mejor_categoria = vapply(g, function(d) d$categoria[which.max(d$puntaje)], character(1))
    )
  }
  por_gen <- por_gen[order(-por_gen$mejor_puntaje), ]
  por_muestra <- if (nrow(x) == 0) tibble::tibble(muestra = character(0)) else {
    m <- split(x, x$muestra)
    tibble::tibble(
      muestra = names(m),
      alta = vapply(m, function(d) sum(d$categoria == "alta"), integer(1)),
      media = vapply(m, function(d) sum(d$categoria == "media"), integer(1)),
      baja = vapply(m, function(d) sum(d$categoria == "baja"), integer(1)),
      principal = vapply(m, function(d) d$gen[which.max(d$puntaje)], character(1))
    )
  }
  acmg <- if ("acmg" %in% names(x)) {
    crit <- trimws(unlist(strsplit(stats::na.omit(x$acmg), ",")))
    t <- sort(table(crit), decreasing = TRUE)
    tibble::tibble(criterio = names(t), variantes = as.integer(t))
  } else tibble::tibble(criterio = character(0), variantes = integer(0))

  structure(list(
    por_categoria = table(cat_f, dnn = "categoria"),
    por_gen = por_gen, por_muestra = por_muestra, acmg = acmg,
    bitacora = attr(object, "ieiprio_bitacora"),
    build = attr(object, "build"), vep_release = attr(object, "vep_release")
  ), class = "summary.iei_resultado")
}

#' @export
print.summary.iei_resultado <- function(x, ...) {
  cli::cli_h2("Resumen de ieiprio")
  cli::cli_text("Genoma {x$build %||% '?'}, Ensembl VEP {x$vep_release %||% '?'}")
  cli::cli_h3("Variantes por prioridad")
  print(x$por_categoria)
  cli::cli_h3("Por muestra")
  print(x$por_muestra)
  cli::cli_h3("Genes con mayor puntaje")
  print(utils::head(x$por_gen, 10))
  if (nrow(x$acmg) > 0) {
    cli::cli_h3("Criterios ACMG sugeridos")
    print(x$acmg)
  }
  if (!is.null(x$bitacora)) {
    cli::cli_h3("Filtrado")
    print(x$bitacora)
  }
  invisible(x)
}

#' Graficar un resultado
#'
#' Devuelve un grafico de ggplot2 que puedes modificar como cualquier otro
#' (agregar temas, cambiar colores, combinar con patchwork). Con
#' `interactivo = TRUE` lo convierte con plotly.
#'
#' @param x Resultado de [iei_priorizar()] o [iei_analizar()].
#' @param tipo Que graficar.
#'   * `"puntaje"` el desglose del puntaje de cada variante por componente.
#'   * `"categorias"` cuantas variantes hay en cada prioridad por muestra.
#'   * `"genes"` los genes con mayor puntaje.
#'   * `"predictores"` REVEL contra AlphaMissense, con los umbrales usados.
#'   * `"mapa"` el mejor puntaje de cada gen en cada muestra, para comparar
#'     varias muestras de una cohorte.
#'   * `"oncoplot"` el tipo de variante de cada gen en cada muestra, al estilo
#'     de maftools, con el total por gen al costado.
#'   * `"cariotipo"` las variantes sobre los cromosomas ([iei_cariograma()]).
#'   * `"circos"` las variantes sobre el genoma en circulo ([iei_circos()]);
#'     dibuja con graficos base y no devuelve un ggplot.
#' @param n Numero maximo de variantes o genes a mostrar.
#' @param interactivo Si `TRUE`, devuelve un grafico de plotly.
#' @param ... No se usa.
#' @return Un objeto ggplot, o plotly si `interactivo = TRUE`.
#' @export
#' @examples
#' \dontrun{
#' res <- iei_analizar(iei_ejemplo_vcf(), sexo = "M", carpeta = NULL)
#' plot(res)
#' plot(res, tipo = "cariotipo")
#' plot(res, tipo = "predictores", interactivo = TRUE)
#' }
plot.iei_resultado <- function(x, tipo = c("puntaje", "categorias", "genes", "predictores", "mapa",
                                           "oncoplot", "cariotipo", "circos"),
                               n = 20, interactivo = FALSE, ...) {
  tipo <- match.arg(tipo)
  if (tipo == "circos") return(iei_circos(x, etiquetas = n))
  rlang::check_installed("ggplot2", reason = "para graficar resultados.")
  d <- as.data.frame(x)
  if (nrow(d) == 0) cli::cli_abort("El resultado no tiene variantes para graficar.")
  g <- switch(tipo,
    puntaje = .grafico_puntaje(d, n),
    categorias = .grafico_categorias(d),
    genes = .grafico_genes(d, n),
    predictores = .grafico_predictores(d),
    mapa = .grafico_mapa(d, n),
    oncoplot = .grafico_oncoplot(d, n),
    cariotipo = iei_cariograma(x, etiquetas = n)
  )
  if (interactivo) {
    rlang::check_installed("plotly", reason = "para graficos interactivos.")
    return(suppressWarnings(plotly::ggplotly(g, tooltip = "text")))
  }
  g
}

.etiqueta_variante <- function(d) {
  p <- sub("^[^:]+:", "", .col(d, "hgvs_p", NA_character_))
  c_ <- sub("^[^:]+:", "", .col(d, "hgvs_c", NA_character_))
  cambio <- ifelse(!is.na(p), p, ifelse(!is.na(c_), c_, paste0(d$chr %||% "", ":", d$pos %||% "")))
  paste0(d$gen, " ", cambio, if (length(unique(d$muestra)) > 1) paste0(" (", d$muestra, ")") else "")
}

.grafico_puntaje <- function(d, n) {
  comp <- intersect(names(.colores_componente), names(d))
  d <- utils::head(d[order(-d$puntaje), ], n)
  d$variante <- factor(.etiqueta_variante(d), levels = rev(unique(.etiqueta_variante(d))))
  largo <- do.call(rbind, lapply(comp, function(cn) data.frame(
    variante = d$variante, componente = .nombres_componente[[cn]], puntos = d[[cn]],
    stringsAsFactors = FALSE)))
  largo <- largo[!is.na(largo$puntos) & largo$puntos != 0, ]
  largo$componente <- factor(largo$componente, levels = unname(.nombres_componente[comp]))
  largo$text <- sprintf("%s\n%s %+g", largo$variante, largo$componente, largo$puntos)
  colores <- stats::setNames(unname(.colores_componente[comp]), unname(.nombres_componente[comp]))
  total <- data.frame(variante = d$variante, puntaje = d$puntaje, categoria = d$categoria,
                      tope = pmax(0, rowSums(pmax(as.matrix(d[, comp, drop = FALSE]), 0), na.rm = TRUE)))
  ggplot2::ggplot(largo, ggplot2::aes(x = .data$puntos, y = .data$variante, text = .data$text)) +
    ggplot2::annotate("rect", xmin = 9, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#fcecec", alpha = 0.6) +
    ggplot2::annotate("rect", xmin = 5, xmax = 9, ymin = -Inf, ymax = Inf, fill = "#fdf3e7", alpha = 0.6) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$componente), width = 0.68,
                      position = ggplot2::position_stack(reverse = TRUE), colour = "white", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = .tinta[["secundaria"]]) +
    ggplot2::geom_label(data = total, ggplot2::aes(x = .data$tope, y = .data$variante, label = .data$puntaje,
                                                   colour = .data$categoria),
                        hjust = -0.25, size = 3.3, fontface = "bold", linewidth = 0, fill = "white",
                        inherit.aes = FALSE, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = colores, drop = FALSE) +
    ggplot2::scale_colour_manual(values = .colores_categoria, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
    ggplot2::labs(x = "Puntos", y = NULL, fill = NULL, title = "De d\u00f3nde sale el puntaje de cada variante",
                  subtitle = "Las bandas marcan la prioridad media (5 a 8) y alta (9 o m\u00e1s). A la izquierda del cero, lo que resta.") +
    .tema_iei() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(colour = .tinta[["primaria"]]))
}

.colores_categoria <- c(alta = "#b42318", media = "#d97706", baja = "#98a2b3")

.grafico_categorias <- function(d) {
  t <- as.data.frame(table(muestra = d$muestra,
                           categoria = factor(d$categoria, levels = .niveles_categoria)))
  t$text <- sprintf("%s\n%s: %d", t$muestra, t$categoria, t$Freq)
  ggplot2::ggplot(t, ggplot2::aes(x = .data$Freq, y = .data$muestra, fill = .data$categoria,
                                  text = .data$text)) +
    ggplot2::geom_col(width = 0.6, colour = "white", linewidth = 0.4,
                      position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$Freq > 0, .data$Freq, "")), colour = "white",
                       fontface = "bold", size = 3.4, position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE)) +
    ggplot2::scale_fill_manual(values = .colores_categoria) +
    ggplot2::labs(x = "Variantes", y = NULL, fill = "Prioridad", title = "Variantes por prioridad") +
    .tema_iei() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

.grafico_genes <- function(d, n) {
  s <- summary.iei_resultado(.como_resultado(d))$por_gen
  s <- utils::head(s, n)
  s$gen <- factor(s$gen, levels = rev(s$gen))
  s$text <- sprintf("%s\n%d variante(s), mejor puntaje %s", s$gen, s$variantes, s$mejor_puntaje)
  s$mejor_categoria <- factor(s$mejor_categoria, levels = .niveles_categoria)
  ggplot2::ggplot(s, ggplot2::aes(x = .data$mejor_puntaje, y = .data$gen, text = .data$text)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$mejor_puntaje, yend = .data$gen,
                                       colour = .data$mejor_categoria), linewidth = 1.1) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$mejor_categoria, size = .data$variantes),
                        shape = 21, colour = "white", stroke = 0.8) +
    ggplot2::scale_fill_manual(values = .colores_categoria, name = "Prioridad") +
    ggplot2::scale_colour_manual(values = .colores_categoria, guide = "none") +
    ggplot2::scale_size_continuous(range = c(4, 9), name = "Variantes", breaks = function(l) unique(round(pretty(l, 3)))) +
    ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(fill = "#a3a29b")),
                    fill = ggplot2::guide_legend(override.aes = list(size = 4.5))) +
    ggplot2::labs(x = "Mejor puntaje", y = NULL, title = "Genes con mayor puntaje",
                  subtitle = "El tama\u00f1o del punto es el n\u00famero de variantes del gen.") +
    .tema_iei() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "bold.italic", colour = .tinta[["primaria"]]))
}

.grafico_predictores <- function(d) {
  d$revel <- .col(d, "revel"); d$alphamissense <- .col(d, "alphamissense")
  d <- d[!is.na(d$revel) | !is.na(d$alphamissense), ]
  if (nrow(d) == 0) cli::cli_abort("Ninguna variante tiene REVEL ni AlphaMissense (son solo para missense).")
  d$etiqueta <- .etiqueta_variante(d)
  d$text <- sprintf("%s\nREVEL %s\nAlphaMissense %s", d$etiqueta,
                    format(round(d$revel, 2)), format(round(d$alphamissense, 2)))
  u <- .umbrales_predictores
  ggplot2::ggplot(d, ggplot2::aes(x = .data$revel, y = .data$alphamissense, text = .data$text)) +
    ggplot2::annotate("rect", xmin = u$revel_favor, xmax = 1, ymin = u$am_favor, ymax = 1, fill = "#fcecec") +
    ggplot2::annotate("rect", xmin = 0, xmax = u$revel_contra, ymin = 0, ymax = u$am_contra, fill = "#eaf2fc") +
    ggplot2::annotate("text", x = u$revel_favor, y = 1.035, label = "ambos a favor de da\u00f1o (PP3)", hjust = 0,
                      size = 3, colour = "#a32626", fontface = "bold") +
    ggplot2::annotate("text", x = 0.01, y = 0.01, label = "ambos benignos (BP4)", hjust = 0, vjust = 0,
                      size = 3, colour = "#184f95", fontface = "bold") +
    ggplot2::geom_vline(xintercept = c(u$revel_contra, u$revel_favor), linetype = "22", linewidth = 0.3,
                        colour = .tinta[["tenue"]]) +
    ggplot2::geom_hline(yintercept = c(u$am_contra, u$am_favor), linetype = "22", linewidth = 0.3,
                        colour = .tinta[["tenue"]]) +
    ggplot2::geom_point(ggplot2::aes(fill = factor(.data$categoria, levels = .niveles_categoria)), shape = 21,
                        colour = "white", size = 4.5, stroke = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$gen), hjust = -0.3, size = 3, fontface = "bold.italic",
                       colour = .tinta[["primaria"]]) +
    ggplot2::scale_fill_manual(values = .colores_categoria, name = "Prioridad") +
    ggplot2::coord_equal(xlim = c(0, 1.05), ylim = c(0, 1.05)) +
    ggplot2::labs(x = "REVEL", y = "AlphaMissense", title = "Predictores de patogenicidad",
                  subtitle = "L\u00edneas punteadas, umbrales de apoyo de ClinGen para PP3 y BP4.") +
    .tema_iei()
}

# Oncoplot al estilo de maftools: tipo de variante por gen y muestra.
.clase_consecuencia <- function(cons) {
  cons <- ifelse(is.na(cons), "", cons)
  ifelse(grepl("stop_gained|stop_lost|start_lost", cons), "sin sentido",
  ifelse(grepl("frameshift|inframe", cons), "frameshift o indel",
  ifelse(grepl("splice", cons), "splicing",
  ifelse(grepl("missense", cons), "missense",
  ifelse(grepl("synonymous", cons), "sinonima", "otra")))))
}

.grafico_oncoplot <- function(d, n) {
  d$clase <- factor(.clase_consecuencia(.col(d, "consecuencia", NA_character_)), levels = names(.colores_tipo))
  orden <- names(sort(tapply(d$puntaje, d$gen, max), decreasing = TRUE))
  orden <- utils::head(orden, n)
  d <- d[d$gen %in% orden, ]
  d <- d[order(d$gen, d$muestra, -d$puntaje), ]
  celda <- d[!duplicated(paste(d$gen, d$muestra)), ]
  celda$gen <- factor(celda$gen, levels = rev(orden))
  n_gen <- as.data.frame(table(gen = factor(d$gen, levels = rev(orden)), clase = d$clase))
  muestras_con <- tapply(celda$muestra, celda$gen, function(m) length(unique(m)))
  total_muestras <- length(unique(d$muestra))
  celda$text <- sprintf("%s en %s\n%s, puntaje %s", celda$gen, celda$muestra, celda$clase, celda$puntaje)
  rotulo <- data.frame(gen = factor(names(muestras_con), levels = levels(celda$gen)),
                       pct = sprintf("%d %%", round(100 * as.numeric(muestras_con) / total_muestras)))
  rotulo <- rotulo[!is.na(rotulo$gen), ]
  base <- ggplot2::ggplot(celda, ggplot2::aes(x = .data$muestra, y = .data$gen, text = .data$text)) +
    ggplot2::geom_tile(data = expand.grid(muestra = unique(d$muestra), gen = levels(celda$gen)),
                       fill = "#f1f1ee", colour = "white", linewidth = 1.2, inherit.aes = FALSE,
                       mapping = ggplot2::aes(x = .data$muestra, y = .data$gen)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$clase), colour = "white", linewidth = 1.2,
                       height = 0.82) +
    ggplot2::scale_fill_manual(values = .colores_tipo, drop = TRUE, name = NULL) +
    ggplot2::scale_y_discrete(position = "left") +
    ggplot2::labs(x = NULL, y = NULL) +
    .tema_iei() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "bold.italic", colour = .tinta[["primaria"]]),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    return(base + ggplot2::labs(title = "Variantes por gen y muestra"))
  }
  lateral <- ggplot2::ggplot(n_gen, ggplot2::aes(x = .data$Freq, y = .data$gen, fill = .data$clase)) +
    ggplot2::geom_col(width = 0.75, show.legend = FALSE) +
    ggplot2::geom_text(data = rotulo, ggplot2::aes(x = 0, y = .data$gen, label = .data$pct), hjust = 1.2,
                       size = 3, colour = .tinta[["secundaria"]], inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = .colores_tipo, drop = TRUE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.45, 0.05)),
                                breaks = function(l) { b <- pretty(c(0, l[2])); b[b >= 0 & b == round(b)] }) +
    ggplot2::labs(x = "Variantes", y = NULL) +
    .tema_iei() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
  patchwork::wrap_plots(base, lateral, widths = c(4, 1.3)) +
    patchwork::plot_annotation(title = "Variantes por gen y muestra",
                               subtitle = "Color, tipo de variante. A la derecha, variantes por gen y porcentaje de muestras con alguna.",
                               theme = .tema_iei())
}

#' Tabla simple de un resultado
#'
#' Devuelve el resultado como un tibble comun, sin la clase ni los metadatos
#' de ieiprio, para exportarlo o trabajarlo con otros paquetes.
#'
#' @param x Resultado de [iei_priorizar()] o [iei_analizar()].
#' @param columnas `"principales"` para las columnas mas usadas o `"todas"`.
#' @return Un tibble.
#' @export
iei_tabla <- function(x, columnas = c("principales", "todas")) {
  columnas <- match.arg(columnas)
  d <- as.data.frame(x)
  if (columnas == "principales") {
    keep <- intersect(c("muestra", "gen", "chr", "pos", "ref", "alt", "genotipo", "hgvs_c",
                        "hgvs_p", "consecuencia", "af_gnomad", "clinvar", "revel",
                        "alphamissense", "spliceai", "puntaje", "categoria", "acmg",
                        "acmg_puntos", "acmg_clase", "razones", "nota"), names(d))
    d <- d[, keep, drop = FALSE]
  }
  tibble::as_tibble(d)
}

#' Metadatos de un analisis
#'
#' Reune en una lista todo lo que se registro durante el analisis: build del
#' genoma, release de Ensembl, archivo VCF y su huella, sexo declarado,
#' parametros de filtrado, pesos, fenotipo usado, bitacora y control de
#' calidad.
#'
#' @param x Resultado de cualquier etapa del flujo.
#' @return Una lista con nombres.
#' @export
iei_metadatos <- function(x) {
  a <- .attrs_ieiprio(x)
  names(a) <- sub("^ieiprio_", "", names(a))
  a
}
