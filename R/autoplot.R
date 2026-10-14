# Graficos con ggplot2 para los objetos del paquete. Todos devuelven un
# objeto ggplot que se puede modificar como cualquier otro.

.necesita_ggplot <- function() rlang::check_installed("ggplot2", reason = "para graficar.")

# ---- iei_resultado ---------------------------------------------------------

#' Graficos con autoplot
#'
#' Metodos de `ggplot2::autoplot()` para los objetos de ieiprio. Devuelven un
#' grafico de ggplot2 que puedes modificar (temas, colores, facetas) o
#' combinar con otros. `plot()` hace lo mismo.
#'
#' * Un resultado de [iei_priorizar()] acepta los tipos de
#'   [plot.iei_resultado()]: `"puntaje"`, `"categorias"`, `"genes"`,
#'   `"predictores"` y `"mapa"` (genes por muestra, util para cohortes).
#' * Una ficha de [iei_gen()] muestra las variantes patogenicas de ClinVar
#'   sobre la proteina con sus dominios (`"proteina"`), cuantas variantes
#'   acumula cada dominio frente al promedio (`"dominios"`) o el reparto de
#'   todas las variantes por clasificacion (`"clinvar"`).
#' * Una busqueda de [iei_buscar_fenotipo()] muestra el ranking de genes
#'   (`"ranking"`), la similitud de cada termino del caso con cada gen
#'   (`"mapa"`) o los candidatos sobre los cromosomas (`"cariotipo"`). Con
#'   `plot()` acepta ademas `"acordes"`, un diagrama de cuerdas entre los
#'   signos del caso y los genes, y `"circos"`, el genoma en circulo; ambos
#'   usan circlize y dibujan con graficos base.
#' * Un cromatograma de [iei_leer_ab1()] muestra las trazas de un tramo.
#' * Una validacion de [iei_validar_fenotipo()] muestra, para cada metodo, el
#'   porcentaje de casos con el gen correcto entre los primeros puestos.
#'
#' @param object Objeto de ieiprio.
#' @param tipo Que graficar; depende del objeto.
#' @param n Numero maximo de elementos a mostrar.
#' @param desde,hasta Bases del cromatograma a mostrar; en una validacion,
#'   `hasta` es el ultimo puesto de la curva.
#' @param ... Otros argumentos; en resultados se pasan a [plot.iei_resultado()].
#' @return Un objeto ggplot.
#' @name autoplot.iei
#' @examples
#' \dontrun{
#' library(ggplot2)
#' autoplot(iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente")))
#' autoplot(iei_gen("BTK"))
#' autoplot(prior, tipo = "mapa")
#' }
NULL

#' @rdname autoplot.iei
#' @exportS3Method ggplot2::autoplot
autoplot.iei_resultado <- function(object, tipo = c("puntaje", "categorias", "genes", "predictores", "mapa"),
                                   n = 20, ...) {
  plot.iei_resultado(object, tipo = match.arg(tipo), n = n, interactivo = FALSE)
}

.grafico_mapa <- function(d, n) {
  mejor <- stats::aggregate(puntaje ~ gen + muestra, d, max)
  orden_gen <- utils::head(names(sort(tapply(mejor$puntaje, mejor$gen, max), decreasing = TRUE)), n)
  mejor <- mejor[mejor$gen %in% orden_gen, ]
  mejor$categoria <- factor(ifelse(mejor$puntaje >= 9, "alta", ifelse(mejor$puntaje >= 5, "media", "baja")),
                            levels = .niveles_categoria)
  # Respeta la categoria real (por ejemplo, portadores limitados a media).
  real <- d[match(paste(mejor$gen, mejor$muestra, mejor$puntaje), paste(d$gen, d$muestra, d$puntaje)), "categoria"]
  mejor$categoria <- factor(ifelse(is.na(real), as.character(mejor$categoria), real), levels = .niveles_categoria)
  mejor$gen <- factor(mejor$gen, levels = rev(orden_gen))
  mejor$text <- sprintf("%s en %s\nmejor puntaje %s (%s)", mejor$gen, mejor$muestra, mejor$puntaje, mejor$categoria)
  ggplot2::ggplot(mejor, ggplot2::aes(x = .data$muestra, y = .data$gen, fill = .data$categoria,
                                      text = .data$text)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$puntaje), size = 3.4, colour = "white", fontface = "bold") +
    ggplot2::scale_fill_manual(values = .colores_categoria, drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Prioridad", title = "Mejor puntaje por gen y muestra") +
    .tema_iei() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# ---- iei_gen ---------------------------------------------------------------

# Variantes de ClinVar con posicion, sumadas por posicion y tipo.
.variantes_proteina <- function(object) {
  cv <- object$clinvar
  if (is.null(cv)) cli::cli_abort("La ficha no tiene datos de ClinVar. Vuelve a crearla con {.code iei_gen(\"{object$gen}\", clinvar = TRUE)} y conexi\u00f3n.")
  v <- cv$variantes[!is.na(cv$variantes$pos_proteina), ]
  if (nrow(v) == 0) cli::cli_abort("Ninguna variante patogenica de ClinVar tiene posicion en la proteina.")
  v
}

# Dominio de cada posicion (el primero que la contiene) o NA.
.dominio_de <- function(pos, dom) {
  if (is.null(dom) || nrow(dom) == 0) return(rep(NA_character_, length(pos)))
  vapply(pos, function(p) {
    k <- which(dom$inicio <= p & dom$fin >= p)
    if (length(k)) dom$nombre[k[1]] else NA_character_
  }, character(1))
}

.colores_dominio_lollipop <- c("#8ec1ef", "#f6b38e", "#93d9bd", "#f5d27f", "#efaecb", "#b3a9ef", "#9fd4e3")

.grafico_lollipop <- function(object) {
  v <- .variantes_proteina(object)
  cv <- object$clinvar
  prot <- object$proteina
  dom <- if (!is.null(prot)) prot$dominios else NULL
  largo <- if (!is.null(prot)) prot$largo else max(v$pos_proteina)
  v$tipo <- factor(ifelse(v$tipo %in% names(.colores_tipo), v$tipo, "otra"), levels = names(.colores_tipo))
  conteo <- stats::aggregate(clinvar_id ~ pos_proteina + tipo, v, length)
  names(conteo)[3] <- "n"
  # Una cabeza por posicion: el tipo mas frecuente en ella.
  conteo <- conteo[order(conteo$pos_proteina, -conteo$n), ]
  total <- tapply(conteo$n, conteo$pos_proteina, sum)
  cab <- conteo[!duplicated(conteo$pos_proteina), ]
  cab$n <- as.integer(total[as.character(cab$pos_proteina)])
  ejemplo <- tapply(v$variante, v$pos_proteina, function(x) x[1])
  # De "NM_000061.3(BTK):c.334T>C (p.Tyr112His)" queda "Tyr112His".
  etiqueta <- ejemplo[as.character(cab$pos_proteina)]
  corta <- grepl("\\(p\\.[^)]+\\)", etiqueta)
  etiqueta[corta] <- sub("^.*\\(p\\.([^)]+)\\).*$", "\\1", etiqueta[corta])
  cab$etiqueta <- sub("^p\\.", "", etiqueta)
  cab$text <- sprintf("Posici\u00f3n %d\n%d variante(s), sobre todo %s\n%s", cab$pos_proteina, cab$n, cab$tipo,
                      ejemplo[as.character(cab$pos_proteina)])
  top <- utils::head(cab[order(-cab$n, cab$pos_proteina), ], 4)
  top <- top[top$n >= 2, ]

  v$dominio <- .dominio_de(v$pos_proteina, dom)
  en_dom <- mean(!is.na(v$dominio))
  subt <- sprintf("%d variantes patog\u00e9nicas en %d posiciones de %d amino\u00e1cidos", nrow(v),
                  length(unique(v$pos_proteina)), largo)
  if (!is.null(dom) && nrow(dom) > 0) {
    porc <- sort(table(v$dominio), decreasing = TRUE)
    if (length(porc)) subt <- paste0(subt, sprintf(".\nEl dominio %s concentra el %s.", .acortar(names(porc)[1], 60),
                                                  .pct(porc[[1]] / nrow(v))))
  }
  ymax <- max(cab$n)
  base <- 0
  g <- ggplot2::ggplot(cab, ggplot2::aes(text = .data$text)) +
    ggplot2::annotate("rect", xmin = 1, xmax = largo, ymin = base - 0.22 * ymax / 4 - 0.12,
                      ymax = base + 0.12, fill = "#e7e6e2")
  if (!is.null(dom) && nrow(dom) > 0) {
    dom$color <- rep(.colores_dominio_lollipop, length.out = nrow(dom))
    dom$ancho <- dom$fin - dom$inicio
    dom$etq <- ifelse(dom$ancho / largo * 80 >= nchar(dom$nombre), dom$nombre, "")
    alto_dom <- 0.28 + ymax * 0.06
    g <- g +
      ggplot2::geom_rect(data = dom, ggplot2::aes(xmin = .data$inicio, xmax = .data$fin,
                                                  ymin = base - alto_dom, ymax = base + alto_dom),
                         fill = dom$color, colour = "white", linewidth = 0.6, inherit.aes = FALSE) +
      ggplot2::geom_text(data = dom, ggplot2::aes(x = (.data$inicio + .data$fin) / 2, y = base, label = .data$etq),
                         colour = .marca, fontface = "bold", size = 3.3, inherit.aes = FALSE)
  }
  g <- g +
    ggplot2::geom_segment(data = cab, ggplot2::aes(x = .data$pos_proteina, xend = .data$pos_proteina,
                                                   y = base + 0.3 + ymax * 0.06, yend = .data$n + 0.3 + ymax * 0.06),
                          colour = "#b8b7b1", linewidth = 0.35) +
    ggplot2::geom_point(data = cab, ggplot2::aes(x = .data$pos_proteina, y = .data$n + 0.3 + ymax * 0.06,
                                                 fill = .data$tipo, shape = .data$tipo, size = .data$n),
                        colour = "white", stroke = 0.7, alpha = 0.95) +
    ggplot2::geom_rug(data = v, ggplot2::aes(x = .data$pos_proteina, colour = .data$tipo), sides = "b",
                      length = ggplot2::unit(0.035, "npc"), alpha = 0.7, show.legend = FALSE, inherit.aes = FALSE)
  if (nrow(top)) {
    g <- g + ggplot2::geom_text(data = top, ggplot2::aes(x = .data$pos_proteina, y = .data$n + 0.3 + ymax * 0.06,
                                                         label = .data$etiqueta),
                                vjust = -1.3, size = 3.1, fontface = "bold", colour = .tinta[["primaria"]],
                                check_overlap = TRUE, inherit.aes = FALSE)
  }
  g +
    ggplot2::scale_fill_manual(values = .colores_tipo, drop = TRUE, name = NULL) +
    ggplot2::scale_colour_manual(values = .colores_tipo, drop = TRUE, guide = "none") +
    ggplot2::scale_shape_manual(values = .formas_tipo, drop = TRUE, name = NULL) +
    ggplot2::scale_size_continuous(range = c(2.6, 7), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(1, largo), expand = ggplot2::expansion(mult = 0.015),
                                breaks = pretty(c(1, largo), 6)) +
    ggplot2::scale_y_continuous(breaks = function(l) {
      b <- pretty(c(1, max(1, l[2] - 0.3 - ymax * 0.06)))
      b[b >= 1 & b == round(b)] + 0.3 + ymax * 0.06
    }, labels = function(b) round(b - 0.3 - ymax * 0.06), expand = ggplot2::expansion(mult = c(0.08, 0.16))) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::labs(x = "Posici\u00f3n en la prote\u00edna (amino\u00e1cido)", y = "Variantes en la posici\u00f3n",
                  title = sprintf("%s: d\u00f3nde caen las variantes patog\u00e9nicas", object$gen),
                  subtitle = subt,
                  caption = sprintf("Patog\u00e9nicas y probablemente patog\u00e9nicas de ClinVar, consulta del %s. Dominios de Ensembl.", cv$fecha)) +
    .tema_iei() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

.grafico_dominios <- function(object) {
  v <- .variantes_proteina(object)
  prot <- object$proteina
  if (is.null(prot) || nrow(prot$dominios) == 0) cli::cli_abort("La ficha no tiene dominios de la prote\u00edna.")
  dom <- prot$dominios
  v$dominio <- .dominio_de(v$pos_proteina, dom)
  cubierto <- unique(unlist(lapply(seq_len(nrow(dom)), function(i) dom$inicio[i]:dom$fin[i])))
  d <- data.frame(dominio = c(dom$nombre, "Fuera de dominios"),
                  largo = c(dom$fin - dom$inicio + 1, prot$largo - length(cubierto)))
  d$variantes <- vapply(d$dominio, function(x) {
    if (x == "Fuera de dominios") sum(is.na(v$dominio)) else sum(v$dominio == x, na.rm = TRUE)
  }, numeric(1))
  d <- d[d$largo > 0, ]
  d$densidad <- 100 * d$variantes / d$largo
  esperado <- 100 * nrow(v) / prot$largo
  d$dominio <- factor(d$dominio, levels = d$dominio[order(d$densidad)])
  d$fuera <- d$dominio == "Fuera de dominios"
  d$text <- sprintf("%s\n%d variantes en %d aa\n%.1f por cada 100 aa", d$dominio, d$variantes, d$largo, d$densidad)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$densidad, y = .data$dominio, text = .data$text)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$fuera), width = 0.62, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = esperado, linetype = "22", colour = .tinta[["secundaria"]], linewidth = 0.5) +
    ggplot2::annotate("text", x = esperado, y = Inf, label = sprintf("promedio de la prote\u00edna (%.1f)", esperado),
                      hjust = -0.05, vjust = 1.4, size = 3.1, colour = .tinta[["secundaria"]]) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d en %d aa", .data$variantes, .data$largo)),
                       hjust = -0.12, size = 3.2, colour = .tinta[["secundaria"]]) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = .marca_2, `TRUE` = "#b8b7b1")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.28))) +
    ggplot2::labs(x = "Variantes patog\u00e9nicas por cada 100 amino\u00e1cidos", y = NULL,
                  title = sprintf("%s: qu\u00e9 dominios acumulan m\u00e1s variantes", object$gen),
                  subtitle = "Una barra m\u00e1s larga que la l\u00ednea indica un dominio con m\u00e1s variantes que el promedio.",
                  caption = sprintf("ClinVar, consulta del %s. Dominios de Ensembl.", object$clinvar$fecha)) +
    .tema_iei() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

.grafico_clinvar <- function(object) {
  cv <- object$clinvar
  if (is.null(cv)) cli::cli_abort("La ficha no tiene datos de ClinVar. Vuelve a crearla con {.code iei_gen(\"{object$gen}\", clinvar = TRUE)} y conexi\u00f3n.")
  k <- cv$conteo
  k$color <- unname(.divergente[match(k$clasificacion, names(.divergente))])
  k$color[is.na(k$color)] <- "#c9c8c2"
  k$nombre <- c(Patogenica = "Patog\u00e9nica", `Probablemente patogenica` = "Probablemente patog\u00e9nica",
                `Significado incierto` = "Significado incierto", `Probablemente benigna` = "Probablemente benigna",
                Benigna = "Benigna")[k$clasificacion]
  k$nombre[is.na(k$nombre)] <- k$clasificacion[is.na(k$nombre)]
  k$nombre <- factor(k$nombre, levels = rev(k$nombre))
  tot <- sum(k$variantes)
  k$prop <- k$variantes / max(tot, 1)
  k$text <- sprintf("%s\n%d variantes (%s)", k$nombre, k$variantes, .pct(k$prop))
  pat <- sum(k$variantes[k$clasificacion %in% c("Patogenica", "Probablemente patogenica")])
  ggplot2::ggplot(k, ggplot2::aes(x = .data$variantes, y = .data$nombre, text = .data$text)) +
    ggplot2::geom_col(fill = k$color, width = 0.66) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d  \u00b7  %s", .data$variantes, .pct(.data$prop))),
                       hjust = -0.1, size = 3.3, colour = .tinta[["secundaria"]]) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) +
    ggplot2::labs(x = "Variantes en ClinVar", y = NULL,
                  title = sprintf("%s en ClinVar", object$gen),
                  subtitle = sprintf("%d de %d variantes (%s) son patog\u00e9nicas o probablemente patog\u00e9nicas.",
                                     pat, tot, .pct(pat / max(tot, 1))),
                  caption = sprintf("Clasificaci\u00f3n germinal de ClinVar, consulta del %s.", cv$fecha)) +
    .tema_iei() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 11, colour = .tinta[["primaria"]]))
}

#' @rdname autoplot.iei
#' @exportS3Method ggplot2::autoplot
autoplot.iei_gen <- function(object, tipo = c("proteina", "dominios", "clinvar"), ...) {
  .necesita_ggplot()
  switch(match.arg(tipo),
         proteina = .grafico_lollipop(object),
         dominios = .grafico_dominios(object),
         clinvar = .grafico_clinvar(object))
}

#' @export
plot.iei_gen <- function(x, tipo = c("proteina", "dominios", "clinvar"), ...) autoplot.iei_gen(x, tipo = tipo, ...)

# ---- iei_busqueda ----------------------------------------------------------

#' @rdname autoplot.iei
#' @exportS3Method ggplot2::autoplot
autoplot.iei_busqueda <- function(object, tipo = c("ranking", "mapa", "cariotipo"), n = 20, ...) {
  .necesita_ggplot()
  tipo <- match.arg(tipo)
  if (tipo == "cariotipo") return(iei_cariograma(object, etiquetas = n))
  d <- utils::head(tibble::as_tibble(unclass(object)), n)
  d$gen <- factor(d$gen, levels = rev(d$gen))
  if (tipo == "ranking") {
    d$herencia <- ifelse(is.na(d$herencia), "sin dato", d$herencia)
    d$enfermedad_txt <- if ("enfermedad" %in% names(d)) .acortar(sub("\\s*\\((OMIM|ORPHA)[^)]*\\)\\s*$", "", ifelse(is.na(d$enfermedad), "", d$enfermedad)), 46) else ""
    d$text <- sprintf("%s\nsimilitud %.2f, especificidad %.2f\n%d t\u00e9rminos exactos\n%s", d$gen, d$similitud,
                      d$especificidad, d$exactos, d$enfermedad_txt)
    fen <- attr(object, "fenotipo")
    subt <- if (length(fen)) paste0("Caso: ", .acortar(paste(stats::na.omit(.hpo_nombre(fen)), collapse = ", "), 110)) else NULL
    return(
      ggplot2::ggplot(d, ggplot2::aes(y = .data$gen, text = .data$text)) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$similitud, yend = .data$gen, colour = .data$herencia),
                              linewidth = 1.6, alpha = 0.35, lineend = "round") +
        ggplot2::geom_point(ggplot2::aes(x = .data$similitud, fill = .data$herencia, size = .data$especificidad),
                            shape = 21, colour = "white", stroke = 0.9) +
        ggplot2::geom_text(ggplot2::aes(x = .data$similitud, label = .data$enfermedad_txt), hjust = -0.08,
                           nudge_x = 0.025, size = 2.9, colour = .tinta[["secundaria"]]) +
        ggplot2::scale_fill_manual(values = .colores_herencia, na.value = "#a3a29b", name = "Herencia") +
        ggplot2::scale_colour_manual(values = .colores_herencia, na.value = "#a3a29b", guide = "none") +
        ggplot2::scale_size_continuous(range = c(2.5, 7), name = "Especificidad", limits = c(0, 1),
                                       breaks = c(0.25, 0.5, 0.75)) +
        ggplot2::guides(size = ggplot2::guide_legend(override.aes = list(fill = "#a3a29b")),
                        fill = ggplot2::guide_legend(override.aes = list(size = 4.5))) +
        ggplot2::scale_x_continuous(limits = c(0, 1.75), breaks = seq(0, 1, 0.25),
                                    expand = ggplot2::expansion(mult = c(0, 0))) +
        ggplot2::labs(x = "Similitud con la enfermedad m\u00e1s parecida del gen", y = NULL,
                      title = "Genes candidatos por fenotipo", subtitle = subt) +
        .tema_iei() +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(face = "bold.italic", colour = .tinta[["primaria"]]))
    )
  }
  m <- attr(object, "matriz")[, as.character(rev(d$gen)), drop = FALSE]
  largo <- data.frame(termino = rep(rownames(m), ncol(m)), gen = rep(colnames(m), each = nrow(m)),
                      similitud = as.vector(m))
  largo$termino <- factor(.acortar(.hpo_nombre(largo$termino), 34),
                          levels = unique(.acortar(.hpo_nombre(rownames(m)), 34)))
  largo$gen <- factor(largo$gen, levels = levels(d$gen))
  largo$text <- sprintf("%s\n%s: %.2f", largo$gen, largo$termino, largo$similitud)
  ggplot2::ggplot(largo, ggplot2::aes(x = .data$termino, y = .data$gen, fill = .data$similitud, text = .data$text)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$similitud), colour = .data$similitud >= 0.55),
                       size = 2.9, show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(colours = c("#f7f7f5", .rampa_azul), limits = c(0, 1), name = "Similitud") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = .tinta[["primaria"]])) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL, title = "Similitud de cada signo del caso con cada gen") +
    .tema_iei(leyenda = "right") +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x.top = ggplot2::element_text(angle = 35, hjust = 0, vjust = 0),
                   axis.text.y = ggplot2::element_text(face = "bold.italic", colour = .tinta[["primaria"]]))
}

#' @export
plot.iei_busqueda <- function(x, tipo = c("ranking", "mapa", "acordes", "circos", "cariotipo"), n = 20, ...) {
  tipo <- match.arg(tipo)
  if (tipo == "acordes") return(.grafico_acordes(x, n = min(n, 15)))
  if (tipo == "circos") return(iei_circos(x, etiquetas = n))
  autoplot.iei_busqueda(x, tipo = tipo, n = n)
}

# ---- iei_validacion --------------------------------------------------------

#' @rdname autoplot.iei
#' @exportS3Method ggplot2::autoplot
autoplot.iei_validacion <- function(object, hasta = 20, ...) {
  .necesita_ggplot()
  y <- tibble::as_tibble(unclass(object)[c("metodo", "puesto")])
  k <- seq_len(hasta)
  curva <- do.call(rbind, lapply(split(y$puesto, y$metodo), function(p) {
    data.frame(puesto = k, porcentaje = vapply(k, function(i) 100 * mean(p <= i), numeric(1)))
  }))
  curva$metodo <- rep(names(split(y$puesto, y$metodo)), each = hasta)
  curva$text <- sprintf("%s: %.0f%% de los casos entre los %d primeros", curva$metodo, curva$porcentaje, curva$puesto)
  ggplot2::ggplot(curva, ggplot2::aes(x = .data$puesto, y = .data$porcentaje, colour = .data$metodo,
                                      text = .data$text, group = .data$metodo)) +
    ggplot2::geom_step(linewidth = 1.2) +
    ggplot2::scale_colour_manual(values = c(enfermedad = .categoricos[1], semantica = .categoricos[2], exacta = "#a3a29b"),
                                 breaks = c("enfermedad", "semantica", "exacta"),
                                 labels = c(enfermedad = "Por enfermedad", semantica = "Por gen", exacta = "Exacta")) +
    ggplot2::scale_x_continuous(breaks = unique(c(1, seq(5, max(hasta, 5), by = 5))), expand = ggplot2::expansion(add = 0.3)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 25), labels = function(v) paste0(v, " %")) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::labs(x = "Puesto del gen correcto", y = "Casos con el gen en ese puesto o antes", colour = NULL,
                  title = "Validaci\u00f3n con casos simulados",
                  subtitle = sprintf("%d casos", attr(object, "parametros")$n)) +
    .tema_iei()
}

#' @export
plot.iei_validacion <- function(x, hasta = 20, ...) autoplot.iei_validacion(x, hasta = hasta)

# ---- iei_ab1 ---------------------------------------------------------------

#' @rdname autoplot.iei
#' @exportS3Method ggplot2::autoplot
autoplot.iei_ab1 <- function(object, desde = 1, hasta = min(desde + 59, nchar(object$bases)), ...) {
  .necesita_ggplot()
  n <- nchar(object$bases)
  desde <- max(1, desde); hasta <- min(n, hasta)
  picos <- object$picos[desde:hasta]
  rango <- max(1, min(picos) - 6):min(nrow(object$trazas), max(picos) + 6)
  tr <- object$trazas[rango, ]
  largo <- data.frame(x = rep(rango, 4), senal = unlist(tr, use.names = FALSE),
                      base = rep(c("A", "C", "G", "T"), each = length(rango)))
  et <- data.frame(x = picos, base = strsplit(substr(object$bases, desde, hasta), "")[[1]],
                   calidad = object$calidad[desde:hasta])
  alto <- max(largo$senal)
  ggplot2::ggplot(largo, ggplot2::aes(x = .data$x, y = .data$senal, colour = .data$base)) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_text(data = et, ggplot2::aes(x = .data$x, y = alto * 1.07, label = .data$base),
                       inherit.aes = FALSE, size = 3) +
    ggplot2::scale_colour_manual(values = c(A = "#1b9e3e", C = "#2166ac", G = "#222222", T = "#d7301f")) +
    ggplot2::labs(x = NULL, y = "Se\u00f1al", colour = NULL,
                  title = sprintf("%s, bases %d a %d", object$archivo, desde, hasta)) +
    .tema_iei() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
}

#' @export
plot.iei_ab1 <- function(x, ...) autoplot.iei_ab1(x, ...)
