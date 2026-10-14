# Busqueda de genes por fenotipo y de terminos HPO por texto.

#' Buscar terminos HPO por texto
#'
#' Encuentra terminos de la Human Phenotype Ontology a partir de una palabra
#' o frase, en espanol (para los signos mas frecuentes en inmunodeficiencias)
#' o en ingles (todos los terminos y sus sinonimos). Sirve para armar el
#' `fenotipo` de [iei_priorizar()] o [iei_buscar_fenotipo()].
#'
#' @param texto Palabra o frase, por ejemplo `"neumonia"` o `"candidiasis"`.
#' @param n Numero maximo de terminos a devolver.
#' @return Un tibble con `hpo_id`, `termino`, `genes_panel` (cuantos genes del
#'   panel lo tienen anotado) e `ic` (contenido de informacion: mas alto es mas
#'   especifico).
#' @export
#' @examples
#' iei_buscar_hpo("neumonia")
#' iei_buscar_hpo("candidiasis")
iei_buscar_hpo <- function(texto, n = 15) {
  if (length(texto) != 1 || !nzchar(trimws(texto))) cli::cli_abort("{.arg texto} debe ser una palabra o frase.")
  o <- .ont_textos()
  t <- ieiprio::hpo_terminos
  q <- .sin_tildes(trimws(texto))
  # Del diccionario en espanol: entradas que contienen el texto.
  dic <- .hpo_espanol[grepl(q, names(.hpo_espanol), fixed = TRUE)]
  palabras <- strsplit(q, "\\s+")[[1]]
  contiene <- function(v) Reduce(`&`, lapply(palabras, function(p) grepl(p, v, fixed = TRUE)))
  en_nombre <- which(contiene(o$nombres_norm))
  en_sinonimo <- which(vapply(o$sinon_norm, function(s) length(s) > 0 && any(contiene(s)), logical(1)))
  ids <- unique(c(unname(dic), t$id[en_nombre], t$id[en_sinonimo]))
  if (length(ids) == 0) {
    cli::cli_inform("No encontr\u00e9 t\u00e9rminos HPO con {.val {texto}}. Prueba en ingl\u00e9s o con una palabra m\u00e1s corta.")
    return(tibble::tibble(hpo_id = character(0), termino = character(0), genes_panel = integer(0), ic = numeric(0)))
  }
  cuenta <- table(ieiprio::hpo_genes$hpo_id)
  out <- tibble::tibble(
    hpo_id = ids,
    termino = .hpo_nombre(ids),
    genes_panel = as.integer(ifelse(is.na(cuenta[ids]), 0, cuenta[ids])),
    ic = round(unname(o$ic[o$idx[ids]]), 2)
  )
  # Primero lo del diccionario, luego los mas usados en el panel y los nombres cortos.
  orden <- order(!out$hpo_id %in% dic, -out$genes_panel, nchar(out$termino))
  utils::head(out[orden, ], n)
}

#' Buscar genes candidatos por fenotipo
#'
#' Ordena los genes del panel segun cuanto se parece el fenotipo del caso al
#' fenotipo conocido de cada gen. Sirve antes de tener variantes, para saber
#' que genes mirar primero, o para elegir el panel de una secuenciacion
#' dirigida.
#'
#' Con `similitud = "enfermedad"` (por defecto) el caso se compara con cada
#' enfermedad anotada en HPO por separado, porque un mismo gen puede causar
#' cuadros muy distintos (STAT1 con ganancia o con perdida de funcion, por
#' ejemplo). Para cada signo del caso se busca el signo mas parecido de la
#' enfermedad usando la jerarquia de HPO, y ese parecido se pesa por la
#' frecuencia del signo en la enfermedad: un signo constante cuenta entero y
#' uno que aparece en el 5% de los pacientes cuenta la mitad. Si la
#' enfermedad excluye un signo que el caso tiene (anotacion NOT de HPO), ese
#' signo resta. Los `ausentes` restan segun lo frecuentes que son en la
#' enfermedad. El valor del gen es el de su enfermedad mas parecida, que se
#' muestra en la columna `enfermedad`. Los genes sin anotaciones por
#' enfermedad se comparan con sus terminos de gen.
#'
#' Con `"semantica"` se compara con todos los terminos del gen juntos, como
#' en la version 0.3, y con `"exacta"` solo cuentan los terminos identicos.
#' [iei_validar_fenotipo()] compara los tres metodos con casos simulados.
#'
#' @param fenotipo Codigos HPO o texto en espanol o ingles, como en
#'   [iei_priorizar()].
#' @param n Numero de genes a devolver.
#' @param evidencia Niveles de evidencia de PanelApp a considerar.
#' @param herencia Opcional, limitar a uno o varios modos de herencia (`"AR"`,
#'   `"AD"`, `"XL"`, `"MT"`).
#' @param similitud `"enfermedad"`, `"semantica"` o `"exacta"`.
#' @param ausentes Signos que el caso no tiene y que ayudan a descartar, por
#'   ejemplo `"eccema"` cuando se piensa en Wiskott-Aldrich. Mismo formato que
#'   `fenotipo`.
#' @param p_valor Si `TRUE`, agrega un valor p empirico: la proporcion de
#'   fenotipos al azar del mismo tamano que dan una similitud igual o mayor
#'   con ese gen. Ayuda a ver si un gen sale primero por el caso o porque su
#'   fenotipo es tan amplio que encaja con casi todo.
#' @param permutaciones Fenotipos al azar para el valor p.
#' @return Un tibble de clase `iei_busqueda` con `gen`, `similitud` (0 a 1),
#'   `especificidad` (0 a 1, cuanto del fenotipo conocido del gen se parece al
#'   del caso; desempata), `exactos` (terminos identicos), `terminos_gen`,
#'   `enfermedad` (la mas parecida, con su identificador OMIM u Orphanet),
#'   `herencia`, `categorias_iuis`, `p_valor` si se pidio y `coincidencias`
#'   (para cada termino del caso, el termino mas parecido de la enfermedad o
#'   del gen). Se grafica con [plot()] o `autoplot()`.
#' @export
#' @examples
#' iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente"))
#' iei_buscar_fenotipo("HP:0002728", herencia = "AD")
#' iei_buscar_fenotipo(c("trombocitopenia", "infecciones recurrentes"), ausentes = "eccema")
iei_buscar_fenotipo <- function(fenotipo, n = 20, evidencia = "verde", herencia = NULL,
                                similitud = c("enfermedad", "semantica", "exacta"),
                                ausentes = NULL, p_valor = FALSE, permutaciones = 200) {
  similitud <- match.arg(similitud)
  caso <- .fenotipo_codigos(fenotipo)
  if (is.null(caso)) cli::cli_abort("Indica al menos un t\u00e9rmino de {.arg fenotipo}.")
  ausentes <- .fenotipo_codigos(ausentes)
  panel <- iei_panel(herencia = herencia, evidencia = evidencia)
  hpo <- ieiprio::hpo_genes[ieiprio::hpo_genes$gen %in% panel$gen, ]
  o <- .ontologia()
  hpo <- hpo[o$fenotipico[o$idx[hpo$hpo_id]] %in% TRUE, ]
  genes <- unique(hpo$gen)
  r <- .fenotipo_por_gen(caso, genes, similitud, ausentes, hpo)
  m <- attr(r, "matriz")
  por_gen <- split(hpo$hpo_id, hpo$gen)
  exactos <- vapply(genes, function(g) sum(caso %in% por_gen[[g]]), integer(1))
  info <- panel[match(genes, panel$gen), ]
  # Especificidad: en promedio, cuanto se parece cada termino del gen al caso.
  # Desempata a favor de genes cuyo fenotipo gira en torno al del caso.
  ci <- match(caso, o$ids)
  por_gen_i <- lapply(por_gen, match, o$ids)
  especificidad <- vapply(genes, function(g) mean(.similitud_idx(por_gen_i[[g]], ci)), numeric(1))
  out <- tibble::tibble(
    gen = genes,
    similitud = r$similitud,
    especificidad = round(unname(especificidad), 3),
    exactos = unname(exactos),
    terminos_gen = unname(lengths(por_gen[genes])),
    enfermedad = ifelse(is.na(r$enfermedad), NA_character_,
                        sprintf("%s (%s)", trimws(r$enfermedad), r$enfermedad_id)),
    herencia = info$herencia,
    categorias_iuis = info$categorias_iuis
  )
  out <- out[order(-out$similitud, -out$especificidad, out$gen), ]
  out <- utils::head(out, n)
  if (isTRUE(p_valor)) {
    out$p_valor <- signif(.p_valor_fenotipo(out$similitud, caso, out$gen, similitud,
                                            permutaciones, hpo, ieiprio::hpo_enfermedades), 2)
  }
  enf <- ieiprio::hpo_enfermedades
  id_enf <- r$enfermedad_id[match(out$gen, r$gen)]
  out$coincidencias <- unname(vapply(seq_len(nrow(out)), function(i) {
    g <- out$gen[i]
    if (!is.na(id_enf[i])) {
      e <- enf[enf$gen == g & enf$enfermedad_id == id_enf[i], ]
      .mejor_coincidencia(caso, e$hpo_id[!e$excluido], e$hpo_id[e$excluido])
    } else {
      .mejor_coincidencia(caso, por_gen[[g]])
    }
  }, ""))
  attr(out, "fenotipo") <- caso
  attr(out, "ausentes") <- ausentes
  attr(out, "similitud") <- similitud
  attr(out, "matriz") <- m[, out$gen, drop = FALSE]
  class(out) <- c("iei_busqueda", class(out))
  out
}

# Para cada termino del caso, el termino del gen (o de la enfermedad) con el
# que mas comparte. Marca los que la enfermedad excluye.
.mejor_coincidencia <- function(caso, terminos_gen, excluidos = character(0)) {
  o <- .ontologia()
  partes <- vapply(caso, function(p) {
    k <- o$idx[[p]]
    if (length(excluidos) > 0 && any(o$idx[excluidos] %in% o$anc[[k]]))
      return(sprintf("%s (la enfermedad lo excluye)", .hpo_nombre(p)))
    if (length(terminos_gen) == 0) return(sprintf("%s (sin parecido)", .hpo_nombre(p)))
    if (p %in% terminos_gen) return(sprintf("%s (igual)", .hpo_nombre(p)))
    mejor <- vapply(terminos_gen, function(g) .similitud_terminos(p, g), numeric(1))
    if (max(mejor) <= 0) return(sprintf("%s (sin parecido)", .hpo_nombre(p)))
    k <- which.max(mejor)
    sprintf("%s ~ %s (%.2f)", .hpo_nombre(p), .hpo_nombre(terminos_gen[k]), mejor[k])
  }, "")
  paste(partes, collapse = "; ")
}

#' @export
print.iei_busqueda <- function(x, ...) {
  caso <- attr(x, "fenotipo")
  cli::cli_text("{.strong Genes candidatos por fenotipo} ({attr(x, 'similitud')})")
  cli::cli_text("Fenotipo: {paste(sprintf('%s %s', caso, .hpo_nombre(caso)), collapse = '; ')}")
  aus <- attr(x, "ausentes")
  if (length(aus) > 0) cli::cli_text("Ausentes: {paste(sprintf('%s %s', aus, .hpo_nombre(aus)), collapse = '; ')}")
  cols <- intersect(c("gen", "similitud", "especificidad", "p_valor", "enfermedad", "herencia",
                      "categorias_iuis"), names(x))
  y <- tibble::as_tibble(unclass(x)[cols])
  print(y, n = nrow(y))
  cli::cli_text("{.emph Detalle de cada coincidencia en la columna coincidencias. Grafico: plot(x) o autoplot(x).}")
  invisible(x)
}
