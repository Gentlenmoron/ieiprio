# Validacion de la busqueda por fenotipo con pacientes simulados.

#' Validar la busqueda por fenotipo con pacientes simulados
#'
#' Mide que tan bien [iei_buscar_fenotipo()] encuentra el gen correcto. Para
#' cada caso simulado elige al azar una enfermedad del panel, toma algunos de
#' sus signos (los frecuentes con mas probabilidad), cambia una parte por un
#' termino mas general, como haria un medico que describe "infecciones
#' recurrentes" en lugar de "neumonia por Pneumocystis", y agrega signos que
#' no tienen que ver. Luego anota en que puesto queda el gen verdadero con
#' cada metodo.
#'
#' Es una validacion interna: los casos salen de las mismas anotaciones que
#' usa el puntaje, asi que el resultado es optimista y sirve para comparar
#' metodos y elegir parametros, no como estimacion del rendimiento con
#' pacientes reales. Para eso hace falta una cohorte con diagnostico conocido.
#'
#' @param n Numero de casos simulados.
#' @param terminos Signos tomados de la enfermedad en cada caso.
#' @param imprecision Proporcion de signos que se reemplazan por un termino
#'   padre.
#' @param ruido Signos al azar agregados a cada caso.
#' @param metodos Metodos a comparar: `"enfermedad"`, `"semantica"` y
#'   `"exacta"`.
#' @param semilla Semilla para que la simulacion se pueda repetir.
#' @return Un tibble de clase `iei_validacion` con una fila por caso y metodo:
#'   `caso`, `gen`, `enfermedad`, `metodo` y `puesto` del gen verdadero (1 es
#'   el primero). [print()] resume el porcentaje de casos en el primer
#'   puesto, entre los 5 y entre los 10 primeros; `autoplot()` grafica la
#'   curva acumulada.
#' @export
#' @examples
#' v <- iei_validar_fenotipo(n = 20)
#' v
iei_validar_fenotipo <- function(n = 100, terminos = 5, imprecision = 0.3, ruido = 1,
                                 metodos = c("enfermedad", "semantica", "exacta"), semilla = 1) {
  metodos <- match.arg(metodos, c("enfermedad", "semantica", "exacta"), several.ok = TRUE)
  # Semilla propia sin alterar la del usuario.
  previa <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(previa)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", previa, envir = globalenv()),
          add = TRUE)
  set.seed(semilla)
  o <- .ontologia()
  enf <- ieiprio::hpo_enfermedades
  enf <- enf[!enf$excluido & !enf$frecuencia %in% 0 & o$fenotipico[o$idx[enf$hpo_id]] %in% TRUE, ]
  enf <- enf[enf$gen %in% ieiprio::iuis_panel$gen, ]
  clave <- paste(enf$gen, enf$enfermedad_id, sep = "|")
  validas <- names(which(table(clave) >= terminos))
  if (length(validas) == 0) cli::cli_abort("Ninguna enfermedad tiene {terminos} signos anotados.")
  genes <- unique(ieiprio::iuis_panel$gen)
  pool <- unique(ieiprio::hpo_genes$hpo_id[o$fenotipico[o$idx[ieiprio::hpo_genes$hpo_id]] %in% TRUE])
  padres <- stats::setNames(ieiprio::hpo_terminos$padres, ieiprio::hpo_terminos$id)
  elegidas <- sample(validas, n, replace = length(validas) < n)
  filas <- vector("list", n)
  cli::cli_progress_bar("Simulando casos", total = n)
  for (i in seq_len(n)) {
    e <- enf[clave == elegidas[i], ]
    prob <- ifelse(is.na(e$frecuencia), 0.5, e$frecuencia)
    caso <- sample(e$hpo_id, terminos, prob = prob)
    cambia <- stats::runif(terminos) < imprecision
    caso[cambia] <- vapply(caso[cambia], function(t) {
      pa <- padres[[t]]
      pa <- pa[o$fenotipico[o$idx[pa]] %in% TRUE & pa != "HP:0000118"]
      if (length(pa) == 0) t else pa[1]
    }, "")
    if (ruido > 0) caso <- c(caso, sample(setdiff(pool, caso), ruido))
    caso <- unique(caso)
    filas[[i]] <- do.call(rbind, lapply(metodos, function(m) {
      r <- .fenotipo_por_gen(caso, genes, m)
      s <- r$similitud[match(e$gen[1], r$gen)]
      # Puesto conservador: los empates con el gen verdadero cuentan como
      # mejores que el.
      tibble::tibble(caso = i, gen = e$gen[1], enfermedad = e$enfermedad[1], metodo = m,
                     puesto = sum(r$similitud >= s))
    }))
    cli::cli_progress_update()
  }
  out <- do.call(rbind, filas)
  attr(out, "parametros") <- list(n = n, terminos = terminos, imprecision = imprecision,
                                  ruido = ruido, genes = length(genes))
  class(out) <- c("iei_validacion", class(out))
  out
}

.resumen_validacion <- function(x) {
  y <- tibble::as_tibble(unclass(x)[c("metodo", "puesto")])
  r <- lapply(split(y$puesto, factor(y$metodo, levels = unique(y$metodo))), function(p) {
    c(primero = mean(p <= 1), top5 = mean(p <= 5), top10 = mean(p <= 10), mediana = stats::median(p))
  })
  tibble::tibble(metodo = names(r),
                 primero = round(100 * vapply(r, `[[`, 0, "primero"), 1),
                 top5 = round(100 * vapply(r, `[[`, 0, "top5"), 1),
                 top10 = round(100 * vapply(r, `[[`, 0, "top10"), 1),
                 mediana = vapply(r, `[[`, 0, "mediana"))
}

#' @export
print.iei_validacion <- function(x, ...) {
  p <- attr(x, "parametros")
  cli::cli_text("{.strong Validacion de la busqueda por fenotipo}")
  cli::cli_text("{p$n} casos simulados, {p$terminos} signos por caso, {round(100 * p$imprecision)}% imprecisos, {p$ruido} al azar; {p$genes} genes en competencia.")
  print(.resumen_validacion(x))
  cli::cli_text("{.emph Porcentaje de casos con el gen correcto en el primer puesto, entre los 5 y los 10 primeros. Validacion interna, optimista.}")
  invisible(x)
}
