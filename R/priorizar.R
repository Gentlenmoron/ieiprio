#' Pesos por defecto de la priorizacion
#'
#' Devuelve la lista de puntos que usa [iei_priorizar()]. Puedes modificar
#' cualquier valor y pasar la lista en el argumento `pesos`.
#'
#' @return Una lista con nombres.
#' @export
#' @examples
#' p <- iei_pesos()
#' p$clinvar_patogenica <- 5
iei_pesos <- function() {
  list(
    clinvar_patogenica = 4,
    clinvar_conflicto  = 1,
    impacto_alto       = 3,
    impacto_moderado   = 1,
    ausente_gnomad     = 1,
    herencia_compatible = 2,
    portador_ar        = -2,
    fenotipo_maximo    = 3,
    corte_alta         = 9,
    corte_media        = 5
  )
}

#' Priorizar variantes
#'
#' Asigna a cada variante un puntaje y una categoria (alta, media o baja), con
#' una columna `razones` que explica de donde sale cada punto. Las variantes
#' clasificadas como benignas en ClinVar se descartan y quedan en la
#' bitacora.
#'
#' Los componentes del puntaje son estos.
#' * **ClinVar.** Patogenica o probablemente patogenica suma 4. Con
#'   interpretaciones en conflicto suma 1.
#' * **Impacto.** Alto (stop, frameshift, splicing canonico) suma 3.
#'   Moderado (missense, inframe) suma 1. En genes con enfermedad por
#'   ganancia de funcion una missense cuenta como impacto alto.
#' * **Frecuencia.** Ausente de gnomAD suma 1.
#' * **Herencia.** Genotipo compatible con el modo del gen suma 2. En genes
#'   autosomicos recesivos, dos heterocigotos en el mismo gen y la misma
#'   muestra se marcan como posible heterocigoto compuesto y suman 2; un
#'   heterocigoto solo resta 2, se marca como portador y nunca queda en
#'   prioridad alta, aunque el fenotipo encaje.
#' * **Fenotipo.** Proporcion de los terminos HPO del paciente presentes en
#'   el gen, escalada de 0 a 3. O 3 puntos si el gen esta en la `categoria`
#'   IUIS indicada.
#'
#' @param x Tabla devuelta por [iei_anotar()] o [iei_filtrar_frecuencia()].
#' @param fenotipo Terminos HPO del paciente, por ejemplo
#'   `c("HP:0004313", "HP:0002719")`.
#' @param categoria Categorias IUIS sospechadas (1 a 10), alternativa a
#'   `fenotipo`.
#' @param pesos Lista de pesos, ver [iei_pesos()].
#' @param panel,hpo,categorias Tablas de referencia. Por defecto los datos
#'   del paquete.
#'
#' @return La tabla con el desglose del puntaje (`p_clinvar`, `p_impacto`,
#'   `p_frecuencia`, `p_herencia`, `p_fenotipo`), el `puntaje` total, la
#'   `categoria`, las `razones`, una `nota` clinica y `acmg`, con los
#'   criterios ACMG/AMP sugeridos (PVS1, PM2_Supporting, PM3, PP4, BA1).
#'   Los criterios son orientativos y deben ser revisados por un
#'   especialista. Ordenada por muestra y puntaje.
#' @export
#' @examples
#' \dontrun{
#' prior <- iei_priorizar(anot, fenotipo = c("HP:0004313", "HP:0002719"))
#' }
iei_priorizar <- function(x, fenotipo = NULL, categoria = NULL, pesos = iei_pesos(),
                          panel = ieiprio::iuis_panel, hpo = ieiprio::hpo_genes,
                          categorias = ieiprio::iuis_categorias) {
  necesarias <- c("muestra", "gen", "genotipo", "impacto", "consecuencia",
                  "af_gnomad", "clinvar")
  faltan <- setdiff(necesarias, names(x))
  if (length(faltan) > 0) {
    cli::cli_abort(c("A {.arg x} le faltan columnas: {.val {faltan}}.",
                     "i" = "Pasala antes por {.fn iei_filtrar} y {.fn iei_anotar}."))
  }
  if (!is.null(fenotipo) && !all(grepl("^HP:[0-9]{7}$", fenotipo))) {
    cli::cli_abort("{.arg fenotipo} debe tener terminos HPO como {.val HP:0004313}.")
  }
  pesos <- utils::modifyList(iei_pesos(), pesos)
  a <- .attrs_ieiprio(x)

  # ClinVar: descartar benignas
  cv <- .clasificar_clinvar(x$clinvar)
  benigna <- !is.na(cv) & cv == "benigna"
  x <- x[!benigna, ]
  cv <- cv[!benigna]
  if (!is.null(a$ieiprio_bitacora)) {
    a$ieiprio_bitacora <- rbind(a$ieiprio_bitacora,
                                tibble::tibble(paso = "ClinVar benigna", filas = nrow(x),
                                               descartadas = sum(benigna)))
  }
  a$ieiprio_fenotipo <- fenotipo
  a$ieiprio_categoria <- categoria
  a$ieiprio_pesos <- pesos

  if (nrow(x) == 0) {
    for (v in c("p_clinvar", "p_impacto", "p_frecuencia", "p_herencia", "p_fenotipo", "puntaje"))
      x[[v]] <- numeric(0)
    for (v in c("categoria", "razones", "nota", "acmg")) x[[v]] <- character(0)
    return(.poner_attrs(x, a))
  }

  info <- panel[match(x$gen, panel$gen), ]
  herencia <- info$herencia
  gof <- !is.na(info$ganancia_funcion) & info$ganancia_funcion

  partes <- list()
  nota <- rep("", nrow(x))
  agregar <- function(cond, puntos, texto) {
    list(puntos = ifelse(cond, puntos, 0), texto = ifelse(cond, sprintf("%s (%+g)", texto, puntos), ""))
  }

  partes$clinvar_p <- agregar(cv %in% "patogenica", pesos$clinvar_patogenica, "ClinVar patogenica")
  partes$clinvar_c <- agregar(cv %in% "conflicto", pesos$clinvar_conflicto, "ClinVar en conflicto")

  missense <- grepl("missense", x$consecuencia %||% "", fixed = TRUE)
  alto <- (!is.na(x$impacto) & x$impacto == "HIGH") | (gof & missense)
  moderado <- !alto & !is.na(x$impacto) & x$impacto == "MODERATE"
  partes$alto <- agregar(alto, pesos$impacto_alto,
                         ifelse(gof & missense & x$impacto %in% "MODERATE",
                                "missense en gen de ganancia de funcion", "impacto alto"))
  partes$moderado <- agregar(moderado, pesos$impacto_moderado, "impacto moderado")
  partes$gnomad <- agregar(is.na(x$af_gnomad), pesos$ausente_gnomad, "ausente en gnomAD")

  # Herencia
  h <- .evaluar_herencia(x, herencia)
  partes$herencia <- agregar(h$compatible, pesos$herencia_compatible, h$texto)
  partes$portador <- agregar(h$portador_ar, pesos$portador_ar, "heterocigoto unico en gen AR")
  nota <- h$nota

  # Fenotipo
  f <- .puntaje_fenotipo(x$gen, fenotipo, categoria, hpo, categorias, pesos$fenotipo_maximo)
  partes$fenotipo <- list(puntos = f$puntos, texto = f$texto)

  puntaje <- Reduce(`+`, lapply(partes, `[[`, "puntos"))
  razones <- do.call(paste, c(lapply(partes, `[[`, "texto"), sep = "; "))
  razones <- gsub("(; )+", "; ", razones)
  razones <- gsub("^; |; $", "", razones)
  razones[razones == ""] <- "sin evidencia a favor"

  pts <- function(...) Reduce(`+`, lapply(list(...), function(n) partes[[n]]$puntos))
  x$p_clinvar <- pts("clinvar_p", "clinvar_c")
  x$p_impacto <- pts("alto", "moderado")
  x$p_frecuencia <- pts("gnomad")
  x$p_herencia <- pts("herencia", "portador")
  x$p_fenotipo <- pts("fenotipo")
  x$puntaje <- puntaje
  x$categoria <- ifelse(puntaje >= pesos$corte_alta, "alta",
                        ifelse(puntaje >= pesos$corte_media, "media", "baja"))
  # Un portador en gen recesivo no explica la enfermedad por si solo.
  tope <- h$portador_ar & x$categoria == "alta"
  x$categoria[tope] <- "media"
  razones[tope] <- paste0(razones[tope], "; prioridad limitada a media por ser portador")
  x$razones <- razones
  x$nota <- ifelse(nota == "", NA_character_, nota)
  x$acmg <- .acmg_sugeridos(x, gof, h$compuesto, x$p_fenotipo)

  x <- x[order(x$muestra, -x$puntaje, x$gen), ]
  .poner_attrs(x, a)
}

# Criterios ACMG/AMP que se pueden sugerir con los datos disponibles.
# Siempre requieren revision por un especialista.
.acmg_sugeridos <- function(x, gof, compuesto, p_fenotipo) {
  cons <- ifelse(is.na(x$consecuencia), "", x$consecuencia)
  nula <- grepl("stop_gained|frameshift|splice_acceptor|splice_donor|start_lost", cons)
  af <- x$af_gnomad
  crit <- cbind(
    ifelse(nula & !gof, "PVS1", ""),
    ifelse(is.na(af) | af < 1e-4, "PM2_Supporting", ""),
    ifelse(compuesto, "PM3", ""),
    ifelse(p_fenotipo >= 2, "PP4", ""),
    ifelse(!is.na(af) & af > 0.05, "BA1", "")
  )
  out <- apply(crit, 1, function(f) paste(f[f != ""], collapse = ", "))
  ifelse(out == "", NA_character_, out)
}

# Resume la significancia de ClinVar en: patogenica, conflicto, benigna, otra o NA.
.clasificar_clinvar <- function(cs) {
  vapply(cs, function(s) {
    if (is.na(s) || s == "") return(NA_character_)
    t <- tolower(trimws(unlist(strsplit(s, "[;/,&]"))))
    t <- t[t != ""]
    if (any(grepl("conflicting", t))) return("conflicto")
    pat <- t %in% c("pathogenic", "likely_pathogenic", "pathogenic_low_penetrance",
                    "likely_pathogenic_low_penetrance")
    ben <- t %in% c("benign", "likely_benign")
    if (any(pat) && any(ben)) return("conflicto")
    if (any(pat)) return("patogenica")
    if (any(ben) && all(ben | t %in% c("not_provided", "other"))) return("benigna")
    "otra"
  }, character(1), USE.NAMES = FALSE)
}

.evaluar_herencia <- function(x, herencia) {
  n <- nrow(x)
  compatible <- rep(FALSE, n)
  portador_ar <- rep(FALSE, n)
  texto <- rep("", n)
  nota <- rep("", n)
  g <- x$genotipo
  h <- ifelse(is.na(herencia), "", herencia)

  # Heterocigotos por muestra y gen (para compuestos en AR)
  clave <- paste(x$muestra, x$gen)
  n_het <- stats::ave(as.integer(g == "het"), clave, FUN = sum)

  es <- function(cond, txt) { compatible[cond] <<- TRUE; texto[cond] <<- txt }
  es(h %in% c("AD", "AR/AD") & g %in% c("het", "hom_alt"), "compatible con herencia AD")
  es(h == "AR" & g == "hom_alt", "homocigota en gen AR")
  es(h == "AR/AD" & g == "hom_alt", "homocigota en gen AR/AD")
  es(h == "XL" & g == "hemi", "hemicigota en gen ligado al X")
  es(h == "XL" & g == "hom_alt", "homocigota en gen ligado al X")
  es(h == "MT" & g %in% c("hom_alt", "hemi"), "homoplasmica en gen mitocondrial")

  comp <- h == "AR" & g == "het" & n_het >= 2
  es(comp, "posible heterocigoto compuesto")
  nota[comp] <- "Posible heterocigoto compuesto; confirmar fase (cis o trans) estudiando a los padres."

  unico <- h == "AR" & g == "het" & n_het < 2
  portador_ar[unico] <- TRUE
  nota[unico] <- "Portador en gen recesivo; buscar segunda variante (CNV, intronica) si el fenotipo encaja."

  xl_het <- h == "XL" & g == "het"
  nota[xl_het] <- "Heterocigota en gen ligado al X; portadora si es mujer. Si es varon, revisar sexo indicado."

  sin_h <- h == ""
  nota[sin_h & nota == ""] <- "Gen sin modo de herencia conocido en PanelApp."

  list(compatible = compatible, portador_ar = portador_ar, compuesto = comp,
       texto = texto, nota = nota)
}

.puntaje_fenotipo <- function(genes, fenotipo, categoria, hpo, categorias, maximo) {
  n <- length(genes)
  puntos <- rep(0, n)
  texto <- rep("", n)
  if (!is.null(fenotipo)) {
    terminos <- split(hpo$hpo_id, hpo$gen)
    prop <- vapply(genes, function(g) mean(unique(fenotipo) %in% terminos[[g]]), numeric(1))
    puntos <- round(maximo * prop, 1)
    k <- round(prop * length(unique(fenotipo)))
    texto <- ifelse(puntos > 0,
                    sprintf("fenotipo %d de %d terminos HPO (%+g)", k, length(unique(fenotipo)), puntos),
                    "")
  }
  if (!is.null(categoria)) {
    en_cat <- genes %in% categorias$gen[categorias$categoria_iuis %in% categoria]
    mejora <- en_cat & puntos < maximo
    puntos[mejora] <- maximo
    texto[mejora] <- sprintf("gen en categoria IUIS sospechada (%+g)", maximo)
  }
  list(puntos = unname(puntos), texto = unname(texto))
}
