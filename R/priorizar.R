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
    clinvar_una_estrella = 3,
    clinvar_sin_estrellas = 1,
    clinvar_conflicto  = 1,
    impacto_alto       = 3,
    impacto_moderado   = 1,
    ausente_gnomad     = 1,
    predictor_a_favor  = 1,
    predictor_en_contra = -1,
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
#' * **ClinVar.** Patogenica o probablemente patogenica suma 4 con dos o mas
#'   estrellas de revision (o si no se pudieron verificar), 3 con una estrella
#'   y 1 sin estrellas. Con interpretaciones en conflicto suma 1.
#' * **Impacto.** Alto (stop, frameshift, splicing canonico) suma 3.
#'   Moderado (missense, inframe) suma 1. En genes con enfermedad por
#'   ganancia de funcion una missense cuenta como impacto alto.
#' * **Frecuencia.** Ausente de gnomAD suma 1.
#' * **Predictores.** REVEL, AlphaMissense o SpliceAI por encima de sus
#'   umbrales calibrados suman 1 (PP3). Si todos los disponibles indican
#'   efecto benigno restan 1 (BP4).
#' * **Herencia.** Genotipo compatible con el modo del gen suma 2. En genes
#'   autosomicos recesivos, dos heterocigotos en el mismo gen y la misma
#'   muestra se marcan como posible heterocigoto compuesto y suman 2; un
#'   heterocigoto solo resta 2, se marca como portador y nunca queda en
#'   prioridad alta, aunque el fenotipo encaje.
#' * **Fenotipo.** Parecido entre el fenotipo del caso y el del gen,
#'   escalado de 0 a 3. Cada termino del caso se compara con el mas parecido
#'   usando la jerarquia de HPO: un termino igual o mas especifico vale 1 y
#'   uno emparentado vale menos segun cuanta informacion comparten, asi
#'   "neumonia recurrente" suma para un gen anotado con "infecciones
#'   respiratorias recurrentes". Con `similitud = "enfermedad"` (por
#'   defecto) la comparacion se hace con cada enfermedad del gen por
#'   separado, pesando cada signo por su frecuencia y restando los signos que
#'   la enfermedad excluye y los `fenotipo_ausente`; las razones nombran la
#'   enfermedad mas parecida. Con `"semantica"` se usan todos los terminos del
#'   gen juntos y con `"exacta"` solo los identicos (ver
#'   [iei_buscar_fenotipo()]). Un gen en la `categoria` IUIS indicada recibe
#'   3 puntos.
#'
#' @param x Tabla devuelta por [iei_anotar()] o [iei_filtrar_frecuencia()].
#' @param fenotipo Fenotipo del caso, como codigos HPO
#'   (`c("HP:0004313", "HP:0002719")`) o como texto en espanol o ingles
#'   (`c("agammaglobulinemia", "neumonia recurrente")`). El texto se traduce a
#'   HPO y se informa la traduccion; ver [iei_buscar_hpo()].
#' @param similitud `"enfermedad"`, `"semantica"` o `"exacta"`, como se
#'   compara el fenotipo con el del gen.
#' @param fenotipo_ausente Signos que el caso no tiene, en el mismo formato
#'   que `fenotipo`. Restan a los genes cuyas enfermedades suelen tenerlos.
#' @param categoria Categorias IUIS sospechadas (1 a 10), alternativa a
#'   `fenotipo`.
#' @param pesos Lista de pesos, ver [iei_pesos()].
#' @param panel,hpo,categorias,enfermedades Tablas de referencia. Por
#'   defecto los datos del paquete.
#'
#' @return La tabla con el desglose del puntaje (`p_clinvar`, `p_impacto`,
#'   `p_frecuencia`, `p_predictores`, `p_herencia`, `p_fenotipo`), el `puntaje` total, la
#'   `categoria`, las `razones`, una `nota` clinica, `acmg` con los
#'   criterios ACMG/AMP sugeridos (PVS1 graduado, PS1, PM1, PM2, PM3, PM4,
#'   PM5, PP2, PP3, PP4, BA1, BS1, BS2, BP4 y BP7, segun
#'   los datos disponibles), `acmg_puntos` con su suma en el sistema de
#'   puntos de Tavtigian y colaboradores (2020) y `acmg_clase` con la clase
#'   que corresponde. Los criterios son orientativos y deben ser revisados
#'   por un especialista. Ordenada por muestra y puntaje. Es un objeto
#'   `iei_resultado`, que se comporta como tibble y tiene metodos
#'   [print()], [summary()] y [plot()]; ver tambien [iei_tabla()] y
#'   [iei_metadatos()].
#' @export
#' @examples
#' \dontrun{
#' prior <- iei_priorizar(anot, fenotipo = c("HP:0004313", "HP:0002719"))
#' }
iei_priorizar <- function(x, fenotipo = NULL, categoria = NULL, pesos = iei_pesos(),
                          similitud = c("enfermedad", "semantica", "exacta"), fenotipo_ausente = NULL,
                          panel = ieiprio::iuis_panel, hpo = ieiprio::hpo_genes,
                          categorias = ieiprio::iuis_categorias,
                          enfermedades = ieiprio::hpo_enfermedades) {
  necesarias <- c("muestra", "gen", "genotipo", "impacto", "consecuencia",
                  "af_gnomad", "clinvar")
  faltan <- setdiff(necesarias, names(x))
  if (length(faltan) > 0) {
    cli::cli_abort(c("A {.arg x} le faltan columnas: {.val {faltan}}.",
                     "i" = "Pasala antes por {.fn iei_filtrar} y {.fn iei_anotar}."))
  }
  similitud <- match.arg(similitud)
  fenotipo <- .fenotipo_codigos(fenotipo)
  fenotipo_ausente <- .fenotipo_codigos(fenotipo_ausente)
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
  a$ieiprio_similitud <- if (!is.null(fenotipo)) similitud
  a$ieiprio_fenotipo_ausente <- fenotipo_ausente
  a$ieiprio_pesos <- pesos

  if (nrow(x) == 0) {
    for (v in c("p_clinvar", "p_impacto", "p_frecuencia", "p_predictores", "p_herencia",
                "p_fenotipo", "puntaje"))
      x[[v]] <- numeric(0)
    for (v in c("categoria", "razones", "nota", "acmg", "acmg_clase")) x[[v]] <- character(0)
    x$acmg_puntos <- numeric(0)
    return(.como_resultado(.poner_attrs(x, a)))
  }

  info <- panel[match(x$gen, panel$gen), ]
  herencia <- info$herencia
  gof <- !is.na(info$ganancia_funcion) & info$ganancia_funcion
  # PVS1 solo se descarta cuando la enfermedad del gen es exclusivamente por
  # ganancia de funcion; si tambien hay enfermedad por perdida, una nula vale.
  gof_puro <- if ("mecanismo" %in% names(info)) info$mecanismo %in% "GOF" else gof

  partes <- list()
  nota <- rep("", nrow(x))
  agregar <- function(cond, puntos, texto) {
    list(puntos = ifelse(cond, puntos, 0), texto = ifelse(cond, sprintf("%s (%+g)", texto, puntos), ""))
  }

  est <- .col(x, "clinvar_estrellas", NA_integer_)
  pts_cv <- ifelse(is.na(est) | est >= 2, pesos$clinvar_patogenica,
                   ifelse(est == 1, pesos$clinvar_una_estrella, pesos$clinvar_sin_estrellas))
  txt_cv <- ifelse(is.na(est), "ClinVar patog\u00e9nica",
                   sprintf("ClinVar patog\u00e9nica, %d estrella%s", est, ifelse(est == 1, "", "s")))
  partes$clinvar_p <- list(puntos = ifelse(cv %in% "patogenica", pts_cv, 0),
                           texto = ifelse(cv %in% "patogenica", sprintf("%s (%+g)", txt_cv, pts_cv), ""))
  partes$clinvar_c <- agregar(cv %in% "conflicto", pesos$clinvar_conflicto, "ClinVar en conflicto")

  missense <- grepl("missense", x$consecuencia %||% "", fixed = TRUE)
  alto <- (!is.na(x$impacto) & x$impacto == "HIGH") | (gof & missense)
  moderado <- !alto & !is.na(x$impacto) & x$impacto == "MODERATE"
  partes$alto <- agregar(alto, pesos$impacto_alto,
                         ifelse(gof & missense & x$impacto %in% "MODERATE",
                                "missense en gen de ganancia de funci\u00f3n", "impacto alto"))
  partes$moderado <- agregar(moderado, pesos$impacto_moderado, "impacto moderado")
  partes$gnomad <- agregar(is.na(.af_referencia(x)), pesos$ausente_gnomad, "ausente en gnomAD")

  # Predictores in silico (PP3 y BP4)
  pr <- .evaluar_predictores(x)
  partes$pred_favor <- agregar(pr$favor, pesos$predictor_a_favor, pr$texto_favor)
  partes$pred_contra <- agregar(pr$contra, pesos$predictor_en_contra, "predictores sugieren efecto benigno")

  # Herencia
  h <- .evaluar_herencia(x, herencia)
  partes$herencia <- agregar(h$compatible, pesos$herencia_compatible, h$texto)
  partes$portador <- agregar(h$portador_ar, pesos$portador_ar, "heterocigoto \u00fanico en gen recesivo")
  nota <- h$nota
  # Advertencias que no cambian el puntaje pero si la lectura del resultado.
  dif <- ieiprio::genes_dificiles
  dif <- dif[dif$tipo != "mosaico", ]
  k <- match(x$gen, dif$gen)
  nota <- ifelse(!is.na(k), ifelse(nota == "", dif$motivo[k], paste(nota, dif$motivo[k])), nota)
  if ("mosaico_posible" %in% names(x)) {
    mos <- x$mosaico_posible %in% TRUE
    txt <- sprintf("Posible mosaico (fracci\u00f3n al\u00e9lica %.2f); confirmar con otra t\u00e9cnica o tejido.", .col(x, "ab"))
    nota[mos] <- ifelse(nota[mos] == "", txt[mos], paste(nota[mos], txt[mos]))
  }

  # Fenotipo
  f <- .puntaje_fenotipo(x$gen, fenotipo, categoria, hpo, categorias, pesos$fenotipo_maximo,
                         similitud, fenotipo_ausente, enfermedades)
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
  x$p_predictores <- pts("pred_favor", "pred_contra")
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
  # PP4 pide un fenotipo muy especifico: solo con fenotipo HPO (no con la
  # categoria IUIS, que es amplia) y similitud de 0.85 o mas.
  p_pp4 <- if (is.null(fenotipo)) rep(0, nrow(x)) else x$p_fenotipo / pesos$fenotipo_maximo * 3
  ac <- .acmg_sugeridos(x, gof_puro, h, p_pp4, pr, herencia)
  x$acmg <- ac$criterios
  x$acmg_puntos <- ac$puntos
  x$acmg_clase <- ac$clase

  x <- x[order(x$muestra, -x$puntaje, x$gen), ]
  .como_resultado(.poner_attrs(x, a))
}

# Criterios ACMG/AMP que se pueden sugerir con los datos disponibles, con su
# fuerza, y clasificacion por el sistema de puntos de Tavtigian y
# colaboradores (ClinGen): muy fuerte 8, fuerte 4, moderado 2, de apoyo 1 y
# los benignos en negativo. Siempre requieren revision por un especialista.
.fuerza_puntos <- c(VeryStrong = 8, Strong = 4, Moderate = 2, Supporting = 1)

.puntos_criterio <- function(cr) {
  vapply(cr, USE.NAMES = FALSE, function(c) {
    base <- sub("_.*$", "", c)
    fuerza <- if (grepl("_", c)) sub("^[^_]+_", "", c) else NA
    defecto <- if (startsWith(base, "PVS")) 8 else if (startsWith(base, "PS")) 4 else
      if (startsWith(base, "PM")) 2 else if (startsWith(base, "PP")) 1 else
      if (startsWith(base, "BA")) -8 else if (startsWith(base, "BS")) -4 else if (startsWith(base, "BP")) -1 else 0
    v <- if (!is.na(fuerza) && fuerza %in% names(.fuerza_puntos)) .fuerza_puntos[[fuerza]] else abs(defecto)
    if (startsWith(base, "B")) -v else v
  }, numeric(1))
}

.clase_acmg <- function(puntos, ba1) {
  clase <- ifelse(puntos >= 10, "Patog\u00e9nica",
           ifelse(puntos >= 6, "Probablemente patog\u00e9nica",
           ifelse(puntos >= 0, "Significado incierto",
           ifelse(puntos >= -6, "Probablemente benigna", "Benigna"))))
  clase[ba1] <- "Benigna"
  clase
}

.acmg_sugeridos <- function(x, gof, h, p_fenotipo, pr = NULL, herencia = NULL) {
  n <- nrow(x)
  if (is.null(pr)) pr <- list(favor = rep(FALSE, n), contra = rep(FALSE, n))
  if (is.logical(h)) h <- list(compuesto = h)
  compuesto <- h$compuesto %||% rep(FALSE, n)
  trans <- h$trans %||% rep(FALSE, n)
  her <- if (is.null(herencia)) rep(NA_character_, n) else herencia
  cons <- ifelse(is.na(x$consecuencia), "", x$consecuencia)
  af <- .af_referencia(x)
  exon <- .col(x, "exon", NA_character_)
  ex_k <- suppressWarnings(as.integer(sub("/.*$", "", exon)))
  ex_n <- suppressWarnings(as.integer(sub("^.*/", "", exon)))
  ultimo_exon <- !is.na(ex_k) & !is.na(ex_n) & ex_k == ex_n & ex_n > 1

  # PVS1 graduado (simplificacion del arbol de Abou Tayoun 2018).
  trunc <- grepl("stop_gained|frameshift", cons)
  splice <- grepl("splice_acceptor|splice_donor", cons)
  inicio <- grepl("start_lost", cons)
  pvs1 <- ifelse(gof, "",
          ifelse(trunc & !ultimo_exon, "PVS1",
          ifelse(trunc & ultimo_exon, "PVS1_Moderate",
          ifelse(splice, "PVS1_Strong",
          ifelse(inicio, "PVS1_Moderate", "")))))

  dominante <- her %in% c("AD", "XL")
  pm2 <- (is.na(af) | (dominante & af <= 1e-5) | (!dominante & af <= 1e-4))
  hom <- .col(x, "gnomad_hom", NA_integer_); hemi <- .col(x, "gnomad_hemi", NA_integer_)
  bs2 <- (her %in% "AR" & !is.na(hom) & hom >= 2) | (her %in% "XL" & !is.na(hemi) & hemi >= 2)
  bs1 <- !is.na(af) & ((dominante & af > 1e-3) | (!dominante & af > 0.02))
  ba1 <- !is.na(af) & af > 0.05
  missense <- grepl("missense", cons)
  mis_z <- .col(x, "mis_z")
  sai <- .col(x, "spliceai")
  bp7 <- grepl("synonymous", cons) & !is.na(sai) & sai < 0.1

  crit <- cbind(
    pvs1,
    ifelse(.col(x, "ps1", FALSE) %in% TRUE, "PS1", ""),
    ifelse(.col(x, "pm1", FALSE) %in% TRUE, "PM1_Supporting", ""),
    ifelse(pm2 & !ba1, "PM2_Supporting", ""),
    ifelse(compuesto & trans, "PM3", ifelse(compuesto, "PM3_Supporting", "")),
    ifelse(grepl("inframe_insertion|inframe_deletion|stop_lost", cons), "PM4", ""),
    ifelse(.col(x, "pm5", FALSE) %in% TRUE & !(.col(x, "ps1", FALSE) %in% TRUE), "PM5", ""),
    ifelse(missense & !is.na(mis_z) & mis_z >= 3.09, "PP2", ""),
    ifelse(pr$favor, "PP3", ""),
    ifelse(p_fenotipo >= 0.85 * 3, "PP4", ""),
    ifelse(ba1, "BA1", ""),
    ifelse(bs1 & !ba1, "BS1", ""),
    ifelse(bs2, "BS2", ""),
    ifelse(pr$contra, "BP4", ""),
    ifelse(bp7, "BP7", "")
  )
  criterios <- apply(crit, 1, function(f) paste(f[f != ""], collapse = ", "))
  puntos <- apply(crit, 1, function(f) sum(.puntos_criterio(f[f != ""])))
  list(criterios = ifelse(criterios == "", NA_character_, criterios),
       puntos = unname(puntos), clase = .clase_acmg(unname(puntos), ba1))
}

# Umbrales de predictores para PP3 y BP4 a nivel de apoyo:
# REVEL segun la calibracion de ClinGen (Pejaver et al. 2022), SpliceAI segun
# las recomendaciones de ClinGen SVI para splicing (Walker et al. 2023) y
# AlphaMissense segun los cortes de sus autores (Cheng et al. 2023).
.umbrales_predictores <- list(revel_favor = 0.644, revel_contra = 0.290,
                              am_favor = 0.564, am_contra = 0.34,
                              spliceai_favor = 0.2, spliceai_contra = 0.1)

.col <- function(x, n, tipo = NA_real_) if (n %in% names(x)) x[[n]] else rep(tipo, nrow(x))

.evaluar_predictores <- function(x, u = .umbrales_predictores) {
  revel <- .col(x, "revel"); am <- .col(x, "alphamissense"); sai <- .col(x, "spliceai")
  f_revel <- !is.na(revel) & revel >= u$revel_favor
  f_am <- !is.na(am) & am >= u$am_favor
  f_sai <- !is.na(sai) & sai >= u$spliceai_favor
  favor <- f_revel | f_am | f_sai
  texto <- vapply(seq_len(nrow(x)), function(i) {
    q <- c(if (f_revel[i]) sprintf("REVEL %.2f", revel[i]),
           if (f_am[i]) sprintf("AlphaMissense %.2f", am[i]),
           if (f_sai[i]) sprintf("SpliceAI %.2f", sai[i]))
    if (length(q) == 0) "" else paste0("predictores a favor (", paste(q, collapse = ", "), ")")
  }, character(1))
  hay <- !is.na(revel) | !is.na(am) | !is.na(sai)
  benigno <- (is.na(revel) | revel < u$revel_contra) & (is.na(am) | am < u$am_contra) &
    (is.na(sai) | sai < u$spliceai_contra)
  contra <- hay & benigno & !favor
  list(favor = favor, contra = contra, texto_favor = texto)
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
  g <- ifelse(is.na(x$genotipo), "desconocido", x$genotipo)
  h <- ifelse(is.na(herencia), "", herencia)

  # Variantes que pueden ser la segunda de un compuesto: no basta cualquier
  # heterocigota, tiene que tener algun indicio de dano.
  sai <- .col(x, "spliceai")
  relevante <- .col(x, "impacto", NA_character_) %in% c("HIGH", "MODERATE") |
    grepl("pathogenic", tolower(.col(x, "clinvar", NA_character_))) |
    (!is.na(sai) & sai >= 0.2)
  recesivo <- h %in% c("AR", "AR/AD")
  clave <- paste(x$muestra, x$gen)
  het_rel <- g == "het" & recesivo & relevante
  n_het <- stats::ave(as.integer(het_rel), clave, FUN = sum)

  es <- function(cond, txt) { compatible[cond] <<- TRUE; texto[cond] <<- txt }
  es(h %in% c("AD", "AR/AD") & g %in% c("het", "hom_alt"), "compatible con herencia AD")
  es(h == "AR" & g == "hom_alt", "homocigota en gen AR")
  es(h == "AR/AD" & g == "hom_alt", "homocigota en gen AR/AD")
  es(h == "XL" & g == "hemi", "hemicigota en gen ligado al X")
  es(h == "XL" & g == "hom_alt", "homocigota en gen ligado al X")
  es(h == "MT" & g %in% c("hom_alt", "hemi"), "homopl\u00e1smica en gen mitocondrial")

  # Fase de cada par: por la fase del propio VCF (GT con "|" y PS), o desconocida.
  fase <- .fase_compuestos(x, het_rel & n_het >= 2)
  comp <- het_rel & n_het >= 2 & fase != "cis"
  es(comp & h == "AR", "posible heterocigoto compuesto")
  es(comp & fase == "trans", "heterocigoto compuesto en trans")
  nota[comp & fase == "desconocida"] <- "Posible heterocigoto compuesto; confirmar la fase (cis o trans) con lecturas que cubran ambas variantes o con un VCF faseado."
  nota[comp & fase == "trans"] <- "Heterocigoto compuesto con las variantes en alelos distintos (trans)."
  cis <- het_rel & n_het >= 2 & fase == "cis"
  nota[cis] <- "Las dos variantes del gen est\u00e1n en el mismo alelo (cis); no forman un compuesto."

  unico <- h == "AR" & g == "het" & !comp
  portador_ar[unico] <- TRUE
  nota[unico & nota == ""] <- "Portador en gen recesivo; buscar segunda variante (CNV, intr\u00f3nica) si el fenotipo encaja."

  xl_het <- h == "XL" & g == "het"
  nota[xl_het] <- "Heterocigota en gen ligado al X; portadora si es mujer. Si es var\u00f3n, revisar el sexo indicado."

  sin_h <- h == ""
  nota[sin_h & nota == ""] <- "Gen sin modo de herencia conocido en PanelApp."
  sin_g <- g == "desconocido"
  nota[sin_g & nota == ""] <- "Genotipo no indicado; no se evalu\u00f3 la compatibilidad con la herencia."

  list(compatible = compatible, portador_ar = portador_ar, compuesto = comp,
       trans = comp & fase == "trans", texto = texto, nota = nota)
}

# Fase de las heterocigotas de un mismo gen y muestra segun el propio VCF:
# "trans" si GT y PS las ponen en haplotipos distintos del mismo bloque, "cis"
# si estan en el mismo haplotipo; si no, "desconocida".
.fase_compuestos <- function(x, candidatas) {
  n <- nrow(x)
  fase <- rep("desconocida", n)
  gt <- .col(x, "gt", NA_character_); ps <- .col(x, "ps", NA_character_)
  clave <- paste(x$muestra, x$gen)
  for (cl in unique(clave[candidatas])) {
    idx <- which(clave == cl & candidatas)
    faseado <- !is.na(gt[idx]) & grepl("|", gt[idx], fixed = TRUE) & !is.na(ps[idx])
    if (all(faseado) && length(unique(ps[idx])) == 1) {
      hap <- vapply(strsplit(gt[idx], "|", fixed = TRUE), function(a) which(a != "0")[1], integer(1))
      fase[idx] <- if (length(unique(hap)) > 1) "trans" else "cis"
    }
  }
  fase
}

.puntaje_fenotipo <- function(genes, fenotipo, categoria, hpo, categorias, maximo,
                              similitud = "exacta", ausentes = NULL,
                              enfermedades = ieiprio::hpo_enfermedades) {
  n <- length(genes)
  puntos <- rep(0, n)
  texto <- rep("", n)
  if (!is.null(fenotipo)) {
    fenotipo <- unique(fenotipo)
    terminos <- split(hpo$hpo_id, hpo$gen)
    exactos <- vapply(genes, function(g) sum(fenotipo %in% terminos[[g]]), numeric(1))
    r <- .fenotipo_por_gen(fenotipo, genes, similitud, ausentes, hpo, enfermedades)
    k <- match(genes, r$gen)
    prop <- r$similitud[k]
    prop[is.na(prop)] <- 0
    puntos <- round(maximo * prop, 1)
    enf <- r$enfermedad[k]
    detalle <- ifelse(is.na(enf), "", sprintf(" con %s (%s)", trimws(enf), r$enfermedad_id[k]))
    aus <- if (length(ausentes) > 0) ", descontando los signos ausentes" else ""
    texto <- switch(similitud,
      exacta = ifelse(puntos > 0, sprintf("fenotipo %d de %d t\u00e9rminos HPO%s (%+g)", exactos,
                                          length(fenotipo), aus, puntos), ""),
      ifelse(puntos > 0, sprintf("fenotipo parecido%s, similitud %.2f y %d de %d t\u00e9rminos exactos%s (%+g)",
                                 detalle, prop, exactos, length(fenotipo), aus, puntos), ""))
  }
  if (!is.null(categoria)) {
    en_cat <- genes %in% categorias$gen[categorias$categoria_iuis %in% categoria]
    mejora <- en_cat & puntos < maximo
    puntos[mejora] <- maximo
    texto[mejora] <- sprintf("gen en categor\u00eda IUIS sospechada (%+g)", maximo)
  }
  list(puntos = unname(puntos), texto = unname(texto))
}
