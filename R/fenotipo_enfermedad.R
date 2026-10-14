# Similitud del fenotipo del caso con cada enfermedad anotada en HPO. A
# diferencia de la similitud por gen, separa las enfermedades de un mismo gen,
# pesa cada signo por su frecuencia en la enfermedad, castiga los signos que la
# enfermedad excluye (anotaciones NOT) y los que el caso no tiene pero la
# enfermedad suele tener.

# Descendientes de cada termino (incluido el propio), calculados una vez.
.descendientes <- function() {
  o <- .ontologia()
  if (is.null(o$desc)) {
    n <- length(o$ids)
    o$desc <- unname(split(rep(seq_len(n), lengths(o$anc)),
                           factor(unlist(o$anc, use.names = FALSE), levels = seq_len(n))))
  }
  o$desc
}

# Similitud de un termino del caso (indice k) con todos los terminos de la
# ontologia. Recorre sus ancestros de menor a mayor contenido de informacion y
# marca a sus descendientes, asi cada termino queda con el ancestro comun mas
# informativo.
.similitud_vector <- function(k) {
  o <- .ontologia()
  s <- numeric(length(o$ids))
  if (is.na(k) || o$ic[k] <= 0 || !o$fenotipico[k]) return(s)
  d <- .descendientes()
  a <- o$anc[[k]]
  for (j in a[order(o$ic[a])]) s[d[[j]]] <- o$ic[j]
  s <- pmin(1, s / o$ic[k])
  s[!o$fenotipico] <- 0
  s
}

# Maximo de v dentro de cada grupo (los grupos vienen ordenados).
.max_grupo <- function(v, grupo) {
  r <- vapply(split(v, grupo), max, numeric(1))
  unname(r)
}

# Prepara la tabla de enfermedades una vez por sesion para los datos del
# paquete; con otra tabla se prepara en cada llamada.
.enf_preparar <- function(enf = ieiprio::hpo_enfermedades) {
  o <- .ontologia()
  propio <- identical(enf, ieiprio::hpo_enfermedades)
  if (propio && !is.null(o$enf)) return(o$enf)
  idx <- match(enf$hpo_id, o$ids)
  ok <- !is.na(idx) & o$fenotipico[idx %|% 1L]
  enf <- enf[ok, ]
  idx <- idx[ok]
  clave <- paste(enf$gen, enf$enfermedad_id, sep = "|")
  grupo <- factor(clave, levels = unique(clave))
  excl <- enf$excluido | enf$frecuencia %in% 0
  primero <- !duplicated(clave)
  info <- tibble::tibble(
    gen = enf$gen[primero], enfermedad_id = enf$enfermedad_id[primero],
    enfermedad = enf$enfermedad[primero], herencia = enf$herencia[primero],
    terminos = as.integer(tabulate(as.integer(grupo)[!excl], nbins = nlevels(grupo)))
  )
  p <- list(idx = idx, grupo = grupo, excl = excl,
            # Un signo raro pesa menos que uno constante: 0.5 a 1 segun su
            # frecuencia; sin frecuencia conocida pesa 1.
            peso = ifelse(is.na(enf$frecuencia), 1, 0.5 + 0.5 * enf$frecuencia),
            # Si el caso no tiene un signo, castiga segun lo frecuente que es
            # en la enfermedad; sin frecuencia conocida, a la mitad.
            f_ausente = ifelse(is.na(enf$frecuencia), 0.5, enf$frecuencia),
            info = info)
  if (propio) o$enf <- p
  p
}

`%|%` <- function(x, y) ifelse(is.na(x), y, x)

# Puntaje de cada enfermedad (columna) para el caso. Devuelve la matriz de
# similitud por termino, la penalizacion por ausentes y el puntaje final.
.similitud_enfermedades <- function(caso, ausentes = NULL, enf = ieiprio::hpo_enfermedades) {
  o <- .ontologia()
  p <- .enf_preparar(enf)
  ng <- nlevels(p$grupo)
  ci <- match(caso, o$ids)
  m <- matrix(0, length(ci), ng)
  for (t in seq_along(ci)) {
    s <- .similitud_vector(ci[t])[p$idx] * p$peso
    s[p$excl] <- 0
    fila <- .max_grupo(s, p$grupo)
    # Contradiccion: el signo del caso es igual o mas especifico que uno que la
    # enfermedad excluye.
    if (!is.na(ci[t]) && any(p$excl)) {
      es_anc <- logical(length(o$ids))
      es_anc[o$anc[[ci[t]]]] <- TRUE
      contra <- .max_grupo(as.numeric(p$excl & es_anc[p$idx]), p$grupo) > 0
      fila[contra] <- -1
    }
    m[t, ] <- fila
  }
  penal <- rep(0, ng)
  ai <- match(ausentes, o$ids)
  ai <- ai[!is.na(ai)]
  if (length(ai) > 0) {
    d <- .descendientes()
    pen <- vapply(ai, function(a) {
      en <- logical(length(o$ids))
      en[d[[a]]] <- TRUE
      .max_grupo(ifelse(!p$excl & en[p$idx], p$f_ausente, 0), p$grupo)
    }, numeric(ng))
    penal <- 0.5 * rowMeans(matrix(pen, nrow = ng))
  }
  presente <- if (length(ci) > 0) colMeans(m) else rep(0, ng)
  list(matriz = m, penal = penal, puntaje = pmax(0, pmin(1, presente - penal)), info = p$info)
}

# Puntaje de fenotipo por gen con cualquiera de los tres metodos. Devuelve un
# tibble con gen, similitud, enfermedad_id, enfermedad y herencia_enfermedad, y
# como atributo la matriz termino x gen (para la mejor enfermedad del gen).
.fenotipo_por_gen <- function(caso, genes, metodo = c("enfermedad", "semantica", "exacta"),
                              ausentes = NULL, hpo = ieiprio::hpo_genes,
                              enf = ieiprio::hpo_enfermedades) {
  metodo <- match.arg(metodo)
  genes <- unique(genes[!is.na(genes)])
  o <- .ontologia()
  hpo <- hpo[o$fenotipico[o$idx[hpo$hpo_id]] %in% TRUE, ]
  out <- tibble::tibble(gen = genes, similitud = 0, enfermedad_id = NA_character_,
                        enfermedad = NA_character_, herencia_enfermedad = NA_character_)
  mat <- matrix(0, length(caso), length(genes), dimnames = list(caso, genes))
  if (length(genes) == 0 || length(caso) == 0) {
    attr(out, "matriz") <- mat
    return(out)
  }
  usa_gen <- rep(TRUE, length(genes))
  if (metodo == "enfermedad") {
    r <- .similitud_enfermedades(caso, ausentes, enf)
    en_panel <- r$info$gen %in% genes
    if (any(en_panel)) {
      k <- which(en_panel)
      # Por gen, la enfermedad con mayor puntaje; ante empate, la que explica
      # mas terminos del caso.
      orden <- k[order(r$info$gen[k], -r$puntaje[k], -colSums(r$matriz[, k, drop = FALSE] > 0))]
      mejor <- orden[!duplicated(r$info$gen[orden])]
      j <- match(r$info$gen[mejor], genes)
      out$similitud[j] <- r$puntaje[mejor]
      out$enfermedad_id[j] <- r$info$enfermedad_id[mejor]
      out$enfermedad[j] <- r$info$enfermedad[mejor]
      out$herencia_enfermedad[j] <- r$info$herencia[mejor]
      mat[, j] <- pmax(r$matriz[, mejor, drop = FALSE], 0)
      usa_gen[j] <- FALSE
    }
  }
  if (any(usa_gen)) {
    g <- genes[usa_gen]
    if (metodo == "exacta") {
      por_gen <- split(hpo$hpo_id, hpo$gen)
      m <- vapply(g, function(x) as.numeric(caso %in% por_gen[[x]]), numeric(length(caso)))
      m <- matrix(m, nrow = length(caso))
    } else {
      m <- .similitud_genes(caso, g, hpo)
    }
    sim <- colMeans(m)
    if (length(ausentes) > 0) sim <- pmax(0, sim - .penal_ausentes_gen(ausentes, g, hpo))
    out$similitud[usa_gen] <- sim
    mat[, usa_gen] <- m
  }
  out$similitud <- round(out$similitud, 3)
  attr(out, "matriz") <- mat
  out
}

# Sin anotaciones por enfermedad: castiga la mitad por cada signo ausente en el
# caso que el gen tiene anotado (o uno mas especifico).
.penal_ausentes_gen <- function(ausentes, genes, hpo) {
  o <- .ontologia()
  d <- .descendientes()
  ai <- match(ausentes, o$ids)
  ai <- ai[!is.na(ai)]
  if (length(ai) == 0) return(rep(0, length(genes)))
  por_gen <- split(match(hpo$hpo_id, o$ids), hpo$gen)
  vapply(genes, function(g) {
    t <- por_gen[[g]]
    0.5 * mean(vapply(ai, function(a) as.numeric(any(t %in% d[[a]])), numeric(1)))
  }, numeric(1), USE.NAMES = FALSE)
}

# Valor p empirico: cuantas veces un fenotipo al azar del mismo tamano,
# tomado de los signos anotados en el panel, da un puntaje igual o mayor.
.p_valor_fenotipo <- function(obs, caso, genes, metodo, permutaciones, hpo, enf) {
  o <- .ontologia()
  pool <- unique(hpo$hpo_id[o$fenotipico[o$idx[hpo$hpo_id]] %in% TRUE])
  mayores <- numeric(length(genes))
  for (i in seq_len(permutaciones)) {
    azar <- sample(pool, length(caso))
    s <- .fenotipo_por_gen(azar, genes, metodo, NULL, hpo, enf)$similitud
    mayores <- mayores + (s >= obs - 1e-9)
  }
  (mayores + 1) / (permutaciones + 1)
}
