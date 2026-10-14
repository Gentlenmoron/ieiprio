# Ontologia HPO: terminos, jerarquia, contenido de informacion y similitud.

# ---- Construccion de los datos (data-raw) ---------------------------------

# Lee hp.obo y devuelve una fila por termino vigente con sus padres (is_a),
# sinonimos exactos e identificadores alternativos.
.obo_procesar <- function(ruta) {
  lineas <- readLines(ruta, warn = FALSE, encoding = "UTF-8")
  inicio <- which(lineas == "[Term]")
  fin <- c(inicio[-1] - 1L, length(lineas))
  bloques <- mapply(function(a, b) lineas[a:b], inicio, fin, SIMPLIFY = FALSE)
  valor <- function(b, clave) {
    l <- b[startsWith(b, paste0(clave, ": "))]
    substring(l, nchar(clave) + 3L)
  }
  filas <- lapply(bloques, function(b) {
    if (length(valor(b, "is_obsolete")) > 0) return(NULL)
    id <- valor(b, "id")
    if (length(id) != 1 || !startsWith(id, "HP:")) return(NULL)
    sin <- valor(b, "synonym")
    sin <- sin[grepl(" EXACT", sin, fixed = TRUE)]
    list(id = id, nombre = valor(b, "name")[1],
         padres = sub(" !.*$", "", valor(b, "is_a")),
         sinonimos = sub('^"(.*)" EXACT.*$', "\\1", sin),
         alternos = valor(b, "alt_id"))
  })
  filas <- Filter(Negate(is.null), filas)
  tibble::tibble(
    id = vapply(filas, `[[`, "", "id"),
    nombre = vapply(filas, `[[`, "", "nombre"),
    padres = lapply(filas, `[[`, "padres"),
    sinonimos = lapply(filas, `[[`, "sinonimos"),
    alternos = lapply(filas, `[[`, "alternos")
  )
}

# Ancestros de cada termino (incluido el propio), por propagacion sobre is_a.
.ancestros_todos <- function(ids, padres) {
  todos <- match(unlist(padres, use.names = FALSE), ids)
  grupo <- rep(seq_along(padres), lengths(padres))
  ok <- !is.na(todos)
  par_i <- unname(split(todos[ok], factor(grupo[ok], levels = seq_along(ids))))
  # Orden topologico (Kahn): cada termino se procesa despues de sus padres.
  n <- length(ids)
  pendientes <- lengths(par_i)
  hijos <- split(rep(seq_len(n), lengths(par_i)), factor(unlist(par_i), levels = seq_len(n)))
  orden <- integer(n)
  raices <- which(pendientes == 0)
  orden[seq_along(raices)] <- raices
  leido <- 0L; k <- length(raices)
  while (leido < k) {
    leido <- leido + 1L
    h <- hijos[[orden[leido]]]
    if (length(h) == 0) next
    pendientes[h] <- pendientes[h] - 1L
    listos <- h[pendientes[h] == 0L]
    if (length(listos)) { orden[k + seq_along(listos)] <- listos; k <- k + length(listos) }
  }
  anc <- vector("list", n)
  for (i in orden[seq_len(k)]) {
    anc[[i]] <- if (length(par_i[[i]]) == 0) i else unique(c(i, unlist(anc[par_i[[i]]], use.names = FALSE)))
  }
  anc
}

# Contenido de informacion: -log de la fraccion de genes anotados con el
# termino o alguno de sus descendientes, sobre todas las anotaciones de HPO.
.ic_calcular <- function(terminos, ruta_g2p) {
  g2p <- utils::read.delim(ruta_g2p, quote = "", stringsAsFactors = FALSE,
                           colClasses = "character")[, c("gene_symbol", "hpo_id")]
  g2p <- unique(g2p[g2p$hpo_id %in% terminos$id, ])
  anc <- .ancestros_todos(terminos$id, terminos$padres)
  idx <- match(g2p$hpo_id, terminos$id)
  pares <- unique(data.frame(gen = rep(g2p$gene_symbol, lengths(anc[idx])),
                             termino = unlist(anc[idx])))
  n_genes <- length(unique(g2p$gene_symbol))
  conteo <- tabulate(pares$termino, nbins = nrow(terminos))
  ic <- -log(pmax(conteo, 1) / n_genes)
  ic[conteo == 0] <- -log(1 / n_genes)   # sin anotaciones: tan informativo como el mas raro
  round(ic, 4)
}

# ---- Uso en tiempo de ejecucion -------------------------------------------

.ont <- new.env(parent = emptyenv())

# Indices y ancestros calculados una sola vez por sesion.
.ontologia <- function() {
  if (is.null(.ont$listo)) {
    t <- ieiprio::hpo_terminos
    .ont$ids <- t$id
    .ont$idx <- stats::setNames(seq_along(t$id), t$id)
    .ont$ic <- t$ic
    .ont$nombre <- stats::setNames(t$nombre, t$id)
    .ont$anc <- .ancestros_todos(t$id, t$padres)
    # Solo los terminos bajo "Phenotypic abnormality" describen signos; los de
    # modo de herencia, curso clinico o frecuencia no deben entrar en la
    # similitud.
    raiz <- .ont$idx[["HP:0000118"]]
    .ont$fenotipico <- vapply(.ont$anc, function(a) raiz %in% a, logical(1))
    alt <- unlist(t$alternos)
    .ont$alternos <- stats::setNames(rep(t$id, lengths(t$alternos)), alt)
    .ont$listo <- TRUE
  }
  .ont
}

# Pasa identificadores antiguos a su termino vigente; NA si no existe.
.hpo_vigente <- function(ids) {
  o <- .ontologia()
  ok <- ids %in% o$ids
  ids[!ok] <- unname(o$alternos[ids[!ok]])
  ids
}

.hpo_nombre <- function(ids) unname(.ontologia()$nombre[ids])

# Nombres y sinonimos sin tildes ni mayusculas, para buscar por texto.
.ont_textos <- function() {
  o <- .ontologia()
  if (is.null(o$nombres_norm)) {
    t <- ieiprio::hpo_terminos
    o$nombres_norm <- .sin_tildes(t$nombre)
    sin <- .sin_tildes(unlist(t$sinonimos, use.names = FALSE))
    o$sinon_norm <- unname(split(sin, factor(rep(seq_along(t$id), lengths(t$sinonimos)), levels = seq_along(t$id))))
  }
  o
}

# Para cada termino del caso, la similitud con el mejor termino del gen:
# contenido de informacion del ancestro comun mas informativo dividido por el
# del termino del caso. Vale 1 si el gen tiene el termino o uno mas especifico
# y baja a medida que el parecido es mas general.
.similitud_terminos <- function(caso, terminos_gen) {
  o <- .ontologia()
  .similitud_idx(match(caso, o$ids), match(terminos_gen, o$ids))
}

# Nucleo con indices enteros (evita buscar nombres en cada llamada).
.similitud_idx <- function(ci, gi) {
  o <- .ontologia()
  gi <- gi[!is.na(gi)]
  gi <- gi[o$fenotipico[gi]]
  if (length(gi) == 0) return(rep(0, length(ci)))
  marca <- unlist(o$anc[gi], use.names = FALSE)
  vapply(ci, function(k) {
    if (is.na(k) || o$ic[k] <= 0 || !o$fenotipico[k]) return(0)
    a <- o$anc[[k]]
    a <- a[a %in% marca]
    if (length(a) == 0) 0 else min(1, max(o$ic[a]) / o$ic[k])
  }, numeric(1))
}

# Similitud de un conjunto de terminos con cada gen. Devuelve una matriz con
# una fila por termino del caso y una columna por gen.
.similitud_genes <- function(caso, genes, hpo = ieiprio::hpo_genes) {
  o <- .ontologia()
  ci <- match(caso, o$ids)
  por_gen <- split(match(hpo$hpo_id, o$ids), hpo$gen)
  m <- vapply(genes, function(g) .similitud_idx(ci, por_gen[[g]] %||% integer(0)),
              numeric(length(caso)))
  if (length(caso) == 1) m <- matrix(m, nrow = 1)
  dimnames(m) <- list(caso, genes)
  m
}

# ---- Texto a terminos HPO --------------------------------------------------

# Terminos frecuentes en errores innatos de la inmunidad, en espanol.
.hpo_espanol <- c(
  "infecciones recurrentes" = "HP:0002719",
  "infeccion recurrente" = "HP:0002719",
  "hipogammaglobulinemia" = "HP:0004313",
  "inmunoglobulinas bajas" = "HP:0004313",
  "agammaglobulinemia" = "HP:0004432",
  "neumonia" = "HP:0002090",
  "neumonia recurrente" = "HP:0006532",
  "neumonias recurrentes" = "HP:0006532",
  "infecciones respiratorias recurrentes" = "HP:0002205",
  "otitis media recurrente" = "HP:0000403",
  "sinusitis recurrente" = "HP:0011108",
  "bronquiectasias" = "HP:0002110",
  "linfopenia" = "HP:0001888",
  "neutropenia" = "HP:0001875",
  "trombocitopenia" = "HP:0001873",
  "anemia" = "HP:0001903",
  "eccema" = "HP:0000964",
  "dermatitis atopica" = "HP:0001047",
  "diarrea" = "HP:0002014",
  "diarrea cronica" = "HP:0002028",
  "falla para crecer" = "HP:0001508",
  "retraso del crecimiento" = "HP:0001508",
  "esplenomegalia" = "HP:0001744",
  "hepatomegalia" = "HP:0002240",
  "linfadenopatia" = "HP:0002716",
  "adenopatias" = "HP:0002716",
  "candidiasis" = "HP:0002728",
  "candidiasis mucocutanea cronica" = "HP:0002728",
  "fiebre recurrente" = "HP:0001954",
  "fiebre periodica" = "HP:0001954",
  "autoinmunidad" = "HP:0002960",
  "citopenia autoinmune" = "HP:5210419",
  "anemia hemolitica autoinmune" = "HP:0001890",
  "trombocitopenia autoinmune" = "HP:0001973",
  "enfermedad inflamatoria intestinal" = "HP:0002037",
  "colitis" = "HP:0002583",
  "enfermedad de crohn" = "HP:0100280",
  "abscesos recurrentes" = "HP:0002722",
  "infecciones micobacterianas" = "HP:0011274",
  "micobacterias" = "HP:0011274",
  "ige elevada" = "HP:0003212",
  "hiper ige" = "HP:0003212",
  "igm elevada" = "HP:0003496",
  "iga disminuida" = "HP:0002720",
  "celulas b ausentes" = "HP:0030252",
  "linfocitos b ausentes" = "HP:0030252",
  "celulas b disminuidas" = "HP:0010976",
  "celulas t disminuidas" = "HP:0005403",
  "hemofagocitosis" = "HP:0012156",
  "granulomas" = "HP:0002955",
  "sepsis" = "HP:0100806",
  "meningitis" = "HP:0001287",
  "artritis" = "HP:0001369",
  "vasculitis" = "HP:0002633",
  "verrugas" = "HP:0200043"
)

.sin_tildes <- function(x) {
  x <- as.character(x)
  de <- c("\u00e1", "\u00e9", "\u00ed", "\u00f3", "\u00fa", "\u00fc", "\u00f1",
          "\u00c1", "\u00c9", "\u00cd", "\u00d3", "\u00da", "\u00d1")
  a <- c("a", "e", "i", "o", "u", "u", "n", "a", "e", "i", "o", "u", "n")
  for (k in seq_along(de)) {
    x <- gsub(de[k], a[k], x, fixed = TRUE, useBytes = TRUE)
    # En sesiones sin UTF-8, R puede convertir la letra en el texto "<U+00ED>".
    x <- gsub(sprintf("<U+%04X>", utf8ToInt(de[k])), a[k], x, fixed = TRUE)
  }
  tolower(x)
}

# Convierte terminos HPO o texto (espanol o ingles) en identificadores HPO.
# Devuelve un tibble con lo escrito, el termino elegido y como se encontro.
.interpretar_fenotipo <- function(x) {
  o <- .ontologia()
  t <- ieiprio::hpo_terminos
  if (!all(grepl("^HP:[0-9]{7}$", trimws(x)))) .ont_textos()
  nombres <- o$nombres_norm
  sinon <- o$sinon_norm
  filas <- lapply(x, function(e) {
    e0 <- trimws(e)
    if (grepl("^HP:[0-9]{7}$", e0)) {
      v <- .hpo_vigente(e0)
      return(list(entrada = e, hpo_id = v, como = if (is.na(v)) "no existe" else if (v == e0) "codigo" else "codigo antiguo"))
    }
    q <- .sin_tildes(e0)
    if (!is.na(.hpo_espanol[q])) return(list(entrada = e, hpo_id = unname(.hpo_espanol[q]), como = "diccionario"))
    k <- which(nombres == q)
    if (length(k) == 0) k <- which(vapply(sinon, function(s) q %in% s, logical(1)))
    if (length(k) > 0) return(list(entrada = e, hpo_id = t$id[k[1]], como = "nombre exacto"))
    palabras <- strsplit(q, "\\s+")[[1]]
    cand <- which(Reduce(`&`, lapply(palabras, function(p) grepl(p, nombres, fixed = TRUE))))
    if (length(cand) > 0) {
      # Entre los que contienen todas las palabras, el nombre mas corto (el mas general).
      k <- cand[which.min(nchar(nombres[cand]))]
      return(list(entrada = e, hpo_id = t$id[k], como = "nombre parecido"))
    }
    list(entrada = e, hpo_id = NA_character_, como = "no encontrado")
  })
  out <- tibble::tibble(
    entrada = vapply(filas, `[[`, "", "entrada"),
    hpo_id = vapply(filas, function(f) f$hpo_id %||% NA_character_, ""),
    como = vapply(filas, `[[`, "", "como")
  )
  out$termino <- .hpo_nombre(out$hpo_id)
  out
}

# Version para las funciones: avisa lo que interpreto y devuelve los codigos.
.fenotipo_codigos <- function(fenotipo, avisar = TRUE) {
  if (is.null(fenotipo) || length(fenotipo) == 0) return(NULL)
  i <- .interpretar_fenotipo(fenotipo)
  if (any(is.na(i$hpo_id))) {
    cli::cli_abort(c("No pude interpretar {cli::qty(sum(is.na(i$hpo_id)))}{?este fenotipo/estos fenotipos}: {.val {i$entrada[is.na(i$hpo_id)]}}.",
                     "i" = "Usa codigos HPO como {.val HP:0004313} o buscalos con {.fn iei_buscar_hpo}."))
  }
  o <- .ontologia()
  no_fen <- !o$fenotipico[o$idx[i$hpo_id]]
  if (any(no_fen)) {
    cli::cli_warn(c("Se ignoran t\u00e9rminos que no describen un signo cl\u00ednico: {.val {paste(i$hpo_id[no_fen], i$termino[no_fen])}}.",
                    "i" = "Modo de herencia, curso o frecuencia no se usan para comparar fenotipos."))
    i <- i[!no_fen, ]
    if (nrow(i) == 0) return(NULL)
  }
  texto <- i[i$como %in% c("diccionario", "nombre exacto", "nombre parecido", "codigo antiguo"), ]
  if (avisar && nrow(texto) > 0) {
    cli::cli_inform(c("i" = "Fenotipo interpretado como:",
                      stats::setNames(sprintf("%s \u2192 %s %s", texto$entrada, texto$hpo_id, texto$termino),
                                      rep("*", nrow(texto)))))
  }
  unique(i$hpo_id)
}

# ---- Anotaciones por enfermedad (data-raw) ----------------------------------

# Frecuencia de un signo en una enfermedad, como proporcion. HPO la expresa
# como "n/m", porcentaje o un termino de frecuencia.
.frecuencias_hpo <- c("HP:0040280" = 1, "HP:0040281" = 0.895, "HP:0040282" = 0.545,
                      "HP:0040283" = 0.17, "HP:0040284" = 0.025, "HP:0040285" = 0)

.frecuencia_numero <- function(f) {
  out <- rep(NA_real_, length(f))
  frac <- grepl("^[0-9]+/[0-9]+$", f)
  out[frac] <- vapply(strsplit(f[frac], "/"), function(p) as.numeric(p[1]) / max(as.numeric(p[2]), 1), numeric(1))
  pct <- grepl("^[0-9.]+%$", f)
  out[pct] <- as.numeric(sub("%", "", f[pct])) / 100
  term <- f %in% names(.frecuencias_hpo)
  out[term] <- .frecuencias_hpo[f[term]]
  out
}

.herencias_hpo <- c("HP:0000006" = "AD", "HP:0000007" = "AR", "HP:0001417" = "XL",
                    "HP:0001419" = "XL", "HP:0001423" = "XL", "HP:0001427" = "MT")

# Une phenotype.hpoa con genes_to_disease.txt para los genes indicados.
.hpoa_procesar <- function(ruta_hpoa, ruta_g2d, genes) {
  g2d <- utils::read.delim(ruta_g2d, quote = "", stringsAsFactors = FALSE, colClasses = "character")
  g2d <- unique(g2d[g2d$gene_symbol %in% genes, c("gene_symbol", "disease_id")])
  h <- utils::read.delim(ruta_hpoa, quote = "", comment.char = "#", stringsAsFactors = FALSE,
                         colClasses = "character", header = FALSE)
  names(h) <- c("database_id", "disease_name", "qualifier", "hpo_id", "reference", "evidence",
                "onset", "frequency", "sex", "modifier", "aspect", "biocuration")
  h <- h[h$database_id %in% g2d$disease_id, ]
  her <- h[h$aspect == "I" & h$hpo_id %in% names(.herencias_hpo), c("database_id", "hpo_id")]
  her <- tapply(.herencias_hpo[her$hpo_id], her$database_id, function(v) paste(sort(unique(v)), collapse = "/"))
  p <- h[h$aspect == "P", ]
  p$frecuencia <- .frecuencia_numero(p$frequency)
  p$excluido <- p$qualifier == "NOT"
  # Un mismo signo puede venir de varias referencias: se queda la frecuencia mayor.
  p <- p[order(p$database_id, p$hpo_id, -ifelse(is.na(p$frecuencia), 1, p$frecuencia)), ]
  p <- p[!duplicated(p[c("database_id", "hpo_id", "excluido")]), ]
  m <- merge(g2d, p[c("database_id", "disease_name", "hpo_id", "frecuencia", "excluido")],
             by.x = "disease_id", by.y = "database_id")
  m <- m[order(m$gene_symbol, m$disease_id, m$hpo_id), ]
  # Los nombres en ASCII (Chediak, Behcet), como pide CRAN para los datos.
  nombre <- iconv(m$disease_name, "UTF-8", "ASCII//TRANSLIT", sub = "")
  nombre <- gsub("[`'^\"~]", "", nombre)
  tibble::tibble(gen = m$gene_symbol, enfermedad_id = m$disease_id, enfermedad = nombre,
                 herencia = as.character(unname(her[m$disease_id])), hpo_id = m$hpo_id,
                 frecuencia = round(m$frecuencia, 3), excluido = m$excluido)
}
