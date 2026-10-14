# Analisis de secuencias Sanger: cromatogramas .ab1, FASTA o texto.

# Codigos IUPAC de una o dos bases.
.iupac <- list(A = "A", C = "C", G = "G", T = "T",
               R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"),
               W = c("A", "T"), K = c("G", "T"), M = c("A", "C"))
.iupac_de <- function(b1, b2) {
  par <- paste(sort(c(b1, b2)), collapse = "")
  c(AG = "R", CT = "Y", CG = "S", AT = "W", GT = "K", AC = "M")[[par]]
}
.complemento <- c(A = "T", C = "G", G = "C", T = "A", R = "Y", Y = "R", S = "S", W = "W",
                  K = "M", M = "K", B = "V", V = "B", D = "H", H = "D", N = "N", "-" = "-")
.revcomp <- function(x) {
  b <- rev(strsplit(x, "")[[1]])
  paste(unname(.complemento[b]), collapse = "")
}

# ---- Lectura ---------------------------------------------------------------

# Devuelve una lista comun para cualquier entrada: bases, calidad, picos,
# trazas (NULL si no es .ab1) y nombre.
.sanger_leer <- function(x, umbral_het) {
  if (file.exists(x)) {
    ext <- tolower(tools::file_ext(x))
    if (ext %in% c("ab1", "abi", "ab!")) {
      a <- iei_leer_ab1(x)
      het <- .llamar_heterocigotos(a, umbral_het)
      a$bases <- het$bases
      # Algunos equipos antiguos guardan calidad 0 en todas las bases.
      if (all(a$calidad %in% c(0L, NA))) a$calidad <- rep(NA_integer_, length(a$calidad))
      return(c(a[c("bases", "calidad", "picos", "trazas")], list(nombre = a$archivo, razon = het$razon)))
    }
    if (ext == "scf") {
      a <- iei_leer_scf(x)
      het <- .llamar_heterocigotos(a, umbral_het)
      a$bases <- het$bases
      if (all(a$calidad %in% c(0L, NA))) a$calidad <- rep(NA_integer_, length(a$calidad))
      return(c(a[c("bases", "calidad", "picos", "trazas")], list(nombre = a$archivo, razon = het$razon)))
    }
    lineas <- readLines(x, warn = FALSE)
    cabeceras <- which(startsWith(lineas, ">"))
    if (length(cabeceras) > 1) {
      cli::cli_warn("{.file {basename(x)}} tiene {length(cabeceras)} secuencias; se usa la primera.")
      lineas <- lineas[seq_len(cabeceras[2] - 1)]
    }
    secuencia <- paste(lineas[!startsWith(lineas, ">")], collapse = "")
    nombre <- basename(x)
  } else {
    if (grepl("[/\\\\]|\\.(ab1|abi|fa|fasta|fas|fna|seq|txt)$", x, ignore.case = TRUE)) {
      cli::cli_abort(c("No encuentro el archivo {.path {x}}.",
                       "i" = "Revisa la ruta. Con {.code file.choose()} puedes buscarlo en una ventana."))
    }
    secuencia <- x
    nombre <- "secuencia"
  }
  secuencia <- toupper(gsub("[^A-Za-z]", "", secuencia))
  invalidas <- setdiff(strsplit(secuencia, "")[[1]], names(.complemento))
  if (nchar(secuencia) < 30 || length(invalidas) > 0) {
    detalle <- if (length(invalidas) > 0) {
      c("x" = "Tiene caracteres que no son bases: {.val {invalidas}}.")
    } else {
      c("x" = "Tiene {nchar(secuencia)} bases y se necesitan al menos 30.")
    }
    cli::cli_abort(c("{.val {substr(x, 1, 40)}} no es un archivo ni una secuencia de ADN valida.", detalle))
  }
  list(bases = secuencia, calidad = rep(NA_integer_, nchar(secuencia)),
       picos = NULL, trazas = NULL, nombre = nombre, razon = NULL)
}

# Senal de las cuatro trazas en cada pico (maximo en +-2 puntos).
.senal_picos <- function(a) {
  tr <- as.matrix(a$trazas)
  n <- nrow(tr)
  s <- t(vapply(a$picos, function(p) {
    apply(tr[max(1, p - 2):min(n, p + 2), , drop = FALSE], 2, max)
  }, numeric(4)))
  colnames(s) <- colnames(tr)
  s
}

# Mediana movil que ignora la posicion central.
.mediana_vecinos <- function(x, ancho) {
  n <- length(x)
  vapply(seq_len(n), function(i) {
    v <- x[setdiff(max(1, i - ancho):min(n, i + ancho), i)]
    if (length(v) == 0) NA_real_ else stats::median(v, na.rm = TRUE)
  }, numeric(1))
}

# Revisa la traza en cada pico: si la segunda senal mas alta es al menos
# `umbral` de la primera, la base se cambia por su codigo IUPAC (heterocigoto).
# Para no confundir ruido con heterocigotos se exige que los vecinos esten
# limpios (segunda senal mediana menor a la mitad del umbral) y que el pico
# principal no sea mucho mas bajo que los de alrededor.
.llamar_heterocigotos <- function(a, umbral) {
  s <- .senal_picos(a)
  bases <- strsplit(a$bases, "")[[1]]
  ordenada <- t(apply(s, 1, sort, decreasing = TRUE))
  primero <- ordenada[, 1]
  razon <- ifelse(primero > 0, ordenada[, 2] / primero, 0)
  ruido <- .mediana_vecinos(razon, 10)
  altura <- .mediana_vecinos(primero, 25)
  candidatas <- which(razon >= umbral & ruido < umbral / 2 & primero >= 0.25 * altura)
  for (i in candidatas) {
    o <- order(s[i, ], decreasing = TRUE)
    bases[i] <- .iupac_de(colnames(s)[o[1]], colnames(s)[o[2]])
  }
  list(bases = paste(bases, collapse = ""), razon = razon)
}

# Busca el punto donde la lectura pasa de limpia a picos dobles sostenidos,
# la huella de una insercion o delecion heterocigota. Devuelve el indice de la
# base o NA.
.indel_heterocigota <- function(razon, conservar, ancho = 20) {
  if (length(conservar) < 3 * ancho) return(NA_integer_)
  r <- razon[conservar]
  for (i in seq(ancho + 1, length(r) - ancho)) {
    antes <- stats::median(r[(i - ancho):(i - 1)])
    despues <- stats::median(r[i:(i + ancho - 1)])
    if (antes < 0.15 && despues >= 0.3) {
      # El cambio empieza en el primer pico doble de la ventana.
      return(conservar[i - 1L + which(r[i:(i + ancho - 1)] >= 0.3)[1]])
    }
  }
  NA_integer_
}

# Recorte de extremos de baja calidad (algoritmo de Mott modificado, como en
# phred y Biopython). Devuelve los indices que se conservan.
.recortar <- function(calidad, limite = 0.05) {
  if (all(is.na(calidad))) return(seq_along(calidad))
  q <- ifelse(is.na(calidad), 0, calidad)
  puntaje <- limite - 10^(q / -10)
  acumulado <- 0; mejor <- 0; inicio <- 1; mejor_ini <- 1; mejor_fin <- 0
  for (i in seq_along(puntaje)) {
    acumulado <- acumulado + puntaje[i]
    if (acumulado < 0) { acumulado <- 0; inicio <- i + 1 }
    if (acumulado > mejor) { mejor <- acumulado; mejor_ini <- inicio; mejor_fin <- i }
  }
  if (mejor_fin == 0) integer(0) else mejor_ini:mejor_fin
}

# ---- Referencia ------------------------------------------------------------

# Unica funcion que pide secuencia a Ensembl.
.referencia_consultar <- function(chr, inicio, fin, build) {
  httr2::request(.vep_servidor(build)) |>
    httr2::req_url_path_append("sequence", "region", "human",
                               sprintf("%s:%d..%d:1", chr, inicio, fin)) |>
    httr2::req_headers(Accept = "application/json") |>
    .agente() |>
    httr2::req_timeout(120) |>
    httr2::req_retry(max_tries = 4, is_transient = .transitorio) |>
    .ejecutar_json("Ensembl")
}

.lookup_gen <- function(gen, build) {
  httr2::request(.vep_servidor(build)) |>
    httr2::req_url_path_append("lookup", "symbol", "homo_sapiens", gen) |>
    httr2::req_headers(Accept = "application/json") |>
    .agente() |>
    httr2::req_timeout(60) |>
    httr2::req_retry(max_tries = 3, is_transient = .transitorio) |>
    .ejecutar_json("Ensembl")
}

.parsear_region <- function(region) {
  r <- gsub("[ ,]", "", region)
  p <- regmatches(r, regexec("^(?:chr)?([0-9XYMT]+):([0-9]+)-([0-9]+)$", r, perl = TRUE))[[1]]
  if (length(p) != 4) cli::cli_abort("{.arg region} debe tener la forma {.val X:101349000-101391000}.")
  list(chr = .normalizar_chr(p[2]), inicio = as.integer(p[3]), fin = as.integer(p[4]))
}

# Obtiene la secuencia genomica (hebra directa) del gen o region, con cache.
.referencia <- function(gen, region, build, cache, margen = 1000) {
  if (!is.null(region)) {
    zona <- .parsear_region(region)
    nombre <- region
  } else {
    co <- ieiprio::genes_coordenadas
    co <- co[co$gen == gen & co$build == build, ]
    if (nrow(co) > 0) {
      zona <- list(chr = co$chr[1], inicio = co$inicio[1] - margen, fin = co$fin[1] + margen)
    } else {
      l <- .lookup_gen(gen, build)
      zona <- list(chr = .normalizar_chr(l$seq_region_name), inicio = as.integer(l$start) - margen,
                   fin = as.integer(l$end) + margen)
    }
    nombre <- gen
  }
  zona$inicio <- max(1L, as.integer(zona$inicio))
  if (zona$fin - zona$inicio > 5e6) cli::cli_abort("La regi\u00f3n de referencia supera 5 Mb; indica una {.arg region} m\u00e1s peque\u00f1a.")
  clave <- sprintf("ref_%s_%s_%d_%d.rds", build, zona$chr, zona$inicio, zona$fin)
  ruta <- if (isFALSE(cache)) NULL else file.path(cache, clave)
  if (!is.null(ruta) && file.exists(ruta)) return(readRDS(ruta))
  r <- .referencia_consultar(zona$chr, zona$inicio, zona$fin, build)
  ref <- c(zona, list(secuencia = toupper(r$seq), nombre = nombre))
  if (!is.null(ruta)) {
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    saveRDS(ref, ruta)
  }
  ref
}

# Referencia dada por el usuario: ruta a un FASTA cuya cabecera trae la
# region ("X:101353540-101354439"), o una lista con `secuencia`, `inicio` y,
# si no se puede deducir del gen o la region, `chr`.
.referencia_propia <- function(referencia, gen = NULL, region = NULL, build = "GRCh38") {
  if (is.character(referencia) && length(referencia) == 1 && file.exists(referencia)) {
    lineas <- readLines(referencia, warn = FALSE)
    cab <- lineas[startsWith(lineas, ">")][1]
    zona <- regmatches(cab, regexpr("(chr)?[0-9XYMT]+:[0-9,]+-[0-9,]+", cab))
    if (length(zona) == 0) {
      cli::cli_abort(c("La cabecera de {.path {basename(referencia)}} no trae la regi\u00f3n.",
                       "i" = "Escr\u00edbela como {.val >BTK X:101353540-101354439}."))
    }
    z <- .parsear_region(zona)
    referencia <- list(secuencia = paste(lineas[!startsWith(lineas, ">")], collapse = ""),
                       chr = z$chr, inicio = z$inicio)
  }
  if (!is.list(referencia) || is.null(referencia$secuencia) || is.null(referencia$inicio)) {
    cli::cli_abort("{.arg referencia} debe ser un FASTA con la regi\u00f3n en la cabecera o una lista con {.field secuencia} e {.field inicio}.")
  }
  chr <- referencia$chr
  if (is.null(chr) && !is.null(region)) chr <- .parsear_region(region)$chr
  if (is.null(chr) && !is.null(gen)) {
    co <- ieiprio::genes_coordenadas
    chr <- co$chr[co$gen == gen & co$build == build][1]
  }
  if (is.null(chr) || is.na(chr)) cli::cli_abort("Indica el cromosoma de la referencia en {.field chr}.")
  secuencia <- toupper(gsub("[^A-Za-z]", "", referencia$secuencia))
  inicio <- as.integer(referencia$inicio)
  list(chr = .normalizar_chr(chr), inicio = inicio, fin = inicio + nchar(secuencia) - 1L,
       secuencia = secuencia, nombre = gen %||% region %||% "referencia")
}

# ---- Alineamiento ----------------------------------------------------------

# Ubica la lectura en la referencia con semillas de k-meros, probando ambas
# hebras. Devuelve la hebra y la ventana de referencia a alinear.
.ubicar <- function(lectura, ref, k = 15) {
  n_ref <- nchar(ref)
  kmeros_ref <- substring(ref, 1:(n_ref - k + 1), k:n_ref)
  probar <- function(q) {
    n <- nchar(q)
    if (n < k) return(list(n = 0))
    km <- substring(q, 1:(n - k + 1), k:n)
    ok <- !grepl("[^ACGT]", km)
    pos_ref <- match(km[ok], kmeros_ref)
    pos_q <- which(ok)
    hits <- !is.na(pos_ref)
    if (!any(hits)) return(list(n = 0))
    desfase <- pos_ref[hits] - pos_q[hits]
    moda <- as.integer(names(which.max(table(desfase))))
    list(n = sum(abs(desfase - moda) <= 30), desfase = moda)
  }
  directa <- probar(lectura)
  inversa <- probar(.revcomp(lectura))
  mejor <- if (inversa$n > directa$n) c(inversa, hebra = "-") else c(directa, hebra = "+")
  if (mejor$n < 3) return(NULL)
  ini <- max(1L, mejor$desfase - 60L)
  fin <- min(n_ref, mejor$desfase + nchar(lectura) + 60L)
  list(hebra = mejor$hebra, inicio = ini, fin = fin, semillas = mejor$n)
}

# Alineamiento semiglobal: toda la lectura se alinea, los extremos de la
# referencia son libres. Programacion dinamica con huecos lineales, vectorizada
# por fila (la propagacion horizontal se resuelve con cummax).
.alinear <- function(q, r, coincide = 2, distinta = -3, hueco = -5) {
  qv <- strsplit(q, "")[[1]]; rv <- strsplit(r, "")[[1]]
  n <- length(qv); m <- length(rv)
  compat <- outer(names(.complemento), c("A", "C", "G", "T"), function(a, b) {
    mapply(function(x, y) y %in% (.iupac[[x]] %||% character(0)), a, b)
  })
  dimnames(compat) <- list(names(.complemento), c("A", "C", "G", "T"))
  r_idx <- match(rv, c("A", "C", "G", "T"))
  puntos <- ifelse(compat, coincide, distinta)
  puntos["N", ] <- 0
  H_prev <- numeric(m + 1)                      # fila 0: extremos de la referencia libres
  ptr <- matrix(as.raw(0), n, m)                # 1 diagonal, 2 arriba, 3 izquierda
  k_g <- hueco * (0:m)
  for (i in seq_len(n)) {
    s <- if (qv[i] %in% rownames(puntos)) puntos[qv[i], r_idx] else rep(distinta, m)
    s[is.na(s)] <- distinta
    diag <- H_prev[1:m] + s
    arriba <- H_prev[2:(m + 1)] + hueco
    temp <- pmax(diag, arriba)
    t0 <- c(H_prev[1] + hueco, temp)
    H <- k_g + cummax(t0 - k_g)
    fila <- ifelse(H[-1] > temp, 3L, ifelse(temp == diag, 1L, 2L))
    ptr[i, ] <- as.raw(fila)
    H_prev <- H
  }
  j <- which.max(H_prev[-1]); i <- n
  aq <- character(0); ar <- character(0); iq <- integer(0); jr <- integer(0)
  while (i > 0) {
    p <- if (j == 0) 2L else as.integer(ptr[i, j])
    if (p == 1L) { aq <- c(qv[i], aq); ar <- c(rv[j], ar); iq <- c(i, iq); jr <- c(j, jr); i <- i - 1; j <- j - 1 }
    else if (p == 2L) { aq <- c(qv[i], aq); ar <- c("-", ar); iq <- c(i, iq); jr <- c(NA, jr); i <- i - 1 }
    else { aq <- c("-", aq); ar <- c(rv[j], ar); iq <- c(NA, iq); jr <- c(j, jr); j <- j - 1 }
  }
  tibble::tibble(q = aq, r = ar, iq = iq, jr = jr)
}

# ---- Variantes -------------------------------------------------------------

# `calidad_het` es la calidad del entorno de cada base: el programa del
# secuenciador baja la calidad de un pico doble aunque sea un heterocigoto
# real, asi que para esos se juzga la calidad de las bases vecinas.
.variantes_de_alineamiento <- function(al, ref_ini, chr, calidad, min_calidad, calidad_het = calidad) {
  vars <- list()
  agregar <- function(...) vars[[length(vars) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
  pos_gen <- function(j) ref_ini + j - 1L
  n <- nrow(al)
  k <- 1
  while (k <= n) {
    a <- al[k, ]
    if (a$q != "-" && a$r != "-") {
      q <- a$q; r <- a$r
      opciones <- .iupac[[q]]
      cal <- if (length(opciones) == 2) calidad_het[a$iq] else calidad[a$iq]
      if (!(q %in% c(r, "N")) && (is.na(cal) || cal >= min_calidad)) {
        if (length(opciones) == 2 && r %in% opciones) {
          agregar(pos = pos_gen(a$jr), ref = r, alt = setdiff(opciones, r), genotipo = "het",
                  iq = a$iq, calidad = cal)
        } else if (length(opciones) == 1) {
          agregar(pos = pos_gen(a$jr), ref = r, alt = q, genotipo = "hom_alt", iq = a$iq, calidad = cal)
        }
      }
      k <- k + 1
    } else {
      # Bloque de huecos: delecion (hueco en la lectura) o insercion (en la referencia).
      tipo <- if (a$q == "-") "del" else "ins"
      fin <- k
      while (fin + 1 <= n && ((tipo == "del" && al$q[fin + 1] == "-") || (tipo == "ins" && al$r[fin + 1] == "-"))) fin <- fin + 1
      ancla <- k - 1
      if (ancla >= 1 && fin < n && !is.na(al$jr[ancla])) {
        bloque <- al[k:fin, ]
        cal_vec <- calidad[stats::na.omit(c(al$iq[ancla], bloque$iq, al$iq[fin + 1]))]
        cal <- if (all(is.na(cal_vec))) NA_real_ else min(cal_vec, na.rm = TRUE)
        if (is.na(cal) || cal >= min_calidad) {
          base_ancla <- al$r[ancla]
          if (tipo == "del") {
            agregar(pos = pos_gen(al$jr[ancla]), ref = paste0(base_ancla, paste(bloque$r, collapse = "")),
                    alt = base_ancla, genotipo = "hom_alt", iq = al$iq[ancla], calidad = cal)
          } else {
            agregar(pos = pos_gen(al$jr[ancla]), ref = base_ancla,
                    alt = paste0(base_ancla, paste(bloque$q, collapse = "")), genotipo = "hom_alt",
                    iq = al$iq[ancla], calidad = cal)
          }
        }
      }
      k <- fin + 1
    }
  }
  out <- do.call(rbind, vars)
  if (is.null(out)) out <- data.frame(pos = integer(0), ref = character(0), alt = character(0),
                                      genotipo = character(0), iq = integer(0), calidad = numeric(0))
  out$chr <- rep(chr, nrow(out))
  out
}

# Fraccion de senal del alelo alternativo en el pico (solo .ab1, solo SNV).
.fraccion_traza <- function(lect, iq_original, ref, alt, hebra) {
  if (is.null(lect$trazas) || nchar(ref) != 1 || nchar(alt) != 1) return(NA_real_)
  if (hebra == "-") { ref <- .complemento[[ref]]; alt <- .complemento[[alt]] }
  tr <- as.matrix(lect$trazas); p <- lect$picos[iq_original]
  v <- apply(tr[max(1, p - 2):min(nrow(tr), p + 2), , drop = FALSE], 2, max)
  if (v[[ref]] + v[[alt]] == 0) return(NA_real_)
  v[[alt]] / (v[[ref]] + v[[alt]])
}

# ---- Funcion principal -----------------------------------------------------

#' Analizar secuencias Sanger
#'
#' Ubica una o varias lecturas Sanger en el genoma, las alinea contra la
#' secuencia de referencia y devuelve las diferencias como una tabla de
#' variantes lista para [iei_anotar()] e [iei_priorizar()]. Acepta
#' cromatogramas `.ab1`, archivos FASTA o la secuencia escrita directamente.
#'
#' Acepta cromatogramas `.ab1` y `.scf` (ver [iei_leer_scf()]). Con
#' cromatogramas, cada pico se revisa en las cuatro trazas: si la segunda
#' senal alcanza `umbral_het` de la principal, la posicion se considera
#' heterocigota. Los extremos de baja calidad se recortan antes de alinear.
#' La lectura puede venir en cualquier hebra; se prueban las dos.
#'
#' Si das varias lecturas del mismo amplicon (por ejemplo directa y reversa),
#' las variantes que se repiten se juntan y la columna `lecturas` dice en
#' cuantas se vieron.
#'
#' @param x Rutas a archivos `.ab1`, FASTA o texto, o una secuencia de ADN.
#' @param gen Gen secuenciado. Se usa para bajar la referencia de Ensembl.
#' @param region Alternativa a `gen`: region genomica como
#'   `"X:101349000-101391000"`.
#' @param muestra Codigo de la muestra.
#' @param sexo `"M"` o `"F"`. En un varon, una variante homocigota en el X
#'   fuera de la region pseudoautosomica pasa a ser hemicigota.
#' @param build `"GRCh38"` o `"GRCh37"`.
#' @param min_calidad Calidad Phred minima de la base para llamar una variante.
#' @param umbral_het Proporcion minima del segundo pico para llamar un
#'   heterocigoto.
#' @param recorte Error maximo por base del recorte de extremos (algoritmo de
#'   Mott). 0.01 equivale a exigir en promedio Phred 20; 0.05 recorta menos.
#' @param cebadores Opcional, secuencias de los cebadores de la PCR o de la
#'   secuenciacion (5' a 3'). Se buscan en la lectura (tolerando 2
#'   diferencias) y se recorta el cebador con lo que queda por fuera, porque
#'   esas bases copian al oligonucleotido y no al paciente.
#' @param indels_het Si `TRUE`, cuando una lectura pasa de limpia a picos
#'   dobles sostenidos se intenta llamar la insercion o delecion heterocigota
#'   probando cada desfase posible (hasta 30 bases) y eligiendo el que explica
#'   los picos dobles, como Indigo o Poly Peak Parser. Si ningun desfase
#'   explica al menos el 80% de los picos se avisa para revision manual.
#' @param cache Carpeta para guardar la secuencia de referencia, o `FALSE`.
#' @param referencia Opcional, la secuencia de referencia sin consultar
#'   Ensembl. Puede ser la ruta a un FASTA cuya cabecera trae la region
#'   (`>BTK X:101353540-101354439`) o una lista con `secuencia`, `inicio`
#'   (coordenada de la primera base) y `chr`. El paquete trae la de BTK que
#'   usan los ejemplos en
#'   `system.file("extdata", "btk_referencia.fa", package = "ieiprio")`.
#'
#' @return Una tabla como la de [iei_leer_vcf()] con las columnas extra
#'   `lecturas` y `fuente`. En `qual` va la calidad Phred de la base y en
#'   `ab` la fraccion de senal del alelo alternativo en el cromatograma. El
#'   atributo `ieiprio_sanger` guarda, por lectura, el alineamiento, la hebra,
#'   la identidad y las trazas, que usa [iei_cromatograma()].
#' @export
#' @examples
#' ab1 <- iei_ejemplo_sanger("het")
#' ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
#' v <- iei_sanger(ab1, gen = "BTK", sexo = "F", referencia = ref)
#' v
#' \donttest{
#' iei_cromatograma(v)
#' }
iei_sanger <- function(x, gen = NULL, region = NULL, muestra = "sanger", sexo = NULL,
                       build = c("GRCh38", "GRCh37"), min_calidad = 20, umbral_het = 0.33,
                       recorte = 0.01, cebadores = NULL, indels_het = TRUE,
                       cache = tools::R_user_dir("ieiprio", "cache"), referencia = NULL) {
  build <- match.arg(build)
  if (is.null(gen) && is.null(region) && is.null(referencia)) {
    cli::cli_abort("Indica el {.arg gen} secuenciado, una {.arg region} o una {.arg referencia}.")
  }
  if (!is.null(gen)) {
    g <- .normalizar_simbolo(gen, ieiprio::iuis_panel$gen, ieiprio::iuis_alias)
    if (!is.na(g)) gen <- g
  }
  # Primero se leen las lecturas: un error de ruta aparece antes de ir a Ensembl.
  crudas <- lapply(x, .sanger_leer, umbral_het = umbral_het)
  ref <- if (is.null(referencia)) .referencia(gen, region, build, cache) else
    .referencia_propia(referencia, gen, region, build)

  lecturas <- lapply(crudas, function(lect) {
    conservar <- intersect(.recortar(lect$calidad, recorte), .recortar_cebadores(lect$bases, cebadores))
    if (length(conservar) < 30) {
      cli::cli_warn("{lect$nombre}: menos de 30 bases de buena calidad; se omite.")
      return(NULL)
    }
    q <- substr(lect$bases, min(conservar), max(conservar))
    ubic <- .ubicar(q, ref$secuencia)
    if (is.null(ubic)) {
      cli::cli_warn("{lect$nombre}: no se pudo ubicar en {ref$nombre}. Revisa que sea el gen correcto.")
      return(NULL)
    }
    q_al <- if (ubic$hebra == "-") .revcomp(q) else q
    al <- .alinear(q_al, substr(ref$secuencia, ubic$inicio, ubic$fin))
    # Indice de cada base alineada en la lectura original (para trazas y calidad).
    idx <- min(conservar) - 1L + if (ubic$hebra == "-") nchar(q) - al$iq + 1L else al$iq
    cal_al <- cal_ctx <- rep(NA_real_, nchar(q_al))
    contexto <- .mediana_vecinos(lect$calidad, 5)
    cal_al[stats::na.omit(al$iq)] <- lect$calidad[idx[!is.na(al$iq)]]
    cal_ctx[stats::na.omit(al$iq)] <- contexto[idx[!is.na(al$iq)]]
    ref_ini <- ref$inicio + ubic$inicio - 1L
    vars <- .variantes_de_alineamiento(al, ref_ini, ref$chr, cal_al, min_calidad, cal_ctx)
    norm <- .normalizar_indels(vars, ref)
    vars$pos <- norm$pos; vars$ref <- norm$ref; vars$alt <- norm$alt
    vars$iq_original <- min(conservar) - 1L + if (ubic$hebra == "-") nchar(q) - vars$iq + 1L else vars$iq
    vars$ab <- vapply(seq_len(nrow(vars)), function(k) {
      .fraccion_traza(lect, vars$iq_original[k], vars$ref[k], vars$alt[k], ubic$hebra)
    }, numeric(1))
    vars$lectura <- rep(lect$nombre, nrow(vars))
    alineadas <- sum(!is.na(al$iq) & !is.na(al$jr))
    iguales <- sum(al$q == al$r)
    al$pos <- ifelse(is.na(al$jr), NA_integer_, ref_ini + al$jr - 1L)
    al$iq_original <- idx
    al$calidad <- cal_al[al$iq]
    indel <- if (is.null(lect$razon)) NA_integer_ else .indel_heterocigota(lect$razon, conservar)
    pos_indel <- if (is.na(indel)) NA_integer_ else {
      cerca <- al$pos[!is.na(al$pos) & !is.na(idx)][which.min(abs(idx[!is.na(al$pos) & !is.na(idx)] - indel))]
      as.integer(cerca)
    }
    dec <- if (!is.na(pos_indel) && isTRUE(indels_het)) .deconvolucionar_indel(lect, al, ubic$hebra, ref, indel) else NULL
    if (!is.null(dec)) {
      # Lo que se llamo despues del punto de quiebre sale de una mezcla de dos
      # alelos desfasados y no es confiable: se reemplaza por la indel.
      fuera <- !is.na(vars$iq_original) & vars$iq_original >= dec$iq
      vars <- vars[!fuera, , drop = FALSE]
      nueva <- vars[0, , drop = FALSE]
      nueva[1, c("pos", "ref", "alt", "genotipo", "chr", "lectura")] <-
        list(dec$pos, dec$ref, dec$alt, "het", ref$chr, lect$nombre)
      nueva$iq_original <- dec$iq
      vars <- rbind(vars, nueva)
      cli::cli_inform(c("v" = "{lect$nombre}: indel heterocigota {ref$chr}:{dec$pos} {dec$ref}>{dec$alt} por deconvoluci\u00f3n de picos dobles ({round(100 * dec$puntaje)}% de los picos explicados).",
                        "i" = "Confirmala con la lectura de la otra hebra o con {.fn iei_cromatograma}."))
    } else if (!is.na(pos_indel)) {
      cli::cli_warn(c("{lect$nombre}: desde {ref$chr}:{pos_indel} la lectura muestra picos dobles sostenidos.",
                      "i" = "Suele ser una inserci\u00f3n o deleci\u00f3n heterocigota que no se pudo resolver. Rev\u00edsala con {.fn iei_cromatograma}."))
    }
    list(variantes = vars,
         info = list(lectura = lect$nombre, hebra = ubic$hebra, bases = nchar(lect$bases),
                     recortadas = nchar(q), alineadas = alineadas,
                     identidad = round(100 * iguales / nrow(al), 1),
                     region = sprintf("%s:%d-%d", ref$chr, min(al$pos, na.rm = TRUE), max(al$pos, na.rm = TRUE)),
                     indel_heterocigota = pos_indel,
                     indel_llamada = if (is.null(dec)) NA_character_ else sprintf("%s:%d %s>%s", ref$chr, dec$pos, dec$ref, dec$alt),
                     alineamiento = al, picos = lect$picos, trazas = lect$trazas))
  })
  lecturas <- Filter(Negate(is.null), lecturas)
  if (length(lecturas) == 0) cli::cli_abort("Ninguna lectura pudo alinearse.")

  todas <- do.call(rbind, lapply(lecturas, `[[`, "variantes"))
  tabla <- .juntar_lecturas(todas)
  sexo_m <- identical(toupper(sexo %||% ""), "M")
  tabla$genotipo[sexo_m & tabla$chr == "X" & !.en_par(tabla$pos, build) & tabla$genotipo == "hom_alt"] <- "hemi"

  out <- tibble::tibble(
    muestra = rep(muestra, nrow(tabla)), chr = tabla$chr, pos = as.integer(tabla$pos), id = NA_character_, ref = tabla$ref, alt = tabla$alt,
    genotipo = tabla$genotipo, qual = tabla$calidad, dp = NA_real_, gq = NA_real_, ab = tabla$ab,
    filtro = "PASS", build = build, lecturas = as.integer(tabla$lecturas), fuente = rep("sanger", nrow(tabla))
  )
  out <- suppressMessages(.genes_de_variantes(out, build))
  # Si se indico el gen, las variantes de los flancos tambien se le asignan.
  if (!is.null(gen)) out$gen[is.na(out$gen)] <- gen
  info <- lapply(lecturas, `[[`, "info")
  names(info) <- vapply(info, `[[`, "", "lectura")
  attr(out, "build") <- build
  attr(out, "ieiprio_sexo") <- stats::setNames(if (is.null(sexo)) NA_character_ else toupper(sexo), muestra)
  attr(out, "ieiprio_vcf") <- list(archivo = paste(vapply(info, `[[`, "", "lectura"), collapse = ", "),
                                   md5 = .md5_texto(vapply(lecturas, function(l) paste(l$info$alineamiento$q, collapse = ""), "")))
  attr(out, "ieiprio_sanger") <- info
  class(out) <- c("iei_sanger", class(out))
  attr(out, "ieiprio_bitacora") <- tibble::tibble(
    paso = c("Lecturas alineadas", "Variantes encontradas"),
    filas = c(length(lecturas), nrow(out)), descartadas = c(length(x) - length(lecturas), 0L))
  for (l in info) {
    cli::cli_inform("{l$lectura}: hebra {l$hebra}, {l$alineadas} bases alineadas en {l$region}, identidad {l$identidad}%.")
  }
  out
}

# Normalizacion a la izquierda de inserciones y deleciones, como la espera
# VCF: el hueco que el alineamiento puso dentro de una repeticion se corre a
# su posicion mas a la izquierda.
.normalizar_indels <- function(vars, ref) {
  base_en <- function(pos) substr(ref$secuencia, pos - ref$inicio + 1, pos - ref$inicio + 1)
  for (k in seq_len(nrow(vars))) {
    r <- vars$ref[k]; a <- vars$alt[k]; p <- vars$pos[k]
    if (nchar(r) == nchar(a)) next
    repeat {
      if (nchar(r) > 0 && nchar(a) > 0 && substring(r, nchar(r)) == substring(a, nchar(a))) {
        r <- substr(r, 1, nchar(r) - 1); a <- substr(a, 1, nchar(a) - 1)
      } else if ((nchar(r) == 0 || nchar(a) == 0) && p > ref$inicio) {
        p <- p - 1L; b <- base_en(p); r <- paste0(b, r); a <- paste0(b, a)
      } else break
    }
    while (nchar(r) > 1 && nchar(a) > 1 && substr(r, 1, 1) == substr(a, 1, 1)) {
      r <- substring(r, 2); a <- substring(a, 2); p <- p + 1L
    }
    vars$ref[k] <- r; vars$alt[k] <- a; vars$pos[k] <- p
  }
  vars[c("pos", "ref", "alt")]
}

# Junta las variantes de varias lecturas: misma posicion y alelo es la misma
# variante. Si una lectura la ve heterocigota y otra no, se conserva het.
.juntar_lecturas <- function(v) {
  if (nrow(v) == 0) return(cbind(v, lecturas = integer(0)))
  clave <- paste(v$chr, v$pos, v$ref, v$alt)
  grupos <- split(seq_len(nrow(v)), clave)
  filas <- lapply(grupos, function(idx) {
    d <- v[idx, , drop = FALSE]
    mejor <- d[which.max(ifelse(is.na(d$calidad), 0, d$calidad)), , drop = FALSE]
    mejor$genotipo <- if (any(d$genotipo == "het")) "het" else mejor$genotipo
    mejor$lecturas <- length(unique(d$lectura))
    mejor
  })
  out <- do.call(rbind, filas)
  out[order(out$chr, out$pos), ]
}

# ---- Grafico ---------------------------------------------------------------

#' Graficar el cromatograma alrededor de una variante
#'
#' Dibuja las cuatro trazas de cada lectura `.ab1` en una ventana alrededor de
#' una variante encontrada por [iei_sanger()], con la base llamada arriba y la
#' base de referencia abajo. Las lecturas reversas se muestran ya
#' complementadas y en el sentido del genoma, para compararlas directamente
#' con las directas.
#'
#' @param x Resultado de [iei_sanger()], o de [iei_anotar()] o
#'   [iei_priorizar()] aplicados sobre el.
#' @param variante Fila de `x` a mostrar.
#' @param posicion Alternativa a `variante`: posicion genomica a mostrar,
#'   util para revisar una zona sin variante o la marcada como posible indel.
#' @param lectura Nombre de la lectura a mostrar. Por defecto, todas las que
#'   cubren la posicion.
#' @param ventana Bases a cada lado de la posicion.
#' @param estilo `"clasico"` (por defecto) dibuja las trazas sobre fondo
#'   blanco, como un visor de cromatogramas. `"iei"` usa el tema del resto de
#'   los graficos del paquete.
#' @return Un grafico de ggplot2.
#' @export
#' @examples
#' ab1 <- iei_ejemplo_sanger("het")
#' ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
#' v <- iei_sanger(ab1, gen = "BTK", referencia = ref)
#' \donttest{
#' iei_cromatograma(v)
#' iei_cromatograma(v, estilo = "iei")
#' }
iei_cromatograma <- function(x, variante = 1, posicion = NULL, lectura = NULL, ventana = 12,
                             estilo = c("clasico", "iei")) {
  estilo <- match.arg(estilo)
  rlang::check_installed("ggplot2", reason = "para dibujar el cromatograma.")
  info <- attr(x, "ieiprio_sanger")
  if (is.null(info)) cli::cli_abort("{.arg x} no viene de {.fn iei_sanger}.")
  if (is.null(posicion)) {
    if (nrow(x) == 0) cli::cli_abort("No hay variantes; indica una {.arg posicion}.")
    if (variante > nrow(x)) cli::cli_abort("{.arg x} solo tiene {nrow(x)} variante{?s}.")
    posicion <- x$pos[variante]
    titulo <- sprintf("%s:%d %s>%s (%s)", x$chr[variante], posicion, x$ref[variante], x$alt[variante],
                      x$genotipo[variante])
  } else {
    titulo <- sprintf("Posici\u00f3n %d", posicion)
  }
  if (!is.null(lectura)) info <- info[intersect(names(info), lectura)]
  colores <- c(A = "#1b9e3e", C = "#2166ac", G = "#222222", T = "#d7301f")

  partes <- lapply(info, function(l) {
    if (is.null(l$trazas)) return(NULL)
    al <- l$alineamiento
    fila <- which(al$pos == posicion & !is.na(al$iq_original))
    if (length(fila) == 0) return(NULL)
    fila <- fila[1]
    filas <- max(1, fila - ventana):min(nrow(al), fila + ventana)
    sub <- al[filas, ]
    # Las bases borradas no tienen pico: se intercalan entre los vecinos.
    con_pico <- !is.na(sub$iq_original)
    picos <- rep(NA_real_, nrow(sub))
    picos[con_pico] <- l$picos[sub$iq_original[con_pico]]
    if (any(!con_pico) && sum(con_pico) >= 2) {
      picos <- stats::approx(which(con_pico), picos[con_pico], xout = seq_along(picos), rule = 2)$y
    }
    sub <- sub[!is.na(picos), ]
    picos <- picos[!is.na(picos)]
    desde <- min(picos) - 6; hasta <- max(picos) + 6
    tr <- l$trazas[max(1, desde):min(nrow(l$trazas), hasta), ]
    puntos <- seq(max(1, desde), length.out = nrow(tr))
    nombres <- c("A", "C", "G", "T")
    if (l$hebra == "-") nombres <- unname(.complemento[nombres])
    # En la reversa el eje se invierte para leer en el sentido del genoma.
    eje <- if (l$hebra == "-") -puntos else puntos
    largo <- data.frame(
      x = rep(eje, 4), senal = unlist(tr, use.names = FALSE),
      base = rep(nombres, each = nrow(tr)), lectura = sprintf("%s (%s)", l$lectura, l$hebra)
    )
    centro <- l$picos[al$iq_original[fila]]
    etiquetas <- data.frame(
      x = if (l$hebra == "-") -picos else picos, llamada = sub$q, referencia = sub$r,
      distinta = sub$q != sub$r, lectura = largo$lectura[1]
    )
    list(trazas = largo, etiquetas = etiquetas,
         centro = data.frame(x = if (l$hebra == "-") -centro else centro, lectura = largo$lectura[1]))
  })
  partes <- Filter(Negate(is.null), partes)
  if (length(partes) == 0) cli::cli_abort("Ninguna lectura con cromatograma cubre la posicion {posicion}.")
  trazas <- do.call(rbind, lapply(partes, `[[`, "trazas"))
  etiquetas <- do.call(rbind, lapply(partes, `[[`, "etiquetas"))
  centro <- do.call(rbind, lapply(partes, `[[`, "centro"))
  tope <- stats::aggregate(senal ~ lectura, trazas, max)
  etiquetas$alto <- tope$senal[match(etiquetas$lectura, tope$lectura)]

  ggplot2::ggplot(trazas, ggplot2::aes(x = .data$x, y = .data$senal, colour = .data$base)) +
    ggplot2::geom_vline(data = centro, ggplot2::aes(xintercept = .data$x),
                        colour = "#f4a582", linewidth = 6, alpha = 0.35) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_text(data = etiquetas, ggplot2::aes(x = .data$x, y = .data$alto * 1.08, label = .data$llamada,
                                                       fontface = ifelse(.data$distinta, "bold", "plain")),
                       inherit.aes = FALSE, size = 3.2) +
    ggplot2::geom_text(data = etiquetas, ggplot2::aes(x = .data$x, y = -.data$alto * 0.07, label = .data$referencia),
                       inherit.aes = FALSE, size = 2.8, colour = "grey45") +
    ggplot2::scale_colour_manual(values = colores) +
    ggplot2::facet_wrap(~lectura, ncol = 1, scales = "free") +
    ggplot2::labs(title = titulo, subtitle = "Arriba la base le\u00edda, abajo la referencia",
                  x = NULL, y = "Se\u00f1al", colour = NULL) +
    (if (estilo == "clasico") ggplot2::theme_minimal(base_size = 11) else .tema_iei(base_size = 11)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                   legend.position = "top")
}
