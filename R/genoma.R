# Datos del genoma para los graficos de cromosomas: longitud de cada
# cromosoma y posicion aproximada del centromero (solo para el dibujo).

.cromosomas <- c(as.character(1:22), "X", "Y")

.longitudes <- list(
  GRCh38 = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
             145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
             101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
             156040895, 57227415),
  GRCh37 = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
             146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
             102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
             155270560, 59373566)
)

# Centromeros en megabases, aproximados a partir de las bandas acen de UCSC.
.centromeros_mb <- c(123.4, 93.9, 90.9, 50.0, 48.8, 59.8, 60.1, 45.2, 43.0, 39.8, 53.4, 35.5,
                     17.7, 17.2, 19.0, 36.8, 25.1, 18.5, 26.2, 28.1, 12.0, 15.0, 60.6, 10.4)

.genoma <- function(build = c("GRCh38", "GRCh37")) {
  build <- match.arg(build)
  data.frame(chr = .cromosomes_factor(.cromosomas), largo = .longitudes[[build]],
             centromero = .centromeros_mb * 1e6, stringsAsFactors = FALSE)
}

.cromosomes_factor <- function(x) factor(sub("^chr", "", as.character(x)), levels = .cromosomas)

# Nombres cortos de las categorias IUIS para leyendas.
.nombres_iuis <- c(`1` = "Combinadas", `2` = "Combinadas sindr\u00f3micas", `3` = "Anticuerpos",
                   `4` = "Desregulaci\u00f3n inmune", `5` = "Fagocitos", `6` = "Inmunidad innata",
                   `7` = "Autoinflamatorias", `8` = "Complemento", `9` = "Falla medular",
                   `10` = "Fenocopias")

.colores_iuis <- c(`1` = "#2a78d6", `2` = "#6250d6", `3` = "#1baf7a", `4` = "#e34948",
                   `5` = "#eb6834", `6` = "#eda100", `7` = "#e87ba4", `8` = "#00a3c4",
                   `9` = "#8a6d3b", `10` = "#a3a29b")

# Genes del panel con su posicion y su primera categoria IUIS.
.genes_panel_posicion <- function(build = "GRCh38", genes = NULL) {
  g <- ieiprio::genes_coordenadas
  g <- g[g$build == build & g$chr %in% .cromosomas, ]
  if (!is.null(genes)) g <- g[g$gen %in% genes, ]
  cat_tab <- ieiprio::iuis_categorias
  cat_tab <- cat_tab[!duplicated(cat_tab$gen), ]
  g$categoria <- as.character(cat_tab$categoria_iuis[match(g$gen, cat_tab$gen)])
  g$color <- unname(.colores_iuis[g$categoria])
  g$color[is.na(g$color)] <- "#c9c8c2"
  g$medio <- (g$inicio + g$fin) / 2
  g
}

# Posiciones de lo que se quiere marcar en el genoma, segun el objeto.
.marcas_genoma <- function(x, build) {
  if (inherits(x, "iei_resultado")) {
    d <- as.data.frame(x)
    d <- d[d$chr %in% .cromosomas & !is.na(d$pos), ]
    if (nrow(d) == 0) cli::cli_abort("El resultado no tiene variantes con posici\u00f3n en cromosomas.")
    d <- d[order(-d$puntaje), ]
    return(data.frame(gen = d$gen, chr = d$chr, pos = d$pos, valor = d$puntaje,
                      grupo = as.character(d$categoria), etiqueta = .etiqueta_variante(d),
                      stringsAsFactors = FALSE))
  }
  if (inherits(x, "iei_busqueda")) {
    d <- tibble::as_tibble(unclass(x))
    g <- .genes_panel_posicion(build, d$gen)
    k <- match(d$gen, g$gen)
    d <- d[!is.na(k), ]; k <- k[!is.na(k)]
    return(data.frame(gen = d$gen, chr = g$chr[k], pos = g$medio[k], valor = d$similitud,
                      grupo = ifelse(is.na(d$herencia), "sin dato", d$herencia), etiqueta = d$gen,
                      stringsAsFactors = FALSE))
  }
  NULL
}
