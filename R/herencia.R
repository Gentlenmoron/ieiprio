# Traduce el texto de modo de herencia de PanelApp a una etiqueta corta.
# MONOALLELIC -> AD, BIALLELIC -> AR, BOTH -> AR/AD, X-LINKED -> XL,
# MITOCHONDRIAL -> MT. Cualquier otro valor (Unknown, Other, vacio) -> NA.
.mapear_herencia <- function(x) {
  x2 <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x2))
  out[which(startsWith(x2, "MONOALLELIC"))]   <- "AD"
  out[which(startsWith(x2, "BIALLELIC"))]     <- "AR"
  out[which(startsWith(x2, "BOTH"))]          <- "AR/AD"
  out[which(startsWith(x2, "X-LINKED"))]      <- "XL"
  out[which(startsWith(x2, "MITOCHONDRIAL"))] <- "MT"
  out
}

# Nivel de confianza de PanelApp: 3 verde, 2 ambar, 0 o 1 rojo.
.mapear_evidencia <- function(x) {
  x <- suppressWarnings(as.integer(x))
  out <- rep(NA_character_, length(x))
  out[which(x == 3)] <- "verde"
  out[which(x == 2)] <- "ambar"
  out[which(x %in% c(0, 1))] <- "rojo"
  out
}
