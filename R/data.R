#' Panel de genes de errores innatos de la inmunidad (IUIS)
#'
#' Genes asociados a errores innatos de la inmunidad segun la clasificacion
#' de la International Union of Immunological Societies (IUIS). La version
#' actual es un panel semilla para desarrollo; se reemplazara por la
#' clasificacion IUIS 2024 completa.
#'
#' Un gen puede aparecer en mas de una fila si causa enfermedades distintas
#' segun el modo de herencia o el mecanismo (por ejemplo STAT1).
#'
#' @format Un tibble con una fila por combinacion de gen y enfermedad:
#' \describe{
#'   \item{gen}{Simbolo oficial HGNC.}
#'   \item{categoria_iuis}{Categoria IUIS, entero de 1 a 10.}
#'   \item{categoria_nombre}{Nombre de la categoria.}
#'   \item{enfermedad}{Enfermedad asociada.}
#'   \item{herencia}{Modo de herencia: AR, AD, XL, XLR o AR/AD.}
#'   \item{ganancia_funcion}{`TRUE` si el mecanismo es ganancia de funcion.}
#' }
#' @source Clasificacion IUIS de errores innatos de la inmunidad,
#'   Journal of Clinical Immunology.
"iuis_panel"

#' Alias de simbolos de genes del panel
#'
#' Tabla para traducir nombres antiguos o alternativos al simbolo oficial
#' HGNC usado en [iuis_panel].
#'
#' @format Un tibble con dos columnas:
#' \describe{
#'   \item{alias}{Nombre alternativo o anterior.}
#'   \item{gen}{Simbolo oficial HGNC.}
#' }
#' @source HUGO Gene Nomenclature Committee, <https://www.genenames.org>.
"iuis_alias"
