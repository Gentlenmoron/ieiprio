#' Fuentes y versiones de los datos del paquete
#'
#' Muestra de donde salio cada dataset incluido, con su version, URL,
#' fecha de descarga y huella md5. Sirve para citar y reproducir un analisis.
#'
#' @return Un tibble, ver [fuentes_datos].
#' @export
#' @examples
#' iei_fuentes()
iei_fuentes <- function() {
  ieiprio::fuentes_datos
}
