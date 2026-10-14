#' Aplicacion interactiva
#'
#' Abre en el navegador una aplicacion Shiny para usar ieiprio sin escribir
#' codigo: buscar genes candidatos por fenotipo, revisar la ficha de un gen,
#' interpretar una variante escrita, analizar un archivo (VCF, cromatograma
#' Sanger, informe en PDF o Word) con descarga del reporte y de la tabla en
#' Excel, y validar la busqueda por fenotipo con casos simulados.
#'
#' La aplicacion corre en tu computadora. Como el resto del paquete, envia a
#' Ensembl, NCBI y gnomAD solo las coordenadas de las variantes; puedes
#' activar el modo sin conexion desde la propia aplicacion.
#'
#' @param ... Argumentos para [shiny::runApp()], por ejemplo `port` o
#'   `launch.browser`.
#' @return No devuelve nada; corre hasta que cierres la aplicacion.
#' @export
#' @examples
#' \dontrun{
#' iei_app()
#' }
iei_app <- function(...) {
  rlang::check_installed(c("shiny", "bslib", "ggplot2"), reason = "para la aplicacion interactiva.")
  carpeta <- system.file("app", package = "ieiprio")
  if (!nzchar(carpeta)) cli::cli_abort("No encuentro la aplicacion dentro del paquete; reinstala ieiprio.")
  shiny::runApp(carpeta, ...)
}
