#' Analizar un VCF de principio a fin
#'
#' Corre todo el flujo en una sola llamada. Lee el VCF, filtra por genes del
#' panel y calidad, anota con Ensembl VEP, filtra por frecuencia, prioriza y
#' genera un reporte HTML por muestra.
#'
#' @param vcf Ruta al archivo `.vcf` o `.vcf.gz`.
#' @param sexo Sexo de las muestras, ver [iei_leer_vcf()].
#' @param fenotipo Terminos HPO del paciente, ver [iei_priorizar()].
#' @param categoria Categorias IUIS sospechadas, ver [iei_priorizar()].
#' @param muestras Muestras a analizar. `NULL` analiza todas.
#' @param evidencia Niveles de evidencia de PanelApp, ver [iei_filtrar()].
#' @param carpeta Carpeta donde guardar los reportes. `NULL` no genera
#'   reportes.
#' @param responsable Nombre que aparece en los reportes.
#' @param ... Otros argumentos para [iei_filtrar()], por ejemplo `min_dp`.
#'
#' @return La tabla priorizada, como [iei_priorizar()]. Si se generaron
#'   reportes, sus rutas quedan en el atributo `reportes`.
#' @export
#' @examples
#' \dontrun{
#' res <- iei_analizar(iei_ejemplo_vcf(), sexo = "M",
#'                     fenotipo = c("HP:0004313", "HP:0002719"),
#'                     carpeta = "reportes")
#' }
iei_analizar <- function(vcf, sexo = NULL, fenotipo = NULL, categoria = NULL,
                         muestras = NULL, evidencia = "verde", carpeta = ".",
                         responsable = NULL, ...) {
  cli::cli_h1("ieiprio")
  cli::cli_alert_info("Leyendo {.file {basename(vcf)}}")
  v <- iei_leer_vcf(vcf, muestras = muestras, sexo = sexo)
  cli::cli_alert_info("Filtrando {nrow(v)} filas")
  f <- iei_filtrar(v, evidencia = evidencia, ...)
  cli::cli_alert_info("Anotando {nrow(f)} variante{?s} con Ensembl VEP")
  a <- iei_filtrar_frecuencia(iei_anotar(f))
  p <- iei_priorizar(a, fenotipo = fenotipo, categoria = categoria)

  n <- table(factor(p$categoria, levels = c("alta", "media", "baja")))
  cli::cli_alert_success("Prioridad alta {n[['alta']]}, media {n[['media']]}, baja {n[['baja']]}")

  if (!is.null(carpeta)) {
    rutas <- vapply(unique(v$muestra), function(m) {
      ruta <- file.path(carpeta, paste0(m, ".html"))
      suppressMessages(iei_reporte(p, ruta, paciente = m, responsable = responsable))
      ruta
    }, character(1))
    cli::cli_alert_success("Reporte{?s} en {.path {rutas}}")
    attr(p, "reportes") <- unname(rutas)
  }
  p
}
