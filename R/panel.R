.categorias_iuis <- c(
  "1"  = "Inmunodeficiencias que afectan la inmunidad celular y humoral",
  "2"  = "Inmunodeficiencias combinadas con rasgos sindromicos",
  "3"  = "Deficiencias predominantemente de anticuerpos",
  "4"  = "Enfermedades de desregulacion inmune",
  "5"  = "Defectos congenitos en numero o funcion de fagocitos",
  "6"  = "Defectos de la inmunidad intrinseca e innata",
  "7"  = "Enfermedades autoinflamatorias",
  "8"  = "Deficiencias del complemento",
  "9"  = "Falla medular",
  "10" = "Fenocopias de errores innatos de la inmunidad"
)

#' Consultar el panel de genes
#'
#' Devuelve el panel de genes de errores innatos de la inmunidad, filtrado
#' por nivel de evidencia, categoria IUIS y modo de herencia.
#'
#' @param categoria Vector de enteros entre 1 y 10, o `NULL` para todas.
#'   Solo devuelve genes con categoria IUIS curada.
#' @param herencia Vector con `"AR"`, `"AD"`, `"XL"` o `"MT"`, o `NULL` para
#'   todos. Un gen `"AR/AD"` aparece si se pide `"AR"` o `"AD"`.
#' @param evidencia Niveles de evidencia de PanelApp a incluir: `"verde"`,
#'   `"ambar"`, `"rojo"`. Por defecto solo verde, que son los genes con
#'   asociacion diagnostica establecida.
#'
#' @return Un tibble con la estructura de [iuis_panel].
#' @export
#' @examples
#' iei_panel()
#' iei_panel(categoria = 3)
#' iei_panel(herencia = "XL", evidencia = c("verde", "ambar"))
iei_panel <- function(categoria = NULL, herencia = NULL, evidencia = "verde") {
  panel <- ieiprio::iuis_panel

  evidencia <- tolower(evidencia)
  niveles <- c("verde", "ambar", "rojo")
  if (!all(evidencia %in% niveles)) {
    cli::cli_abort(c("{.arg evidencia} solo acepta {.val {niveles}}.",
                     "x" = "Recibi {.val {setdiff(evidencia, niveles)}}."))
  }
  panel <- panel[panel$evidencia %in% evidencia, ]

  if (!is.null(categoria)) {
    if (!is.numeric(categoria) || any(!categoria %in% 1:10)) {
      cli::cli_abort(c("{.arg categoria} debe ser un numero entre 1 y 10.",
                       "x" = "Recibi {.val {categoria}}."))
    }
    cats <- ieiprio::iuis_categorias
    panel <- panel[panel$gen %in% cats$gen[cats$categoria_iuis %in% categoria], ]
  }

  if (!is.null(herencia)) {
    herencia <- toupper(herencia)
    validas <- c("AR", "AD", "XL", "MT")
    invalidas <- setdiff(herencia, validas)
    if (length(invalidas) > 0) {
      cli::cli_abort(c("Modo de herencia no reconocido: {.val {invalidas}}.",
                       "i" = "Usa alguno de {.val {validas}}."))
    }
    partes <- strsplit(ifelse(is.na(panel$herencia), "", panel$herencia), "/",
                       fixed = TRUE)
    coincide <- vapply(partes, function(p) any(p %in% herencia), logical(1))
    panel <- panel[coincide, ]
  }

  tibble::as_tibble(panel)
}

#' Buscar genes en el panel
#'
#' Busca uno o mas genes por su simbolo oficial o por un alias o simbolo
#' anterior (por ejemplo `"TACI"` para TNFRSF13B). No distingue mayusculas.
#' Busca en todo el panel, sin filtrar por evidencia.
#'
#' @param genes Vector de caracteres con simbolos o alias.
#'
#' @return Un tibble con la columna `consulta` seguida de las columnas de
#'   [iuis_panel]. Los genes que no estan en el panel se omiten con un aviso.
#' @export
#' @examples
#' iei_buscar_gen("BTK")
#' iei_buscar_gen(c("taci", "STAT1"))
iei_buscar_gen <- function(genes) {
  if (!is.character(genes) || length(genes) == 0) {
    cli::cli_abort("{.arg genes} debe ser un vector de texto con al menos un gen.")
  }

  panel <- ieiprio::iuis_panel
  consulta <- trimws(genes)
  simbolo <- .normalizar_simbolo(consulta, panel$gen, ieiprio::iuis_alias)

  no_encontrados <- consulta[is.na(simbolo)]
  if (length(no_encontrados) > 0) {
    cli::cli_warn("No estan en el panel: {.val {no_encontrados}}.")
  }

  idx <- match(simbolo[!is.na(simbolo)], panel$gen)
  out <- panel[idx, ]
  tibble::add_column(tibble::as_tibble(out),
                     consulta = consulta[!is.na(simbolo)], .before = 1)
}

# Traduce simbolos o alias al simbolo oficial. NA si no lo encuentra.
# El simbolo oficial tiene prioridad sobre un alias igual de otro gen.
.normalizar_simbolo <- function(x, oficiales, alias) {
  x_may <- toupper(x)
  directo <- oficiales[match(x_may, toupper(oficiales))]
  por_alias <- alias$gen[match(x_may, toupper(alias$alias))]
  ifelse(is.na(directo), por_alias, directo)
}
