.herencias_validas <- c("AR", "AD", "XL", "XLR", "AR/AD")

#' Consultar el panel de genes IUIS
#'
#' Devuelve el panel completo o filtrado por categoria IUIS y modo de
#' herencia.
#'
#' @param categoria Vector de enteros entre 1 y 10, o `NULL` para todas.
#' @param herencia Vector de modos de herencia (`"AR"`, `"AD"`, `"XL"`,
#'   `"XLR"`). Un gen marcado como `"AR/AD"` aparece si se pide `"AR"` o
#'   `"AD"`. `NULL` para todos.
#'
#' @return Un tibble con la estructura de [iuis_panel].
#' @export
#' @examples
#' iei_panel()
#' iei_panel(categoria = 3)
#' iei_panel(herencia = "XL")
iei_panel <- function(categoria = NULL, herencia = NULL) {
  panel <- ieiprio::iuis_panel

  if (!is.null(categoria)) {
    if (!is.numeric(categoria) || any(!categoria %in% 1:10)) {
      cli::cli_abort(c(
        "{.arg categoria} debe ser un numero entre 1 y 10.",
        "x" = "Recibi {.val {categoria}}."
      ))
    }
    panel <- panel[panel$categoria_iuis %in% categoria, ]
  }

  if (!is.null(herencia)) {
    herencia <- toupper(herencia)
    validas <- .herencias_validas
    invalidas <- setdiff(herencia, validas)
    if (length(invalidas) > 0) {
      cli::cli_abort(c(
        "Modo de herencia no reconocido: {.val {invalidas}}.",
        "i" = "Usa alguno de {.val {validas}}."
      ))
    }
    partes <- strsplit(panel$herencia, "/", fixed = TRUE)
    coincide <- vapply(partes, function(p) any(p %in% herencia), logical(1))
    panel <- panel[coincide | panel$herencia %in% herencia, ]
  }

  tibble::as_tibble(panel)
}

#' Buscar genes en el panel IUIS
#'
#' Busca uno o mas genes por su simbolo oficial o por un alias conocido
#' (por ejemplo `"TACI"` para TNFRSF13B). La busqueda no distingue
#' mayusculas de minusculas.
#'
#' @param genes Vector de caracteres con simbolos o alias.
#'
#' @return Un tibble con la columna `consulta` (lo que se busco) seguida de
#'   las columnas de [iuis_panel]. Un gen con varias enfermedades devuelve
#'   varias filas. Los genes que no estan en el panel se omiten con un aviso.
#' @export
#' @examples
#' iei_buscar_gen("BTK")
#' iei_buscar_gen(c("taci", "STAT1"))
iei_buscar_gen <- function(genes) {
  if (!is.character(genes) || length(genes) == 0) {
    cli::cli_abort("{.arg genes} debe ser un vector de texto con al menos un gen.")
  }

  panel <- ieiprio::iuis_panel
  alias <- ieiprio::iuis_alias

  consulta <- trimws(genes)
  simbolo <- .normalizar_simbolo(consulta, panel$gen, alias)

  no_encontrados <- consulta[is.na(simbolo)]
  if (length(no_encontrados) > 0) {
    cli::cli_warn("No estan en el panel: {.val {no_encontrados}}.")
  }

  ok <- !is.na(simbolo)
  filas <- lapply(which(ok), function(i) {
    sub <- panel[panel$gen == simbolo[i], ]
    tibble::add_column(sub, consulta = consulta[i], .before = 1)
  })

  if (length(filas) == 0) {
    return(tibble::add_column(panel[0, ], consulta = character(0), .before = 1))
  }
  tibble::as_tibble(do.call(rbind, filas))
}

# Traduce simbolos o alias al simbolo oficial. Devuelve NA si no lo encuentra.
.normalizar_simbolo <- function(x, oficiales, alias) {
  x_may <- toupper(x)
  directo <- unique(oficiales)[match(x_may, toupper(unique(oficiales)))]
  por_alias <- alias$gen[match(x_may, toupper(alias$alias))]
  ifelse(is.na(directo), por_alias, directo)
}
