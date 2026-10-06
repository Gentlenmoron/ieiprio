# Genera los datasets `iuis_panel` e `iuis_alias`.
# Correr con source("data-raw/iuis_panel.R") desde la raiz del paquete.
#
# Por ahora usa un panel semilla de 16 filas para poder desarrollar.
# Cuando tengas la tabla completa IUIS 2024 y el archivo de HGNC,
# activa la seccion "Panel completo" de abajo.

categorias <- c(
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

semilla <- utils::read.csv("data-raw/iuis_semilla.csv", stringsAsFactors = FALSE)

iuis_panel <- tibble::tibble(
  gen              = semilla$gen,
  categoria_iuis   = as.integer(semilla$categoria_iuis),
  categoria_nombre = unname(categorias[as.character(semilla$categoria_iuis)]),
  enfermedad       = semilla$enfermedad,
  herencia         = semilla$herencia,
  ganancia_funcion = as.logical(semilla$ganancia_funcion)
)

iuis_alias <- tibble::as_tibble(
  utils::read.csv("data-raw/alias_semilla.csv", stringsAsFactors = FALSE)
)

# ---- Panel completo (activar cuando tengas los archivos) -------------------
# 1. Guarda el suplemento de la clasificacion IUIS 2024 como
#    data-raw/iuis_2024.xlsx y revisa como se llaman sus columnas.
# 2. Descarga hgnc_complete_set.txt desde genenames.org a data-raw/.
#
# hgnc <- readr::read_tsv("data-raw/hgnc_complete_set.txt")
# iuis_alias <- hgnc |>
#   dplyr::select(gen = symbol, alias_symbol, prev_symbol) |>
#   tidyr::pivot_longer(c(alias_symbol, prev_symbol), values_to = "alias") |>
#   tidyr::separate_longer_delim(alias, "|") |>
#   dplyr::filter(!is.na(alias), alias != "") |>
#   dplyr::distinct(alias, gen)
#
# crudo <- readxl::read_excel("data-raw/iuis_2024.xlsx")
# iuis_panel <- crudo |> ... (mapear columnas a la misma estructura de arriba)
# Luego reemplazar cada simbolo viejo por el oficial usando iuis_alias.
# ----------------------------------------------------------------------------

stopifnot(
  all(iuis_panel$categoria_iuis %in% 1:10),
  all(iuis_panel$herencia %in% c("AR", "AD", "XL", "XLR", "AR/AD")),
  all(iuis_alias$gen %in% iuis_panel$gen)
)

usethis::use_data(iuis_panel, iuis_alias, overwrite = TRUE)
