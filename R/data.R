#' Panel de genes de errores innatos de la inmunidad
#'
#' Genes del panel 398 de PanelApp (Genomics England), "Primary
#' immunodeficiency or monogenic inflammatory bowel disease", en la version
#' indicada en [iei_fuentes()], unidos a una curacion propia de categorias
#' IUIS. Se genera con `data-raw/construir_datos.R`.
#'
#' @format Un tibble con una fila por gen:
#' \describe{
#'   \item{gen}{Simbolo HGNC.}
#'   \item{hgnc_id}{Identificador HGNC.}
#'   \item{herencia}{AD, AR, AR/AD, XL, MT o NA.}
#'   \item{herencia_panelapp}{Texto original de PanelApp.}
#'   \item{evidencia}{verde, ambar o rojo segun PanelApp.}
#'   \item{fenotipos}{Fenotipos listados en PanelApp, separados por punto y coma.}
#'   \item{categorias_iuis}{Categorias IUIS curadas, separadas por punto y coma. NA si no estan curadas.}
#'   \item{ganancia_funcion}{`TRUE` si el gen tiene una enfermedad por ganancia de funcion.}
#' }
#' @source <https://panelapp.genomicsengland.co.uk/panels/398/>
"iuis_panel"

#' Alias y simbolos anteriores de los genes del panel
#'
#' @format Un tibble con columnas `alias` y `gen` (simbolo oficial).
#' @source PanelApp, que a su vez toma los datos de HGNC.
"iuis_alias"

#' Categorias IUIS de los genes del panel
#'
#' Curacion manual de la clasificacion IUIS de errores innatos de la
#' inmunidad. Un gen puede estar en mas de una categoria.
#'
#' @format Un tibble con columnas `gen`, `categoria_iuis` (1 a 10) y
#'   `categoria_nombre`.
#' @source Clasificacion IUIS, Journal of Clinical Immunology.
"iuis_categorias"

#' Fenotipos HPO de los genes del panel
#'
#' @format Un tibble con columnas `gen`, `hpo_id` y `hpo_nombre`.
#' @source Human Phenotype Ontology, archivo genes_to_phenotype.txt,
#'   <https://hpo.jax.org>.
"hpo_genes"

#' Registro de fuentes de los datos
#'
#' @format Un tibble con columnas `fuente`, `detalle`, `version`, `url`,
#'   `fecha` y `md5`.
"fuentes_datos"
