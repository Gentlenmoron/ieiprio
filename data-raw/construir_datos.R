# Construye todos los datasets del paquete a partir de fuentes con API.
# Correr desde la raiz del paquete:  source("data-raw/construir_datos.R")
#
# Fuentes:
#   - PanelApp (Genomics England), panel y version fijados en config.R
#   - HPO genes_to_phenotype.txt, release fijado en config.R
#   - data-raw/iuis_curado.csv, categorias IUIS curadas a mano (en git)

devtools::load_all()
source("data-raw/config.R")

dir_crudos <- "data-raw/crudos"
dir.create(dir_crudos, showWarnings = FALSE)
hoy <- format(Sys.Date())

verificar_md5 <- function(ruta, esperado, nombre) {
  real <- unname(tools::md5sum(ruta))
  if (!is.null(esperado) && !identical(real, esperado)) {
    stop(nombre, ": la huella md5 no coincide. Esperado ", esperado,
         ", obtenido ", real, ". El archivo de origen cambio.")
  }
  if (is.null(esperado)) message(nombre, " md5 = ", real, "  (copiar en config.R)")
  real
}

# ---- PanelApp ---------------------------------------------------------------
pa_cfg <- config$panelapp
ruta_pa <- file.path(dir_crudos, sprintf("panelapp_%s_v%s.json",
                                         pa_cfg$panel_id, pa_cfg$version))
if (!file.exists(ruta_pa)) {
  .panelapp_descargar(pa_cfg$panel_id, pa_cfg$version, destino = ruta_pa)
}
md5_pa <- verificar_md5(ruta_pa, pa_cfg$md5, "PanelApp")
pa <- .panelapp_procesar(jsonlite::read_json(ruta_pa))
stopifnot(identical(as.character(pa$version), pa_cfg$version))

# ---- Curacion IUIS ----------------------------------------------------------
curado <- utils::read.csv("data-raw/iuis_curado.csv", stringsAsFactors = FALSE)
faltan <- setdiff(curado$gen, pa$panel$gen)
if (length(faltan) > 0) {
  warning("Genes curados que no estan en PanelApp: ", paste(faltan, collapse = ", "))
}

iuis_panel      <- .construir_panel(pa, curado)
iuis_alias      <- pa$alias
genes_coordenadas <- pa$coordenadas
sin_coord <- setdiff(iuis_panel$gen, genes_coordenadas$gen[genes_coordenadas$build == "GRCh38"])
if (length(sin_coord) > 0) {
  message(length(sin_coord), " genes sin coordenadas GRCh38: ",
          paste(utils::head(sin_coord, 10), collapse = ", "),
          if (length(sin_coord) > 10) "...")
}
iuis_categorias <- .construir_categorias(curado, iuis_panel$gen)

# ---- HPO --------------------------------------------------------------------
hpo_cfg <- config$hpo
ruta_hpo <- file.path(dir_crudos, sprintf("genes_to_phenotype_%s.txt", hpo_cfg$tag))
if (!file.exists(ruta_hpo)) .hpo_descargar(hpo_cfg$tag, destino = ruta_hpo)
md5_hpo <- verificar_md5(ruta_hpo, hpo_cfg$md5, "HPO")
hpo_genes <- .hpo_procesar(ruta_hpo, iuis_panel$gen)

# ---- Registro de fuentes ----------------------------------------------------
fuentes_datos <- tibble::tibble(
  fuente  = c("PanelApp", "HPO", "Curacion IUIS"),
  detalle = c(pa$nombre, "genes_to_phenotype.txt", "data-raw/iuis_curado.csv"),
  version = c(pa_cfg$version, hpo_cfg$tag, "git"),
  url     = c(.panelapp_url_panel(pa_cfg$panel_id, pa_cfg$version),
              .hpo_url(hpo_cfg$tag),
              "https://github.com/Gentlenmoron/ieiprio/blob/main/data-raw/iuis_curado.csv"),
  fecha   = hoy,
  md5     = c(md5_pa, md5_hpo, unname(tools::md5sum("data-raw/iuis_curado.csv")))
)

usethis::use_data(iuis_panel, iuis_alias, iuis_categorias, hpo_genes,
                  genes_coordenadas, fuentes_datos, overwrite = TRUE)

message("Listo. ", nrow(iuis_panel), " genes, ", nrow(iuis_alias), " alias, ",
        nrow(hpo_genes), " pares gen fenotipo, ",
        sum(genes_coordenadas$build == "GRCh38"), " genes con coordenadas GRCh38 y ",
        sum(genes_coordenadas$build == "GRCh37"), " con GRCh37.")
