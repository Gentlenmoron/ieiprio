# Construye todos los datasets del paquete a partir de fuentes con API.
# Correr desde la raiz del paquete:  source("data-raw/construir_datos.R")
#
# Fuentes:
#   - PanelApp (Genomics England), panel y version fijados en config.R
#   - HPO genes_to_phenotype.txt, release fijado en config.R
#   - HPO hp.obo (ontologia), mismo release
#   - data-raw/iuis_curado.csv, categorias IUIS agregadas a mano (en git)
#   - data-raw/iuis_mecanismo.csv, genes con enfermedad por ganancia de funcion (en git)
#
# Las categorias IUIS salen sobre todo de los fenotipos de PanelApp, que
# nombran la tabla IUIS de cada gen; ver .categorias_iuis_derivar().

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
mecanismo <- utils::read.csv("data-raw/iuis_mecanismo.csv", stringsAsFactors = FALSE)
faltan <- setdiff(curado$gen, pa$panel$gen)
if (length(faltan) > 0) {
  warning("Genes curados que no estan en PanelApp: ", paste(faltan, collapse = ", "))
}

iuis_panel      <- .construir_panel(pa, curado, mecanismo)
iuis_alias      <- pa$alias
genes_coordenadas <- pa$coordenadas
sin_coord <- setdiff(iuis_panel$gen, genes_coordenadas$gen[genes_coordenadas$build == "GRCh38"])
if (length(sin_coord) > 0) {
  message(length(sin_coord), " genes sin coordenadas GRCh38: ",
          paste(utils::head(sin_coord, 10), collapse = ", "),
          if (length(sin_coord) > 10) "...")
}
iuis_categorias <- .categorias_iuis_derivar(iuis_panel$gen, iuis_panel$fenotipos, curado)

# ---- HPO --------------------------------------------------------------------
hpo_cfg <- config$hpo
ruta_hpo <- file.path(dir_crudos, sprintf("genes_to_phenotype_%s.txt", hpo_cfg$tag))
if (!file.exists(ruta_hpo)) .hpo_descargar(hpo_cfg$tag, destino = ruta_hpo)
md5_hpo <- verificar_md5(ruta_hpo, hpo_cfg$md5, "HPO")
hpo_genes <- .hpo_procesar(ruta_hpo, iuis_panel$gen)

ruta_obo <- file.path(dir_crudos, sprintf("hp_%s.obo", hpo_cfg$tag))
if (!file.exists(ruta_obo)) .hpo_descargar(hpo_cfg$tag, destino = ruta_obo, archivo = "hp.obo")
md5_obo <- verificar_md5(ruta_obo, hpo_cfg$md5_obo, "HPO ontologia")
hpo_terminos <- .obo_procesar(ruta_obo)
hpo_terminos$ic <- .ic_calcular(hpo_terminos, ruta_hpo)

# Anotaciones por enfermedad: frecuencia de cada signo, signos excluidos (NOT)
# y modo de herencia de cada enfermedad.
ruta_hpoa <- file.path(dir_crudos, sprintf("phenotype_%s.hpoa", hpo_cfg$tag))
if (!file.exists(ruta_hpoa)) .hpo_descargar(hpo_cfg$tag, destino = ruta_hpoa, archivo = "phenotype.hpoa")
md5_hpoa <- verificar_md5(ruta_hpoa, hpo_cfg$md5_hpoa, "HPO enfermedades")
ruta_g2d <- file.path(dir_crudos, sprintf("genes_to_disease_%s.txt", hpo_cfg$tag))
if (!file.exists(ruta_g2d)) .hpo_descargar(hpo_cfg$tag, destino = ruta_g2d, archivo = "genes_to_disease.txt")
md5_g2d <- verificar_md5(ruta_g2d, hpo_cfg$md5_g2d, "HPO genes a enfermedad")
hpo_enfermedades <- .hpoa_procesar(ruta_hpoa, ruta_g2d, iuis_panel$gen)

# ---- Genes con dificultades tecnicas -----------------------------------------
genes_dificiles <- tibble::as_tibble(utils::read.csv("data-raw/genes_dificiles.csv",
                                                     stringsAsFactors = FALSE, encoding = "UTF-8"))

# ---- Registro de fuentes ----------------------------------------------------
md5_curado <- unname(tools::md5sum("data-raw/iuis_curado.csv"))
md5_mec <- unname(tools::md5sum("data-raw/iuis_mecanismo.csv"))
md5_dif <- unname(tools::md5sum("data-raw/genes_dificiles.csv"))
repo <- "https://github.com/Gentlenmoron/ieiprio/blob/master/"
fuentes_datos <- tibble::tibble(
  fuente  = c("PanelApp", "HPO", "HPO ontologia", "HPO enfermedades", "HPO genes a enfermedad",
              "Curacion IUIS", "Mecanismo GOF", "Genes dificiles"),
  detalle = c(pa$nombre, "genes_to_phenotype.txt", "hp.obo", "phenotype.hpoa", "genes_to_disease.txt",
              "data-raw/iuis_curado.csv", "data-raw/iuis_mecanismo.csv", "data-raw/genes_dificiles.csv"),
  version = c(pa_cfg$version, rep(hpo_cfg$tag, 4), "git", "git", "git"),
  url     = c(.panelapp_url_panel(pa_cfg$panel_id, pa_cfg$version),
              .hpo_url(hpo_cfg$tag), .hpo_url(hpo_cfg$tag, "hp.obo"),
              .hpo_url(hpo_cfg$tag, "phenotype.hpoa"), .hpo_url(hpo_cfg$tag, "genes_to_disease.txt"),
              paste0(repo, "data-raw/iuis_curado.csv"),
              paste0(repo, "data-raw/iuis_mecanismo.csv"),
              paste0(repo, "data-raw/genes_dificiles.csv")),
  fecha   = hoy,
  md5     = c(md5_pa, md5_hpo, md5_obo, md5_hpoa, md5_g2d, md5_curado, md5_mec, md5_dif)
)

usethis::use_data(iuis_panel, iuis_alias, iuis_categorias, hpo_genes, hpo_terminos,
                  hpo_enfermedades, genes_dificiles, genes_coordenadas, fuentes_datos, overwrite = TRUE, compress = "xz")

message("Listo. ", nrow(iuis_panel), " genes, ", nrow(iuis_alias), " alias, ",
        nrow(hpo_genes), " pares gen fenotipo, ",
        sum(genes_coordenadas$build == "GRCh38"), " genes con coordenadas GRCh38 y ",
        sum(genes_coordenadas$build == "GRCh37"), " con GRCh37.")
