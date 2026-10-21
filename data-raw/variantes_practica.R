# Banco corto de variantes reales para practicar con el paquete. Cada fila es
# una variante descrita en ClinVar y en la literatura de su enfermedad, con su
# transcrito MANE Select, su identificador de ClinVar (VCV), un genotipo y tres
# signos HPO tomados de las anotaciones de esa enfermedad en hpo_enfermedades.
# Dos filas son contraejemplos a proposito. RAG1 aparece en heterocigosis en
# un gen recesivo (portadora) y TNFRSF13B c.310T>C es un factor de riesgo con
# clasificaciones en conflicto. Se curo a mano en
# data-raw/variantes_practica.csv; los transcritos y los VCV se verificaron
# contra ClinVar. Este script valida y guarda.
devtools::load_all()

v <- utils::read.csv("data-raw/variantes_practica.csv", encoding = "UTF-8", stringsAsFactors = FALSE,
                     na.strings = "")
v <- tibble::as_tibble(v)

stopifnot(all(v$gen %in% iuis_panel$gen))
stopifnot(all(v$genotipo %in% c("het", "hom_alt", "hemi")))
stopifnot(all(grepl("^NM_[0-9]+\\.[0-9]+$", v$transcrito)))
stopifnot(all(grepl("^VCV[0-9]{9}$", v$clinvar_id)), anyDuplicated(v$clinvar_id) == 0)
for (i in seq_len(nrow(v))) {
  anot <- hpo_enfermedades[hpo_enfermedades$gen == v$gen[i] &
                             hpo_enfermedades$enfermedad_id == v$enfermedad_id[i] &
                             !hpo_enfermedades$excluido, ]
  if (nrow(anot) == 0) stop(v$gen[i], ": ", v$enfermedad_id[i], " no esta en hpo_enfermedades")
  faltan <- setdiff(strsplit(v$fenotipo[i], ";")[[1]], anot$hpo_id)
  if (length(faltan)) stop(v$gen[i], ": signos sin anotar en ", v$enfermedad_id[i], ": ", toString(faltan))
}

variantes_practica <- v
usethis::use_data(variantes_practica, overwrite = TRUE, compress = "xz")
