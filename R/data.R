#' Panel de genes de errores innatos de la inmunidad
#'
#' Genes del panel 398 de PanelApp (Genomics England), "Primary
#' immunodeficiency or monogenic inflammatory bowel disease", en la version
#' indicada en [iei_fuentes()], con su categoria IUIS y el mecanismo de
#' enfermedad. Se genera con `data-raw/construir_datos.R`.
#'
#' @format Un tibble con una fila por gen:
#' \describe{
#'   \item{gen}{Simbolo HGNC.}
#'   \item{hgnc_id}{Identificador HGNC.}
#'   \item{herencia}{AD, AR, AR/AD, XL, MT o NA.}
#'   \item{herencia_panelapp}{Texto original de PanelApp.}
#'   \item{evidencia}{verde, ambar o rojo segun PanelApp.}
#'   \item{fenotipos}{Fenotipos listados en PanelApp, separados por punto y coma.}
#'   \item{categorias_iuis}{Categorias IUIS separadas por punto y coma, ver [iuis_categorias]. NA si no tiene.}
#'   \item{mecanismo}{`"GOF"` si su enfermedad es por ganancia de funcion, `"LOF y GOF"` si tiene enfermedades por los dos mecanismos, NA en el resto (perdida de funcion o sin especificar).}
#'   \item{ganancia_funcion}{`TRUE` si `mecanismo` incluye ganancia de funcion.}
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
#' Categoria de la clasificacion IUIS de errores innatos de la inmunidad de
#' cada gen. Un gen puede estar en mas de una. La mayoria sale de los
#' fenotipos de PanelApp, que nombran la tabla IUIS del gen; cuando no la
#' nombran, de palabras clave de esos fenotipos (por ejemplo SCID o
#' agammaglobulinemia); y unas pocas de una curacion manual. La columna
#' `fuente` dice de donde salio cada una.
#'
#' @format Un tibble con columnas `gen`, `categoria_iuis` (1 a 10),
#'   `categoria_nombre` y `fuente` (`"PanelApp, tabla IUIS"`,
#'   `"PanelApp, fenotipo"` o `"curado"`).
#' @source PanelApp panel 398 y clasificacion IUIS (Poli y colaboradores).
"iuis_categorias"

#' Fenotipos HPO de los genes del panel
#'
#' @format Un tibble con columnas `gen`, `hpo_id` y `hpo_nombre`.
#' @source Human Phenotype Ontology, archivo genes_to_phenotype.txt,
#'   <https://hpo.jax.org>.
"hpo_genes"

#' Ontologia HPO
#'
#' Todos los terminos vigentes de la Human Phenotype Ontology con su
#' jerarquia, para la similitud semantica de [iei_priorizar()] y
#' [iei_buscar_fenotipo()] y la busqueda por texto de [iei_buscar_hpo()].
#'
#' @format Un tibble con una fila por termino:
#' \describe{
#'   \item{id}{Codigo HPO.}
#'   \item{nombre}{Nombre en ingles.}
#'   \item{padres}{Lista con los codigos de sus padres (relacion is_a).}
#'   \item{sinonimos}{Lista con sus sinonimos exactos.}
#'   \item{alternos}{Lista con codigos antiguos que hoy apuntan a este termino.}
#'   \item{ic}{Contenido de informacion: menos logaritmo de la fraccion de genes
#'     anotados con el termino o uno mas especifico. Mas alto es mas especifico.}
#' }
#' @source Archivos hp.obo y genes_to_phenotype.txt del release indicado en
#'   [iei_fuentes()], <https://github.com/obophenotype/human-phenotype-ontology>.
"hpo_terminos"

#' Registro de fuentes de los datos
#'
#' @format Un tibble con columnas `fuente`, `detalle`, `version`, `url`,
#'   `fecha` y `md5`.
"fuentes_datos"

#' Coordenadas de los genes del panel
#'
#' Ubicacion de cada gen en GRCh37 (Ensembl 82) y GRCh38 (Ensembl 90) segun
#' PanelApp. La usa [iei_filtrar()] para quedarse con variantes dentro de
#' genes del panel.
#'
#' @format Un tibble con columnas `gen`, `build`, `chr` (sin prefijo chr),
#'   `inicio`, `fin` y `ensembl_id`.
#' @source <https://panelapp.genomicsengland.co.uk/panels/398/>
"genes_coordenadas"

#' Genes dificiles de analizar con lecturas cortas
#'
#' Genes del panel en los que un VCF de lecturas cortas puede no mostrar la
#' variante causal: por pseudogenes, duplicaciones segmentarias o cambios en
#' el numero de copias, o porque la enfermedad suele deberse a variantes
#' somaticas o en mosaico. [iei_priorizar()], [iei_gen()] e
#' [iei_buscar_fenotipo()] muestran la advertencia, y [iei_filtrar()] no
#' descarta los posibles mosaicos de los genes de tipo `"mosaico"`.
#'
#' @format Un tibble con columnas `gen`, `tipo` (`"pseudogen"`,
#'   `"duplicacion segmentaria"`, `"numero de copias"` o `"mosaico"`) y
#'   `motivo`.
#' @source Curacion propia en `data-raw/genes_dificiles.csv`, a partir de la
#'   literatura de cada gen.
"genes_dificiles"

#' Fenotipos por enfermedad de los genes del panel
#'
#' Signos HPO de cada enfermedad (OMIM u Orphanet) asociada a un gen del
#' panel, con su frecuencia y los signos que la enfermedad excluye. Es la base
#' de la similitud por enfermedad de [iei_buscar_fenotipo()] e
#' [iei_priorizar()], y de los casos simulados de [iei_validar_fenotipo()].
#'
#' @format Un tibble con columnas `gen`, `enfermedad_id`, `enfermedad`,
#'   `herencia` (de la anotacion de la enfermedad; `NA` si HPO no la indica),
#'   `hpo_id`, `frecuencia` (proporcion de pacientes con el signo, `NA` si no
#'   se conoce) y `excluido` (`TRUE` cuando la enfermedad no presenta ese
#'   signo, anotacion NOT).
#' @source `phenotype.hpoa` y `genes_to_disease.txt` de la Human Phenotype
#'   Ontology, version indicada en [fuentes_datos].
#'   <https://github.com/obophenotype/human-phenotype-ontology/releases>
"hpo_enfermedades"

#' Variantes reales para practicar
#'
#' Doce variantes descritas en ClinVar y en la literatura de su enfermedad,
#' en genes conocidos y en otros menos conocidos del panel. Sirven para
#' practicar sin tener una secuencia propia. Cada fila trae la variante en
#' el formato que acepta [iei_variante()], el transcrito MANE Select al que
#' se refiere, su identificador de ClinVar, un genotipo y tres signos HPO de
#' la enfermedad tomados de [hpo_enfermedades]. La columna `que_ensena` dice
#' que aspecto del analisis muestra cada una.
#'
#' Dos filas son contraejemplos a proposito. RAG1 c.1186C>T aparece en
#' heterocigosis en un gen recesivo, asi que representa a una portadora y no
#' deberia subir a prioridad alta. TNFRSF13B c.310T>C es un factor de riesgo
#' que tambien se encuentra en personas sanas y tiene clasificaciones en
#' conflicto en ClinVar.
#'
#' El banco sirve para practicar el flujo del paquete. Doce casos no alcanzan
#' para medir sensibilidad ni especificidad.
#'
#' @format Un tibble con una fila por variante y columnas `gen`, `variante`
#'   (gen y cambio en HGVS), `transcrito` (RefSeq MANE Select del cambio),
#'   `clinvar_id` (identificador de variacion de ClinVar, `VCV` seguido de
#'   nueve digitos), `proteina` (cambio en la proteina, `NA` si no aplica),
#'   `tipo` (consecuencia), `genotipo` (`"het"`, `"hom_alt"` o `"hemi"`),
#'   `herencia`, `enfermedad`, `enfermedad_id` (OMIM), `fenotipo` (codigos
#'   HPO separados por `;`) y `que_ensena`.
#' @source Curacion propia en `data-raw/variantes_practica.csv`. El transcrito
#'   y el identificador de ClinVar de cada fila se verificaron contra ClinVar
#'   el 2026-10-20. La clasificacion de ClinVar puede cambiar;
#'   [iei_variante()] la consulta al momento. Los signos salen de la Human
#'   Phenotype Ontology. <https://www.ncbi.nlm.nih.gov/clinvar/>
#' @examples
#' variantes_practica[, c("gen", "variante", "clinvar_id")]
#' v <- variantes_practica[variantes_practica$gen == "STAT1", ]
#' strsplit(v$fenotipo, ";")[[1]]
#' \dontrun{
#' iei_variante(v$variante, genotipo = v$genotipo,
#'              fenotipo = strsplit(v$fenotipo, ";")[[1]])
#' }
"variantes_practica"
