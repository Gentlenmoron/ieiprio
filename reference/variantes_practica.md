# Variantes reales para practicar

Doce variantes descritas en ClinVar y en la literatura de su enfermedad,
en genes conocidos y en otros menos conocidos del panel. Sirven para
practicar sin tener una secuencia propia. Cada fila trae la variante en
el formato que acepta
[`iei_variante()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variante.md),
el transcrito MANE Select al que se refiere, su identificador de
ClinVar, un genotipo y tres signos HPO de la enfermedad tomados de
[hpo_enfermedades](https://gentlenmoron.github.io/ieiprio/reference/hpo_enfermedades.md).
La columna `que_ensena` dice que aspecto del analisis muestra cada una.

## Usage

``` r
variantes_practica
```

## Format

Un tibble con una fila por variante y columnas `gen`, `variante` (gen y
cambio en HGVS), `transcrito` (RefSeq MANE Select del cambio),
`clinvar_id` (identificador de variacion de ClinVar, `VCV` seguido de
nueve digitos), `proteina` (cambio en la proteina, `NA` si no aplica),
`tipo` (consecuencia), `genotipo` (`"het"`, `"hom_alt"` o `"hemi"`),
`herencia`, `enfermedad`, `enfermedad_id` (OMIM), `fenotipo` (codigos
HPO separados por `;`) y `que_ensena`.

## Source

Curacion propia en `data-raw/variantes_practica.csv`. El transcrito y el
identificador de ClinVar de cada fila se verificaron contra ClinVar el
2026-10-20. La clasificacion de ClinVar puede cambiar;
[`iei_variante()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variante.md)
la consulta al momento. Los signos salen de la Human Phenotype Ontology.
<https://www.ncbi.nlm.nih.gov/clinvar/>

## Details

Dos filas son contraejemplos a proposito. RAG1 c.1186C\>T aparece en
heterocigosis en un gen recesivo, asi que representa a una portadora y
no deberia subir a prioridad alta. TNFRSF13B c.310T\>C es un factor de
riesgo que tambien se encuentra en personas sanas y tiene
clasificaciones en conflicto en ClinVar.

El banco sirve para practicar el flujo del paquete. Doce casos no
alcanzan para medir sensibilidad ni especificidad.

## Examples

``` r
variantes_practica[, c("gen", "variante", "clinvar_id")]
#> # A tibble: 12 × 3
#>    gen       variante           clinvar_id  
#>    <chr>     <chr>              <chr>       
#>  1 BTK       BTK c.1632-2A>G    VCV000265455
#>  2 RAG1      RAG1 c.1186C>T     VCV000013144
#>  3 ADA       ADA c.646G>A       VCV000001968
#>  4 STAT3     STAT3 c.1144C>T    VCV000018304
#>  5 PIK3CD    PIK3CD c.3061G>A   VCV000088675
#>  6 TNFRSF13B TNFRSF13B c.310T>C VCV000005302
#>  7 FOXP3     FOXP3 c.1150G>A    VCV000011410
#>  8 AIRE      AIRE c.769C>T      VCV000003307
#>  9 CYBB      CYBB c.676C>T      VCV000010929
#> 10 NCF1      NCF1 c.75_76del    VCV000002249
#> 11 STAT1     STAT1 c.820C>T     VCV000030083
#> 12 MEFV      MEFV c.2080A>G     VCV000002538
v <- variantes_practica[variantes_practica$gen == "STAT1", ]
strsplit(v$fenotipo, ";")[[1]]
#> [1] "HP:0002728" "HP:0009098" "HP:0012203"
if (FALSE) { # \dontrun{
iei_variante(v$variante, genotipo = v$genotipo,
             fenotipo = strsplit(v$fenotipo, ";")[[1]])
} # }
```
