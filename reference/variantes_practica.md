# Variantes reales para practicar

Doce variantes descritas en ClinVar y en la literatura de su enfermedad,
en genes conocidos y en otros menos conocidos del panel. Sirven para
practicar sin tener una secuencia propia: cada fila trae la variante en
el formato que acepta
[`iei_variante()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variante.md),
el genotipo con que suele causar la enfermedad y tres signos HPO de esa
enfermedad, tomados de
[hpo_enfermedades](https://gentlenmoron.github.io/ieiprio/reference/hpo_enfermedades.md).
La columna `que_ensena` dice que aspecto del analisis muestra cada una.

## Usage

``` r
variantes_practica
```

## Format

Un tibble con una fila por variante y columnas `gen`, `variante` (gen y
cambio en HGVS), `proteina` (cambio en la proteina, `NA` si no aplica),
`tipo` (consecuencia), `genotipo` (`"het"`, `"hom_alt"` o `"hemi"`),
`herencia`, `enfermedad`, `enfermedad_id` (OMIM), `fenotipo` (codigos
HPO separados por `;`) y `que_ensena`.

## Source

Curacion propia en `data-raw/variantes_practica.csv`. Las variantes se
pueden revisar en ClinVar y los signos salen de la Human Phenotype
Ontology. La clasificacion de ClinVar puede cambiar;
[`iei_variante()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variante.md)
la consulta al momento.

## Examples

``` r
variantes_practica[, c("gen", "variante", "enfermedad")]
#> # A tibble: 12 × 3
#>    gen       variante           enfermedad                                      
#>    <chr>     <chr>              <chr>                                           
#>  1 BTK       BTK c.1632-2A>G    Agammaglobulinemia ligada al X                  
#>  2 RAG1      RAG1 c.1186C>T     Inmunodeficiencia combinada grave T-B-NK+       
#>  3 ADA       ADA c.646G>A       Inmunodeficiencia combinada grave por déficit d…
#>  4 STAT3     STAT3 c.1144C>T    Síndrome de hiper-IgE autosómico dominante      
#>  5 PIK3CD    PIK3CD c.3061G>A   Síndrome de PI3K delta activada (APDS1)         
#>  6 TNFRSF13B TNFRSF13B c.310T>C Inmunodeficiencia común variable 2              
#>  7 FOXP3     FOXP3 c.1150G>A    Síndrome IPEX                                   
#>  8 AIRE      AIRE c.769C>T      Poliendocrinopatía autoinmune tipo 1 (APECED)   
#>  9 CYBB      CYBB c.676C>T      Enfermedad granulomatosa crónica ligada al X    
#> 10 NCF1      NCF1 c.75_76del    Enfermedad granulomatosa crónica autosómica rec…
#> 11 STAT1     STAT1 c.820C>T     Candidiasis mucocutánea crónica por ganancia de…
#> 12 MEFV      MEFV c.2080A>G     Fiebre mediterránea familiar                    
v <- variantes_practica[variantes_practica$gen == "STAT1", ]
strsplit(v$fenotipo, ";")[[1]]
#> [1] "HP:0002728" "HP:0009098" "HP:0012203"
if (FALSE) { # \dontrun{
iei_variante(v$variante, genotipo = v$genotipo,
             fenotipo = strsplit(v$fenotipo, ";")[[1]])
} # }
```
