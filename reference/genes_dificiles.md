# Genes dificiles de analizar con lecturas cortas

Genes del panel en los que un VCF de lecturas cortas puede no mostrar la
variante causal: por pseudogenes, duplicaciones segmentarias o cambios
en el numero de copias, o porque la enfermedad suele deberse a variantes
somaticas o en mosaico.
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md),
[`iei_gen()`](https://gentlenmoron.github.io/ieiprio/reference/iei_gen.md)
e
[`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md)
muestran la advertencia, y
[`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md)
no descarta los posibles mosaicos de los genes de tipo `"mosaico"`.

## Usage

``` r
genes_dificiles
```

## Format

Un tibble con columnas `gen`, `tipo` (`"pseudogen"`,
`"duplicacion segmentaria"`, `"numero de copias"` o `"mosaico"`) y
`motivo`.

## Source

Curacion propia en `data-raw/genes_dificiles.csv`, a partir de la
literatura de cada gen.
