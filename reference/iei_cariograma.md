# Cariograma del panel

Dibuja los 24 cromosomas uno al lado del otro con los genes del panel
marcados a la izquierda en el color de su categoria IUIS. Con un
resultado de
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
marca las variantes a la derecha, del color de su prioridad, y con una
busqueda de
[`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md)
marca los genes candidatos, del color de su herencia. Los mejores se
rotulan.

## Usage

``` r
iei_cariograma(
  x = NULL,
  build = c("GRCh38", "GRCh37"),
  etiquetas = 12,
  titulo = NULL
)
```

## Arguments

- x:

  `NULL`, un resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  o una busqueda de
  [`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md).

- build:

  Version del genoma.

- etiquetas:

  Cuantos genes rotular.

- titulo:

  Titulo del grafico; `NULL` usa uno por defecto.

## Value

Un objeto ggplot.

## See also

[`iei_circos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_circos.md)

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  iei_cariograma()
  iei_cariograma(iei_buscar_fenotipo("candidiasis mucocutanea cronica", n = 10))
}
#> ℹ Fenotipo interpretado como:
#> • candidiasis mucocutanea cronica → HP:0002728 Recurrent mucocutaneous
#>   candidiasis
```
