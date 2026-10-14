# Convertir a GRanges

Devuelve la tabla como un objeto `GRanges` de Bioconductor, con todas
las columnas como metadatos, para seguir el analisis con GenomicRanges,
VariantAnnotation o plyranges.

## Usage

``` r
iei_granges(x)
```

## Arguments

- x:

  Tabla de ieiprio con `chr`, `pos` y `ref`.

## Value

Un objeto `GRanges` de GenomicRanges, que se instala aparte con
`BiocManager::install("GenomicRanges")`.

## Examples

``` r
if (FALSE) { # \dontrun{
gr <- iei_granges(prior)
} # }
```
