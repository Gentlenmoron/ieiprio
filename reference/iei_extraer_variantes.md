# Encontrar variantes escritas en un texto

Busca descripciones de variantes en texto libre: nomenclatura HGVS de
transcrito (`NM_000061.3:c.1632-2A>G`, tambien con el gen entre
parentesis como aparece en ClinVar), HGVS genomica
(`NC_000023.11:g.101353990T>C`), identificadores rsID y el gen seguido
del cambio (`BTK c.1632-2A>G` o `BTK p.Arg562Trp`). Si en la misma linea
aparece heterocigoto, homocigoto o hemicigoto, tambien lo registra.

## Usage

``` r
iei_extraer_variantes(texto)
```

## Arguments

- texto:

  Vector de texto, por ejemplo la salida de
  [`iei_leer_texto()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_texto.md).

## Value

Un tibble con `entrada` (lo encontrado, normalizado para Ensembl),
`tipo`, `genotipo` y `linea`.

## Examples

``` r
iei_extraer_variantes(c(
  "Se identifico la variante NM_000061.3(BTK):c.1632-2A>G en hemicigosis.",
  "Tambien rs121908711 y STAT1 p.Arg274Gln en heterocigosis."
))
#> # A tibble: 3 × 4
#>   entrada                 tipo            genotipo linea
#>   <chr>                   <chr>           <chr>    <int>
#> 1 NM_000061.3:c.1632-2A>G hgvs_transcrito hemi         1
#> 2 rs121908711             rsid            het          2
#> 3 STAT1:p.Arg274Gln       gen_y_cambio    het          2
```
