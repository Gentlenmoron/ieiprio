# Fuentes y versiones de los datos del paquete

Muestra de donde salio cada dataset incluido, con su version, URL, fecha
de descarga y huella md5. Sirve para citar y reproducir un analisis.

## Usage

``` r
iei_fuentes()
```

## Value

Un tibble, ver
[fuentes_datos](https://gentlenmoron.github.io/ieiprio/reference/fuentes_datos.md).

## Examples

``` r
iei_fuentes()
#> # A tibble: 3 × 6
#>   fuente        detalle                                version url   fecha md5  
#>   <chr>         <chr>                                  <chr>   <chr> <chr> <chr>
#> 1 PanelApp      Primary immunodeficiency or monogenic… 8.78    http… 2026… 6fe8…
#> 2 HPO           genes_to_phenotype.txt                 v2026-… http… 2026… 6dac…
#> 3 Curacion IUIS data-raw/iuis_curado.csv               git     http… 2026… d554…
```
