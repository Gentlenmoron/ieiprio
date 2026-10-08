# Buscar genes en el panel

Busca uno o mas genes por su simbolo oficial o por un alias o simbolo
anterior (por ejemplo `"TACI"` para TNFRSF13B). No distingue mayusculas.
Busca en todo el panel, sin filtrar por evidencia.

## Usage

``` r
iei_buscar_gen(genes)
```

## Arguments

- genes:

  Vector de caracteres con simbolos o alias.

## Value

Un tibble con la columna `consulta` seguida de las columnas de
[iuis_panel](https://gentlenmoron.github.io/ieiprio/reference/iuis_panel.md).
Los genes que no estan en el panel se omiten con un aviso.

## Examples

``` r
iei_buscar_gen("BTK")
#> # A tibble: 1 × 9
#>   consulta gen   hgnc_id   herencia herencia_panelapp        evidencia fenotipos
#>   <chr>    <chr> <chr>     <chr>    <chr>                    <chr>     <chr>    
#> 1 BTK      BTK   HGNC:1133 XL       X-LINKED: hemizygous mu… verde     Agammagl…
#> # ℹ 2 more variables: categorias_iuis <chr[1d]>, ganancia_funcion <lgl>
iei_buscar_gen(c("taci", "STAT1"))
#> # A tibble: 2 × 9
#>   consulta gen       hgnc_id    herencia herencia_panelapp   evidencia fenotipos
#>   <chr>    <chr>     <chr>      <chr>    <chr>               <chr>     <chr>    
#> 1 taci     TNFRSF13B HGNC:18153 AR/AD    BOTH monoallelic a… rojo      Immunode…
#> 2 STAT1    STAT1     HGNC:11362 AR/AD    BOTH monoallelic a… verde     Immunode…
#> # ℹ 2 more variables: categorias_iuis <chr[1d]>, ganancia_funcion <lgl>
```
