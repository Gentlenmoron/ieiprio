# Consultar el panel de genes

Devuelve el panel de genes de errores innatos de la inmunidad, filtrado
por nivel de evidencia, categoria IUIS y modo de herencia.

## Usage

``` r
iei_panel(categoria = NULL, herencia = NULL, evidencia = "verde")
```

## Arguments

- categoria:

  Vector de enteros entre 1 y 10, o `NULL` para todas. Solo devuelve
  genes con categoria IUIS curada.

- herencia:

  Vector con `"AR"`, `"AD"`, `"XL"` o `"MT"`, o `NULL` para todos. Un
  gen `"AR/AD"` aparece si se pide `"AR"` o `"AD"`.

- evidencia:

  Niveles de evidencia de PanelApp a incluir: `"verde"`, `"ambar"`,
  `"rojo"`. Por defecto solo verde, que son los genes con asociacion
  diagnostica establecida.

## Value

Un tibble con la estructura de
[iuis_panel](https://gentlenmoron.github.io/ieiprio/reference/iuis_panel.md).

## Examples

``` r
iei_panel()
#> # A tibble: 371 × 8
#>    gen    hgnc_id herencia herencia_panelapp evidencia fenotipos categorias_iuis
#>    <chr>  <chr>   <chr>    <chr>             <chr>     <chr>     <chr[1d]>      
#>  1 ACD    HGNC:2… AR/AD    BOTH monoallelic… verde     Dyskerat… NA             
#>  2 ACP5   HGNC:1… AR       BIALLELIC, autos… verde     Spondylo… NA             
#>  3 ADA    HGNC:1… AR       BIALLELIC, autos… verde     Severe c… 1              
#>  4 ADA2   HGNC:1… AR       BIALLELIC, autos… verde     Vasculit… NA             
#>  5 ADAM17 HGNC:1… AR       BIALLELIC, autos… verde     inflamma… NA             
#>  6 ADAR   HGNC:2… AR       BIALLELIC, autos… verde     Aicardi-… NA             
#>  7 AGR2   HGNC:3… AR       BIALLELIC, autos… verde     Cystic f… NA             
#>  8 AICDA  HGNC:1… AR/AD    BOTH monoallelic… verde     Immunode… NA             
#>  9 AIRE   HGNC:3… AR/AD    BOTH monoallelic… verde     Autoimmu… NA             
#> 10 AK2    HGNC:3… AR       BIALLELIC, autos… verde     Reticula… NA             
#> # ℹ 361 more rows
#> # ℹ 1 more variable: ganancia_funcion <lgl>
iei_panel(categoria = 3)
#> # A tibble: 2 × 8
#>   gen    hgnc_id  herencia herencia_panelapp evidencia fenotipos categorias_iuis
#>   <chr>  <chr>    <chr>    <chr>             <chr>     <chr>     <chr[1d]>      
#> 1 BTK    HGNC:11… XL       X-LINKED: hemizy… verde     Agammagl… 3              
#> 2 PIK3CD HGNC:89… AR/AD    BOTH monoallelic… verde     Immunode… 3              
#> # ℹ 1 more variable: ganancia_funcion <lgl>
iei_panel(herencia = "XL", evidencia = c("verde", "ambar"))
#> # A tibble: 25 × 8
#>    gen    hgnc_id herencia herencia_panelapp evidencia fenotipos categorias_iuis
#>    <chr>  <chr>   <chr>    <chr>             <chr>     <chr>     <chr[1d]>      
#>  1 ATP6A… HGNC:8… XL       X-LINKED: hemizy… verde     Immunode… NA             
#>  2 BTK    HGNC:1… XL       X-LINKED: hemizy… verde     Agammagl… 3              
#>  3 CD40LG HGNC:1… XL       X-LINKED: hemizy… verde     Immunode… 1              
#>  4 CFP    HGNC:8… XL       X-LINKED: hemizy… verde     Properdi… NA             
#>  5 CYBB   HGNC:2… XL       X-LINKED: hemizy… verde     Immunode… 5              
#>  6 DKC1   HGNC:2… XL       X-LINKED: hemizy… verde     Dyskerat… NA             
#>  7 DOCK11 HGNC:2… XL       X-LINKED: hemizy… verde     Autoinfl… NA             
#>  8 ELF4   HGNC:3… XL       X-LINKED: hemizy… verde     Autoinfl… NA             
#>  9 FOXP3  HGNC:6… XL       X-LINKED: hemizy… verde     Immunody… 4              
#> 10 G6PD   HGNC:4… XL       X-LINKED: hemizy… verde     Glucose-… NA             
#> # ℹ 15 more rows
#> # ℹ 1 more variable: ganancia_funcion <lgl>
```
