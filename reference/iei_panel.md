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
#> # A tibble: 371 × 9
#>    gen    hgnc_id herencia herencia_panelapp evidencia fenotipos categorias_iuis
#>    <chr>  <chr>   <chr>    <chr>             <chr>     <chr>     <chr>          
#>  1 ACD    HGNC:2… AR/AD    BOTH monoallelic… verde     Dyskerat… 9              
#>  2 ACP5   HGNC:1… AR       BIALLELIC, autos… verde     Spondylo… 7              
#>  3 ADA    HGNC:1… AR       BIALLELIC, autos… verde     Severe c… 1              
#>  4 ADA2   HGNC:1… AR       BIALLELIC, autos… verde     Vasculit… 7              
#>  5 ADAM17 HGNC:1… AR       BIALLELIC, autos… verde     inflamma… 7              
#>  6 ADAR   HGNC:2… AR       BIALLELIC, autos… verde     Aicardi-… 7              
#>  7 AGR2   HGNC:3… AR       BIALLELIC, autos… verde     Cystic f… NA             
#>  8 AICDA  HGNC:1… AR/AD    BOTH monoallelic… verde     Immunode… NA             
#>  9 AIRE   HGNC:3… AR/AD    BOTH monoallelic… verde     Autoimmu… NA             
#> 10 AK2    HGNC:3… AR       BIALLELIC, autos… verde     Reticula… 1              
#> # ℹ 361 more rows
#> # ℹ 2 more variables: mecanismo <chr>, ganancia_funcion <lgl>
iei_panel(categoria = 3)
#> # A tibble: 36 × 9
#>    gen    hgnc_id herencia herencia_panelapp evidencia fenotipos categorias_iuis
#>    <chr>  <chr>   <chr>    <chr>             <chr>     <chr>     <chr>          
#>  1 ATP6A… HGNC:8… XL       X-LINKED: hemizy… verde     Immunode… 3              
#>  2 BLNK   HGNC:1… AR       BIALLELIC, autos… verde     Agammagl… 3              
#>  3 BTK    HGNC:1… XL       X-LINKED: hemizy… verde     Agammagl… 3              
#>  4 CARD11 HGNC:1… AR/AD    BOTH monoallelic… verde     Immunode… 1;2;3          
#>  5 CD19   HGNC:1… AR       BIALLELIC, autos… verde     Immunode… 3              
#>  6 CD79A  HGNC:1… AR       BIALLELIC, autos… verde     Agammagl… 3              
#>  7 CD79B  HGNC:1… AR       BIALLELIC, autos… verde     Agammagl… 3              
#>  8 CD81   HGNC:1… AR       BIALLELIC, autos… verde     Immunode… 3              
#>  9 CR2    HGNC:2… AR       BIALLELIC, autos… verde     ?Immunod… 3              
#> 10 FNIP1  HGNC:2… AR       BIALLELIC, autos… verde     Primary … 3;5            
#> # ℹ 26 more rows
#> # ℹ 2 more variables: mecanismo <chr>, ganancia_funcion <lgl>
iei_panel(herencia = "XL", evidencia = c("verde", "ambar"))
#> # A tibble: 25 × 9
#>    gen    hgnc_id herencia herencia_panelapp evidencia fenotipos categorias_iuis
#>    <chr>  <chr>   <chr>    <chr>             <chr>     <chr>     <chr>          
#>  1 ATP6A… HGNC:8… XL       X-LINKED: hemizy… verde     Immunode… 3              
#>  2 BTK    HGNC:1… XL       X-LINKED: hemizy… verde     Agammagl… 3              
#>  3 CD40LG HGNC:1… XL       X-LINKED: hemizy… verde     Immunode… 1              
#>  4 CFP    HGNC:8… XL       X-LINKED: hemizy… verde     Properdi… 8              
#>  5 CYBB   HGNC:2… XL       X-LINKED: hemizy… verde     Immunode… 5;6            
#>  6 DKC1   HGNC:2… XL       X-LINKED: hemizy… verde     Dyskerat… 2              
#>  7 DOCK11 HGNC:2… XL       X-LINKED: hemizy… verde     Autoinfl… 4;7            
#>  8 ELF4   HGNC:3… XL       X-LINKED: hemizy… verde     Autoinfl… 7              
#>  9 FOXP3  HGNC:6… XL       X-LINKED: hemizy… verde     Immunody… 4              
#> 10 G6PD   HGNC:4… XL       X-LINKED: hemizy… verde     Glucose-… 5              
#> # ℹ 15 more rows
#> # ℹ 2 more variables: mecanismo <chr>, ganancia_funcion <lgl>
```
