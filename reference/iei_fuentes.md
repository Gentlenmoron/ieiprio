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
#> # A tibble: 8 × 6
#>   fuente                 detalle                       version url   fecha md5  
#>   <chr>                  <chr>                         <chr>   <chr> <chr> <chr>
#> 1 PanelApp               Primary immunodeficiency or … 8.78    http… 2026… 6fe8…
#> 2 HPO                    genes_to_phenotype.txt        v2026-… http… 2026… 6dac…
#> 3 HPO ontologia          hp.obo                        v2026-… http… 2026… 9762…
#> 4 HPO enfermedades       phenotype.hpoa                v2026-… http… 2026… b412…
#> 5 HPO genes a enfermedad genes_to_disease.txt          v2026-… http… 2026… 337e…
#> 6 Curacion IUIS          data-raw/iuis_curado.csv      git     http… 2026… b9d0…
#> 7 Mecanismo GOF          data-raw/iuis_mecanismo.csv   git     http… 2026… c88b…
#> 8 Genes dificiles        data-raw/genes_dificiles.csv  git     http… 2026… 076c…
```
