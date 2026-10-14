# Buscar terminos HPO por texto

Encuentra terminos de la Human Phenotype Ontology a partir de una
palabra o frase, en espanol (para los signos mas frecuentes en
inmunodeficiencias) o en ingles (todos los terminos y sus sinonimos).
Sirve para armar el `fenotipo` de
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
o
[`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md).

## Usage

``` r
iei_buscar_hpo(texto, n = 15)
```

## Arguments

- texto:

  Palabra o frase, por ejemplo `"neumonia"` o `"candidiasis"`.

- n:

  Numero maximo de terminos a devolver.

## Value

Un tibble con `hpo_id`, `termino`, `genes_panel` (cuantos genes del
panel lo tienen anotado) e `ic` (contenido de informacion: mas alto es
mas especifico).

## Examples

``` r
iei_buscar_hpo("neumonia")
#> # A tibble: 15 × 4
#>    hpo_id     termino                                       genes_panel    ic
#>    <chr>      <chr>                                               <int> <dbl>
#>  1 HP:0006532 Recurrent pneumonia                                    83  3.49
#>  2 HP:0002090 Pneumonia                                              68  2.72
#>  3 HP:0020102 Pneumocystis jirovecii pneumonia                       12  6.09
#>  4 HP:0006538 Recurrent bronchopulmonary infections                  10  5.74
#>  5 HP:0031950 Usual interstitial pneumonia                            6  5.86
#>  6 HP:0002100 Recurrent aspiration pneumonia                          4  5.04
#>  7 HP:0005366 Recurrent streptococcus pneumoniae infections           4  7.18
#>  8 HP:0011949 Acute infectious pneumonia                              3  6.78
#>  9 HP:0006527 Lymphocytic interstitial pneumonia                      3  7.18
#> 10 HP:0011945 Bronchiolitis obliterans organizing pneumonia           2  7.88
#> 11 HP:0011951 Aspiration pneumonia                                    1  3.9 
#> 12 HP:0033670 Organizing pneumonia                                    1  8.57
#> 13 HP:0033214 Recurrent viral pneumonia                               1  8.57
#> 14 HP:5210428 HSV pneumonia                                           0  8.57
#> 15 HP:5210302 VZV pneumonitis                                         0  8.57
iei_buscar_hpo("candidiasis")
#> # A tibble: 15 × 4
#>    hpo_id     termino                             genes_panel    ic
#>    <chr>      <chr>                                     <int> <dbl>
#>  1 HP:0002728 Recurrent mucocutaneous candidiasis          28  4.49
#>  2 HP:0009098 Recurrent oral thrush                        36  4.96
#>  3 HP:0005401 Recurrent Candida infection                  16  5.58
#>  4 HP:0012203 Unusual fungal nail infection                12  5.35
#>  5 HP:0033351 Candida esophagitis                           5  6.96
#>  6 HP:0012204 Recurrent vulvovaginal candidiasis            2  7.47
#>  7 HP:0005411 Chronic intestinal candidiasis                1  8.57
#>  8 HP:0033768 Penile thrush                                 0  8.57
#>  9 HP:5210253 Severe oral thrush                            0  8.57
#> 10 HP:5210167 Hepatic candidiasis                           0  8.57
#> 11 HP:5210252 Unusual oral thrush                           0  4.96
#> 12 HP:5210198 Fungal kidney abscess                         0  8.57
#> 13 HP:5210241 Recurrent onychomycosis                       0  8.57
#> 14 HP:6001283 Severe Candida infection                      0  8.57
#> 15 HP:5210102 Unusual Candida infection                     0  4.66
```
