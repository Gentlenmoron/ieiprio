# Cromatogramas Sanger de ejemplo

Genera dos cromatogramas `.ab1` (directo y reverso) de un amplicon de
500 pb de BTK que contiene la variante patogenica real
NM_000061.3:c.1632-2A\>G (ClinVar 265455), la misma del VCF de ejemplo.
Sirven para probar
[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
y
[`iei_cromatograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cromatograma.md)
sin tener archivos propios. La secuencia sale de la referencia GRCh38 de
BTK que trae el paquete (`inst/extdata/btk_referencia.fa`, bajada de
Ensembl), asi que no necesita conexion.

## Usage

``` r
iei_ejemplo_sanger(genotipo = c("hemi", "het"), dir = tempfile("sanger_"))
```

## Arguments

- genotipo:

  `"hemi"` para un varon afectado o `"het"` para una portadora, que se
  ve como pico doble.

- dir:

  Carpeta donde escribir los archivos.

## Value

Las rutas de los dos archivos, de forma invisible.

## Details

Los picos son simulados: el cromatograma es limpio salvo en los
extremos, donde la calidad cae como en una lectura real.

## Examples

``` r
ab1 <- iei_ejemplo_sanger()
#> Cromatogramas de ejemplo en /tmp/RtmpKS5CK3/sanger_1986367a2473.
ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
v <- iei_sanger(ab1, gen = "BTK", sexo = "M", referencia = ref)
#> BTK_c1632-2_F.ab1: hebra +, 461 bases alineadas en X:101353757-101354217,
#> identidad 99.8%.
#> BTK_c1632-2_R.ab1: hebra -, 461 bases alineadas en X:101353758-101354218,
#> identidad 99.8%.
v
#> # A tibble: 1 × 16
#>   muestra gen   chr       pos id    ref   alt   genotipo  qual    dp    gq    ab
#>   <chr>   <chr> <chr>   <int> <chr> <chr> <chr> <chr>    <dbl> <dbl> <dbl> <dbl>
#> 1 sanger  BTK   X      1.01e8 NA    T     C     hemi        45    NA    NA     1
#> # ℹ 4 more variables: filtro <chr>, build <chr>, lecturas <int>, fuente <chr>
```
