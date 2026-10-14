# Ruta del VCF de ejemplo

El paquete incluye un VCF de un paciente ficticio (EJEMPLO01, varon) con
cinco variantes reales tomadas de ClinVar en genes del panel y unas 4000
variantes de fondo sinteticas fuera de esos genes. Se genera con
`data-raw/ejemplo_vcf.R`; las variantes de ClinVar estan fijadas por su
identificador en `data-raw/ejemplo_clinvar.csv`.

## Usage

``` r
iei_ejemplo_vcf()
```

## Value

La ruta al archivo `ejemplo.vcf.gz`.

## Examples

``` r
if (FALSE) { # \dontrun{
v <- iei_leer_vcf(iei_ejemplo_vcf(), sexo = "M")
} # }
```
