# Anotar variantes con Ensembl VEP

Consulta la API REST de Ensembl VEP para cada variante unica y agrega su
consecuencia, nomenclatura HGVS, frecuencia en gnomAD y significancia
clinica en ClinVar. Usa el servidor de GRCh37 o GRCh38 segun el build
del VCF. Las respuestas se guardan en una cache local, asi que volver a
anotar las mismas variantes no consulta la API de nuevo.

## Usage

``` r
iei_anotar(x, cache = tools::R_user_dir("ieiprio", "cache"), lote = 200)
```

## Arguments

- x:

  Tabla devuelta por
  [`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md)
  o
  [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md).

- cache:

  Carpeta de la cache, o `FALSE` para no usarla.

- lote:

  Variantes por consulta. La API acepta hasta 200.

## Value

La tabla `x` con las columnas `transcrito`, `mane`, `consecuencia`,
`impacto`, `hgvs_c`, `hgvs_p`, `af_gnomad`, `clinvar`, `proteina_id` y
`pos_proteina`. Si la columna `id` estaba vacia, se completa con el rsID
de Ensembl. Conserva la bitacora y agrega el atributo `vep_release`.

## Details

Para cada variante se elige un transcrito, en este orden de preferencia.
Primero uno del mismo gen de la columna `gen` (si existe), luego el MANE
Select, luego el canonico de Ensembl y por ultimo el de mayor impacto.

## Examples

``` r
if (FALSE) { # \dontrun{
anot <- iei_anotar(iei_filtrar(iei_leer_vcf("paciente01.vcf.gz", sexo = "M")))
} # }
```
