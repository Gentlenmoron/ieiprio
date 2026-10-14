# Anotar variantes con Ensembl VEP

Consulta Ensembl VEP para cada variante unica y agrega su consecuencia,
nomenclatura HGVS, exon, frecuencias en gnomAD (global y del grupo mas
alto), significancia en ClinVar y predictores. Usa el servidor de GRCh37
o GRCh38 segun el build. Las respuestas se guardan en una cache local
por release de Ensembl, asi que repetir el analisis no vuelve a
consultar y nunca se mezclan anotaciones de releases distintos.

## Usage

``` r
iei_anotar(
  x,
  cache = tools::R_user_dir("ieiprio", "cache"),
  lote = 200,
  normalizar = TRUE,
  clinvar = TRUE,
  gnomad = FALSE,
  vep_json = NULL
)
```

## Arguments

- x:

  Tabla devuelta por
  [`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md),
  [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md),
  [`iei_variantes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variantes.md)
  o
  [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md).

- cache:

  Carpeta de la cache, o `FALSE` para no usarla.

- lote:

  Variantes por consulta. La API acepta hasta 200.

- normalizar:

  Si `TRUE`, normaliza inserciones y deleciones.

- clinvar:

  Si `TRUE`, descarga las variantes patogenicas de ClinVar de cada gen
  afectado para obtener las estrellas de revision y los criterios PS1,
  PM5 y PM1 (ver
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)).

- gnomad:

  Si `TRUE`, consulta en gnomAD los homocigotos, hemicigotos, la
  frecuencia filtrada al 95 % y la restriccion del gen. Es lento (gnomAD
  limita las consultas); conviene usarlo despues de filtrar, o aplicar
  [`iei_gnomad()`](https://gentlenmoron.github.io/ieiprio/reference/iei_gnomad.md)
  solo a las candidatas.

- vep_json:

  Ruta a la salida JSON de un VEP local, en lugar de la API.

## Value

La tabla `x` con las columnas `transcrito`, `mane`, `exon`,
`consecuencia`, `impacto`, `hgvs_c`, `hgvs_p`, `af_gnomad` (global),
`af_grpmax` y `grupo_grpmax` (grupo poblacional mas alto), `clinvar`,
`proteina_id`, `pos_proteina`, los predictores `alphamissense`,
`am_clase`, `revel` y `spliceai`, y `vep_release`. Con `clinvar = TRUE`
agrega `clinvar_estrellas`, `ps1`, `pm5` y `pm1`; con `gnomad = TRUE`,
las columnas de
[`iei_gnomad()`](https://gentlenmoron.github.io/ieiprio/reference/iei_gnomad.md).
REVEL es de uso no comercial.

## Details

Antes de anotar, las inserciones y deleciones se normalizan (recorte de
bases comunes y alineamiento a la izquierda contra la referencia), para
que coincidan con gnomAD y ClinVar aunque el VCF las escriba de otra
forma.

Para cada variante se elige un transcrito, en este orden de preferencia.
Primero uno del mismo gen de la columna `gen` (si existe), luego el MANE
Select, luego el canonico de Ensembl y por ultimo el de mayor impacto.

**Privacidad.** En modo en linea cada variante se envia a Ensembl (Reino
Unido) y, con `clinvar = TRUE`, el gen se consulta en NCBI (Estados
Unidos). Para no enviar nada, corre VEP en tu computadora con el archivo
de
[`iei_vep_entrada()`](https://gentlenmoron.github.io/ieiprio/reference/iei_vep_entrada.md)
y pasa su salida en `vep_json`, con `clinvar = FALSE`. Con
`options(ieiprio.sin_red = TRUE)` el paquete se niega a hacer cualquier
consulta externa.

## Examples

``` r
if (FALSE) { # \dontrun{
anot <- iei_anotar(iei_filtrar(iei_leer_vcf("paciente01.vcf.gz", sexo = "M")))
} # }
```
