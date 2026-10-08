# Analizar un VCF de principio a fin

Corre todo el flujo en una sola llamada. Lee el VCF, filtra por genes
del panel y calidad, anota con Ensembl VEP, filtra por frecuencia,
prioriza y genera un reporte HTML por muestra.

## Usage

``` r
iei_analizar(
  vcf,
  sexo = NULL,
  fenotipo = NULL,
  categoria = NULL,
  muestras = NULL,
  evidencia = "verde",
  carpeta = ".",
  responsable = NULL,
  ...
)
```

## Arguments

- vcf:

  Ruta al archivo `.vcf` o `.vcf.gz`.

- sexo:

  Sexo de las muestras, ver
  [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md).

- fenotipo:

  Terminos HPO del paciente, ver
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).

- categoria:

  Categorias IUIS sospechadas, ver
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).

- muestras:

  Muestras a analizar. `NULL` analiza todas.

- evidencia:

  Niveles de evidencia de PanelApp, ver
  [`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md).

- carpeta:

  Carpeta donde guardar los reportes. `NULL` no genera reportes.

- responsable:

  Nombre que aparece en los reportes.

- ...:

  Otros argumentos para
  [`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md),
  por ejemplo `min_dp`.

## Value

La tabla priorizada, como
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).
Si se generaron reportes, sus rutas quedan en el atributo `reportes`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- iei_analizar(iei_ejemplo_vcf(), sexo = "M",
                    fenotipo = c("HP:0004313", "HP:0002719"),
                    carpeta = "reportes")
} # }
```
