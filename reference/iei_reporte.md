# Generar el reporte HTML de un paciente

Crea un archivo HTML autocontenido, sin dependencias externas, con las
variantes priorizadas de una muestra, una ficha por cada variante de
prioridad alta, la bitacora del filtrado, las versiones de todas las
fuentes de datos y las limitaciones del analisis. Se puede abrir en
cualquier navegador y enviar por correo.

## Usage

``` r
iei_reporte(
  x,
  archivo,
  paciente = NULL,
  responsable = NULL,
  panel = ieiprio::iuis_panel,
  dominios = TRUE,
  cache = tools::R_user_dir("ieiprio", "cache")
)
```

## Arguments

- x:

  Tabla devuelta por
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).

- archivo:

  Ruta del HTML a crear.

- paciente:

  Codigo de la muestra. Obligatorio si `x` tiene varias.

- responsable:

  Nombre de quien firma el analisis (opcional).

- panel:

  Tabla de genes, por defecto
  [iuis_panel](https://gentlenmoron.github.io/ieiprio/reference/iuis_panel.md).

- dominios:

  Si `TRUE`, consulta en Ensembl el largo y los dominios de cada
  proteina para dibujar el diagrama de la variante. Usa cache.

- cache:

  Carpeta de la cache, o `FALSE` para no usarla.

## Value

La ruta del archivo, de forma invisible.

## Details

El reporte usa solo el codigo de la muestra. No incluyas nombres ni
datos que identifiquen al paciente.

## Examples

``` r
if (FALSE) { # \dontrun{
iei_reporte(prior, "P01.html", paciente = "P01")
} # }
```
