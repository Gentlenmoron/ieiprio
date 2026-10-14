# Aplicacion interactiva

Abre en el navegador una aplicacion Shiny para usar ieiprio sin escribir
codigo: buscar genes candidatos por fenotipo, revisar la ficha de un
gen, interpretar una variante escrita, analizar un archivo (VCF,
cromatograma Sanger, informe en PDF o Word) con descarga del reporte y
de la tabla en Excel, y validar la busqueda por fenotipo con casos
simulados.

## Usage

``` r
iei_app(...)
```

## Arguments

- ...:

  Argumentos para
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html), por
  ejemplo `port` o `launch.browser`.

## Value

No devuelve nada; corre hasta que cierres la aplicacion.

## Details

La aplicacion corre en tu computadora. Como el resto del paquete, envia
a Ensembl, NCBI y gnomAD solo las coordenadas de las variantes; puedes
activar el modo sin conexion desde la propia aplicacion.

## Examples

``` r
if (FALSE) { # \dontrun{
iei_app()
} # }
```
