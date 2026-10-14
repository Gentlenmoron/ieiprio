# Leer el texto de un archivo

Extrae el texto de un informe para buscar variantes en el. Acepta texto
plano (`.txt`, `.md`, `.csv`, `.tsv`), PDF con texto (no escaneado) y
Word (`.docx`).

## Usage

``` r
iei_leer_texto(archivo)
```

## Arguments

- archivo:

  Ruta al archivo.

## Value

Un vector de texto, una linea por elemento.

## Examples

``` r
if (FALSE) { # \dontrun{
texto <- iei_leer_texto("informe_secuenciacion.pdf")
iei_extraer_variantes(texto)
} # }
```
