# Exportar resultados

Guarda una tabla de ieiprio (de cualquier etapa: lectura, anotacion o
priorizacion) en el formato que necesite el siguiente paso.

## Usage

``` r
iei_exportar(x, archivo, formato = NULL, muestra = NULL)
```

## Arguments

- x:

  Tabla de ieiprio.

- archivo:

  Ruta del archivo a escribir.

- formato:

  `"tsv"`, `"xlsx"`, `"vcf"`, `"bed"` o `"phenopacket"`. Por defecto se
  deduce de la extension (`.json` es phenopacket).

- muestra:

  Para `"phenopacket"`, la muestra a exportar si hay mas de una.

## Value

La ruta del archivo, invisible.

## Details

- `"tsv"`: texto separado por tabuladores, para Excel, R o Python.

- `"xlsx"`: libro de Excel con una hoja de variantes, una de parametros
  y una de fuentes de datos (necesita el paquete writexl).

- `"vcf"`: VCF 4.2 con una columna por muestra y la anotacion de ieiprio
  en INFO (gen, consecuencia, HGVS, puntaje, prioridad y ACMG), para
  abrir en IGV o pasar a otro programa.

- `"bed"`: coordenadas en formato BED (base 0), para IGV o bedtools.

- `"phenopacket"`: un Phenopacket de GA4GH (version 2) en JSON con el
  fenotipo del caso, los signos ausentes y las variantes de prioridad
  alta y media como interpretaciones genomicas. Es el formato que leen
  Exomiser, LIRICAL y los registros de enfermedades raras.

## Examples

``` r
if (FALSE) { # \dontrun{
iei_exportar(prior, "resultado.xlsx")
iei_exportar(prior, "candidatas.vcf")
iei_exportar(prior, "P01.json", muestra = "P01")
} # }
```
