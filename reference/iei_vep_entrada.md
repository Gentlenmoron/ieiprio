# Archivo de entrada para un VEP local

Escribe las variantes en el formato por defecto de VEP, para anotarlas
en tu computadora sin enviarlas a internet. Despues corre VEP con salida
JSON y los mismos complementos que usa el paquete, y pasa el resultado a
[`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
con `vep_json`.

## Usage

``` r
iei_vep_entrada(x, archivo)
```

## Arguments

- x:

  Tabla con `chr`, `pos`, `ref` y `alt`.

- archivo:

  Ruta del archivo a escribir.

## Value

La ruta, de forma invisible.

## Details

Un comando de ejemplo (requiere VEP con su cache y los complementos
instalados):

`vep -i entrada.txt --json -o salida.json --cache --offline --assembly GRCh38 --hgvs --mane --canonical --protein --numbers --af_gnomade --af_gnomadg --check_existing --plugin AlphaMissense,file=... --plugin REVEL,file=... --plugin SpliceAI,snv=...,indel=...`
