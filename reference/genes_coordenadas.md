# Coordenadas de los genes del panel

Ubicacion de cada gen en GRCh37 (Ensembl 82) y GRCh38 (Ensembl 90) segun
PanelApp. La usa
[`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md)
para quedarse con variantes dentro de genes del panel.

## Usage

``` r
genes_coordenadas
```

## Format

Un tibble con columnas `gen`, `build`, `chr` (sin prefijo chr),
`inicio`, `fin` y `ensembl_id`.

## Source

<https://panelapp.genomicsengland.co.uk/panels/398/>
