# Panel de genes de errores innatos de la inmunidad

Genes del panel 398 de PanelApp (Genomics England), "Primary
immunodeficiency or monogenic inflammatory bowel disease", en la version
indicada en
[`iei_fuentes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_fuentes.md),
unidos a una curacion propia de categorias IUIS. Se genera con
`data-raw/construir_datos.R`.

## Usage

``` r
iuis_panel
```

## Format

Un tibble con una fila por gen:

- gen:

  Simbolo HGNC.

- hgnc_id:

  Identificador HGNC.

- herencia:

  AD, AR, AR/AD, XL, MT o NA.

- herencia_panelapp:

  Texto original de PanelApp.

- evidencia:

  verde, ambar o rojo segun PanelApp.

- fenotipos:

  Fenotipos listados en PanelApp, separados por punto y coma.

- categorias_iuis:

  Categorias IUIS curadas, separadas por punto y coma. NA si no estan
  curadas.

- ganancia_funcion:

  `TRUE` si el gen tiene una enfermedad por ganancia de funcion.

## Source

<https://panelapp.genomicsengland.co.uk/panels/398/>
