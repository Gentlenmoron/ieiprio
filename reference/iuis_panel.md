# Panel de genes de errores innatos de la inmunidad

Genes del panel 398 de PanelApp (Genomics England), "Primary
immunodeficiency or monogenic inflammatory bowel disease", en la version
indicada en
[`iei_fuentes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_fuentes.md),
con su categoria IUIS y el mecanismo de enfermedad. Se genera con
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

  Categorias IUIS separadas por punto y coma, ver
  [iuis_categorias](https://gentlenmoron.github.io/ieiprio/reference/iuis_categorias.md).
  NA si no tiene.

- mecanismo:

  `"GOF"` si su enfermedad es por ganancia de funcion, `"LOF y GOF"` si
  tiene enfermedades por los dos mecanismos, NA en el resto (perdida de
  funcion o sin especificar).

- ganancia_funcion:

  `TRUE` si `mecanismo` incluye ganancia de funcion.

## Source

<https://panelapp.genomicsengland.co.uk/panels/398/>
