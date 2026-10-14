# Graficar un resultado

Devuelve un grafico de ggplot2 que puedes modificar como cualquier otro
(agregar temas, cambiar colores, combinar con patchwork). Con
`interactivo = TRUE` lo convierte con plotly.

## Usage

``` r
# S3 method for class 'iei_resultado'
plot(
  x,
  tipo = c("puntaje", "categorias", "genes", "predictores", "mapa", "oncoplot",
    "cariotipo", "circos"),
  n = 20,
  interactivo = FALSE,
  ...
)
```

## Arguments

- x:

  Resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  o
  [`iei_analizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_analizar.md).

- tipo:

  Que graficar.

  - `"puntaje"` el desglose del puntaje de cada variante por componente.

  - `"categorias"` cuantas variantes hay en cada prioridad por muestra.

  - `"genes"` los genes con mayor puntaje.

  - `"predictores"` REVEL contra AlphaMissense, con los umbrales usados.

  - `"mapa"` el mejor puntaje de cada gen en cada muestra, para comparar
    varias muestras de una cohorte.

  - `"oncoplot"` el tipo de variante de cada gen en cada muestra, al
    estilo de maftools, con el total por gen al costado.

  - `"cariotipo"` las variantes sobre los cromosomas
    ([`iei_cariograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cariograma.md)).

  - `"circos"` las variantes sobre el genoma en circulo
    ([`iei_circos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_circos.md));
    dibuja con graficos base y no devuelve un ggplot.

- n:

  Numero maximo de variantes o genes a mostrar.

- interactivo:

  Si `TRUE`, devuelve un grafico de plotly.

- ...:

  No se usa.

## Value

Un objeto ggplot, o plotly si `interactivo = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- iei_analizar(iei_ejemplo_vcf(), sexo = "M", carpeta = NULL)
plot(res)
plot(res, tipo = "cariotipo")
plot(res, tipo = "predictores", interactivo = TRUE)
} # }
```
