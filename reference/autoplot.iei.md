# Graficos con autoplot

Metodos de
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
para los objetos de ieiprio. Devuelven un grafico de ggplot2 que puedes
modificar (temas, colores, facetas) o combinar con otros.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) hace lo mismo.

## Usage

``` r
# S3 method for class 'iei_resultado'
autoplot(
  object,
  tipo = c("puntaje", "categorias", "genes", "predictores", "mapa"),
  n = 20,
  ...
)

# S3 method for class 'iei_gen'
autoplot(object, tipo = c("proteina", "dominios", "clinvar"), ...)

# S3 method for class 'iei_busqueda'
autoplot(object, tipo = c("ranking", "mapa", "cariotipo"), n = 20, ...)

# S3 method for class 'iei_validacion'
autoplot(object, hasta = 20, ...)

# S3 method for class 'iei_ab1'
autoplot(object, desde = 1, hasta = min(desde + 59, nchar(object$bases)), ...)
```

## Arguments

- object:

  Objeto de ieiprio.

- tipo:

  Que graficar; depende del objeto.

- n:

  Numero maximo de elementos a mostrar.

- ...:

  Otros argumentos; en resultados se pasan a
  [`plot.iei_resultado()`](https://gentlenmoron.github.io/ieiprio/reference/plot.iei_resultado.md).

- desde, hasta:

  Bases del cromatograma a mostrar; en una validacion, `hasta` es el
  ultimo puesto de la curva.

## Value

Un objeto ggplot.

## Details

- Un resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  acepta los tipos de
  [`plot.iei_resultado()`](https://gentlenmoron.github.io/ieiprio/reference/plot.iei_resultado.md):
  `"puntaje"`, `"categorias"`, `"genes"`, `"predictores"` y `"mapa"`
  (genes por muestra, util para cohortes).

- Una ficha de
  [`iei_gen()`](https://gentlenmoron.github.io/ieiprio/reference/iei_gen.md)
  muestra las variantes patogenicas de ClinVar sobre la proteina con sus
  dominios (`"proteina"`), cuantas variantes acumula cada dominio frente
  al promedio (`"dominios"`) o el reparto de todas las variantes por
  clasificacion (`"clinvar"`).

- Una busqueda de
  [`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md)
  muestra el ranking de genes (`"ranking"`), la similitud de cada
  termino del caso con cada gen (`"mapa"`) o los candidatos sobre los
  cromosomas (`"cariotipo"`). Con
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) acepta ademas
  `"acordes"`, un diagrama de cuerdas entre los signos del caso y los
  genes, y `"circos"`, el genoma en circulo; ambos usan circlize y
  dibujan con graficos base.

- Un cromatograma de
  [`iei_leer_ab1()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_ab1.md)
  muestra las trazas de un tramo.

- Una validacion de
  [`iei_validar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_validar_fenotipo.md)
  muestra, para cada metodo, el porcentaje de casos con el gen correcto
  entre los primeros puestos.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
autoplot(iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente")))
autoplot(iei_gen("BTK"))
autoplot(prior, tipo = "mapa")
} # }
```
