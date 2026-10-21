# Validar la busqueda por fenotipo con pacientes simulados

Mide que tan bien
[`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md)
encuentra el gen correcto. Para cada caso simulado elige al azar una
enfermedad del panel, toma algunos de sus signos (los frecuentes con mas
probabilidad), cambia una parte por un termino mas general, como haria
un medico que describe "infecciones recurrentes" en lugar de "neumonia
por Pneumocystis", y agrega signos que no tienen que ver. Luego anota en
que puesto queda el gen verdadero con cada metodo.

## Usage

``` r
iei_validar_fenotipo(
  n = 100,
  terminos = 5,
  imprecision = 0.3,
  ruido = 1,
  metodos = c("enfermedad", "semantica", "exacta"),
  semilla = 1
)
```

## Arguments

- n:

  Numero de casos simulados.

- terminos:

  Signos tomados de la enfermedad en cada caso.

- imprecision:

  Proporcion de signos que se reemplazan por un termino padre.

- ruido:

  Signos al azar agregados a cada caso.

- metodos:

  Metodos a comparar: `"enfermedad"`, `"semantica"` y `"exacta"`.

- semilla:

  Semilla para que la simulacion se pueda repetir.

## Value

Un tibble de clase `iei_validacion` con una fila por caso y metodo:
`caso`, `gen`, `enfermedad`, `metodo` y `puesto` del gen verdadero (1 es
el primero). [`print()`](https://rdrr.io/r/base/print.html) resume el
porcentaje de casos en el primer puesto, entre los 5 y entre los 10
primeros;
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
grafica la curva acumulada.

## Details

Es una validacion interna: los casos salen de las mismas anotaciones que
usa el puntaje, asi que el resultado es optimista y sirve para comparar
metodos y elegir parametros, no como estimacion del rendimiento con
pacientes reales. Para eso hace falta una cohorte con diagnostico
conocido.

## Examples

``` r
v <- iei_validar_fenotipo(n = 20)
#> Simulando casos ■■■■■■■■■■■■■■■                   45% | ETA:  1s
#> Simulando casos ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  100% | ETA:  0s
v
#> Validacion de la busqueda por fenotipo
#> 20 casos simulados, 5 signos por caso, 30% imprecisos, 1 al azar; 583 genes en
#> competencia.
#> # A tibble: 3 × 5
#>   metodo     primero  top5 top10 mediana
#>   <chr>        <dbl> <dbl> <dbl>   <dbl>
#> 1 enfermedad      65    90    90       1
#> 2 semantica       60    85    95       1
#> 3 exacta          45    75    85       2
#> Porcentaje de casos con el gen correcto en el primer puesto, entre los 5 y los
#> 10 primeros. Validacion interna, optimista.
```
