# Tabla simple de un resultado

Devuelve el resultado como un tibble comun, sin la clase ni los
metadatos de ieiprio, para exportarlo o trabajarlo con otros paquetes.

## Usage

``` r
iei_tabla(x, columnas = c("principales", "todas"))
```

## Arguments

- x:

  Resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  o
  [`iei_analizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_analizar.md).

- columnas:

  `"principales"` para las columnas mas usadas o `"todas"`.

## Value

Un tibble.
