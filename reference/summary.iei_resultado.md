# Resumen de un resultado

Resumen de un resultado

## Usage

``` r
# S3 method for class 'iei_resultado'
summary(object, ...)
```

## Arguments

- object:

  Resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  o
  [`iei_analizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_analizar.md).

- ...:

  No se usa.

## Value

Un objeto `summary.iei_resultado`, una lista con las tablas
`por_categoria`, `por_gen`, `por_muestra`, `acmg` y `bitacora`.
