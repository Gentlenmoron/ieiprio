# Pesos por defecto de la priorizacion

Devuelve la lista de puntos que usa
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).
Puedes modificar cualquier valor y pasar la lista en el argumento
`pesos`.

## Usage

``` r
iei_pesos()
```

## Value

Una lista con nombres.

## Examples

``` r
p <- iei_pesos()
p$clinvar_patogenica <- 5
```
