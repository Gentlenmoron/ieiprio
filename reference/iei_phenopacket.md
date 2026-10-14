# Phenopacket de un caso

Arma un Phenopacket de GA4GH (esquema 2.0) como lista de R: el sujeto,
el fenotipo usado en la priorizacion (con los signos ausentes marcados
como `excluded`) y las variantes de prioridad alta y media como
interpretaciones genomicas candidatas, con su clase ACMG sugerida. El
diagnostico propuesto es la enfermedad del gen de la primera variante
que mas se parece al fenotipo. Se guarda como JSON con
[`iei_exportar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_exportar.md).

## Usage

``` r
iei_phenopacket(x, muestra = NULL)
```

## Arguments

- x:

  Resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).

- muestra:

  Muestra a exportar; obligatoria si hay mas de una.

## Value

Una lista con la estructura del Phenopacket.

## Examples

``` r
if (FALSE) { # \dontrun{
pp <- iei_phenopacket(prior, muestra = "P01")
} # }
```
