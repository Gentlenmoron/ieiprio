# Verificar el sexo de las muestras desde el cromosoma X

Estima el sexo de cada muestra con la proporcion de variantes
heterocigotas en el cromosoma X fuera de las regiones pseudoautosomicas.
Un varon casi no tiene heterocigotos ahi. Lo compara con el sexo
declarado en
[`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md),
porque un error en ese dato cambia la interpretacion de todos los genes
ligados al X.

## Usage

``` r
iei_verificar_sexo(vcf, minimo = 20)
```

## Arguments

- vcf:

  Tabla devuelta por
  [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md),
  antes de filtrar.

- minimo:

  Numero minimo de variantes en el X para estimar.

## Value

Un tibble con `muestra`, `sexo_declarado`, `sexo_inferido` (`M`, `F`,
`indeterminado` o `NA` si hay pocas variantes), `frac_het_x`,
`n_variantes_x` y `concordante`.
