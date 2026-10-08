# Filtrar por frecuencia poblacional segun el modo de herencia

Descarta variantes demasiado frecuentes en gnomAD para causar una
enfermedad rara. El umbral depende de la herencia del gen, porque una
variante recesiva puede ser mas frecuente que una dominante. Las
variantes ausentes de gnomAD siempre pasan.

## Usage

``` r
iei_filtrar_frecuencia(
  x,
  max_af = c(AD = 1e-04, XL = 1e-04, AR = 0.01, `AR/AD` = 0.01, MT = 0.01)
)
```

## Arguments

- x:

  Tabla devuelta por
  [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md),
  con columna `gen`.

- max_af:

  Umbrales por modo de herencia. Los genes sin herencia conocida usan el
  de `AR`.

## Value

La tabla filtrada, con el paso "Frecuencia" agregado a la bitacora.
