# Graficos de un resultado de iei_sanger

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) sobre la salida
de
[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
dibuja el cromatograma de la variante
([`iei_cromatograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cromatograma.md))
o el dot plot de la lectura frente a la referencia
([`iei_dotplot()`](https://gentlenmoron.github.io/ieiprio/reference/iei_dotplot.md)).

## Usage

``` r
# S3 method for class 'iei_sanger'
plot(x, tipo = c("cromatograma", "dotplot"), ...)
```

## Arguments

- x:

  Resultado de
  [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md).

- tipo:

  `"cromatograma"` o `"dotplot"`.

- ...:

  Otros argumentos para
  [`iei_cromatograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cromatograma.md)
  o
  [`iei_dotplot()`](https://gentlenmoron.github.io/ieiprio/reference/iei_dotplot.md).

## Value

Un grafico de ggplot2.

## Examples

``` r
ab1 <- iei_ejemplo_sanger()
#> Cromatogramas de ejemplo en /tmp/RtmpEJUhwz/sanger_19cf1b19e77b.
ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
v <- iei_sanger(ab1, gen = "BTK", sexo = "M", referencia = ref)
#> BTK_c1632-2_F.ab1: hebra +, 461 bases alineadas en X:101353757-101354217,
#> identidad 99.8%.
#> BTK_c1632-2_R.ab1: hebra -, 461 bases alineadas en X:101353758-101354218,
#> identidad 99.8%.
plot(v, tipo = "dotplot")
```
