# Graficar el cromatograma alrededor de una variante

Dibuja las cuatro trazas de cada lectura `.ab1` en una ventana alrededor
de una variante encontrada por
[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md),
con la base llamada arriba y la base de referencia abajo. Las lecturas
reversas se muestran ya complementadas y en el sentido del genoma, para
compararlas directamente con las directas.

## Usage

``` r
iei_cromatograma(
  x,
  variante = 1,
  posicion = NULL,
  lectura = NULL,
  ventana = 12,
  estilo = c("clasico", "iei")
)
```

## Arguments

- x:

  Resultado de
  [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md),
  o de
  [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
  o
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  aplicados sobre el.

- variante:

  Fila de `x` a mostrar.

- posicion:

  Alternativa a `variante`: posicion genomica a mostrar, util para
  revisar una zona sin variante o la marcada como posible indel.

- lectura:

  Nombre de la lectura a mostrar. Por defecto, todas las que cubren la
  posicion.

- ventana:

  Bases a cada lado de la posicion.

- estilo:

  `"clasico"` (por defecto) dibuja las trazas sobre fondo blanco, como
  un visor de cromatogramas. `"iei"` usa el tema del resto de los
  graficos del paquete.

## Value

Un grafico de ggplot2.

## Examples

``` r
ab1 <- iei_ejemplo_sanger("het")
#> Cromatogramas de ejemplo en /tmp/RtmphYR5wI/sanger_1a083706e8f6.
ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
v <- iei_sanger(ab1, gen = "BTK", referencia = ref)
#> BTK_c1632-2_F.ab1: hebra +, 461 bases alineadas en X:101353757-101354217,
#> identidad 99.8%.
#> BTK_c1632-2_R.ab1: hebra -, 461 bases alineadas en X:101353758-101354218,
#> identidad 99.8%.
# \donttest{
iei_cromatograma(v)

iei_cromatograma(v, estilo = "iei")

# }
```
