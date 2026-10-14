# Leer un cromatograma SCF

Lee un archivo SCF (versiones 2 y 3), el formato abierto que exportan
muchos secuenciadores antiguos y programas como Staden, Chromas o
4Peaks. Devuelve el mismo objeto que
[`iei_leer_ab1()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_ab1.md),
asi que sirve para
[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
y para graficar con
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Usage

``` r
iei_leer_scf(archivo)
```

## Arguments

- archivo:

  Ruta al archivo `.scf`.

## Value

Un objeto `iei_ab1` con `bases`, `calidad`, `picos`, `trazas`, `muestra`
y `archivo`.

## Examples

``` r
if (FALSE) { # \dontrun{
scf <- iei_leer_scf("BTK_exon15_F.scf")
plot(scf)
} # }
```
