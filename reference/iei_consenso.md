# Secuencia consenso de las lecturas Sanger

Junta las lecturas alineadas por
[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
en una sola secuencia en coordenadas del genoma. En cada posicion gana
la base con mas calidad sumada entre las lecturas; si dos lecturas
buenas no coinciden, la posicion se marca como discordante y se escribe
con el codigo IUPAC de las dos bases. Las inserciones respecto de la
referencia no entran en el consenso; estan en la tabla de variantes.

## Usage

``` r
iei_consenso(x, archivo = NULL, min_calidad = 20)
```

## Arguments

- x:

  Resultado de
  [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
  (o de
  [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
  sobre el).

- archivo:

  Opcional, ruta de un FASTA donde guardar el consenso.

- min_calidad:

  Calidad Phred minima para que una lectura cuente como discordante.

## Value

Un tibble con `pos`, `ref`, `consenso`, `lecturas` (cuantas cubren la
posicion) y `discordante`. El atributo `secuencia` tiene el consenso
como texto.

## Examples

``` r
if (FALSE) { # \dontrun{
v <- iei_sanger(c("F.ab1", "R.ab1"), gen = "BTK")
cons <- iei_consenso(v, "BTK_consenso.fasta")
} # }
```
