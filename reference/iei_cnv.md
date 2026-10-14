# CNV y variantes estructurales en genes del panel

Devuelve los CNV y variantes estructurales que traia el VCF (registros
con ALT simbolico como `<DEL>` o `<DUP>`, o con `SVTYPE` en INFO), con
los genes del panel que tocan. ieiprio no llama CNV; esta funcion
muestra los que ya llamo otro programa (por ejemplo un llamador de CNV
de exoma), para no perder las deleciones grandes de genes como BTK, CYBB
o la region 22q11.2.

## Usage

``` r
iei_cnv(x, genes = NULL)
```

## Arguments

- x:

  Tabla de
  [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md)
  o cualquier tabla derivada que conserve sus atributos.

- genes:

  Genes a cruzar. Por defecto, todo el panel.

## Value

Un tibble con `muestra`, `chr`, `inicio`, `fin`, `tamano` (pb), `tipo`,
`cn` (numero de copias si el VCF lo trae), `genes` del panel que toca y
`nota`.

## Examples

``` r
if (FALSE) { # \dontrun{
iei_cnv(iei_leer_vcf("paciente_con_cnv.vcf.gz"))
} # }
```
