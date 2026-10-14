# Consultar una variante al vuelo

Atajo para revisar una o pocas variantes sin VCF: las ubica, las anota
con Ensembl VEP y las prioriza. Imprime una ficha por variante.

## Usage

``` r
iei_variante(
  x,
  genotipo = NA,
  sexo = NULL,
  fenotipo = NULL,
  build = c("GRCh38", "GRCh37")
)
```

## Arguments

- x:

  Un vector con variantes (`"BTK c.1632-2A>G"`, `"rs699"`,
  `"NM_000061.3:c.1632-2A>G"`), un texto libre, o la ruta a un archivo
  de texto, PDF o Word.

- genotipo:

  Genotipo de las variantes (`"het"`, `"hom_alt"`, `"hemi"`). Puede ser
  un valor o uno por variante. Si es `NA` se usa el que diga el texto, y
  si no dice nada queda sin indicar.

- sexo:

  `"M"` o `"F"`. En un varon, una variante homocigota en el X fuera de
  la region pseudoautosomica pasa a ser hemicigota.

- fenotipo:

  Terminos HPO del paciente, opcional.

- build:

  `"GRCh38"` o `"GRCh37"`.

## Value

Un `iei_resultado`, de forma invisible.

## Examples

``` r
if (FALSE) { # \dontrun{
iei_variante("BTK c.1632-2A>G", genotipo = "hemi")
} # }
```
