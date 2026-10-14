# Crear una tabla de variantes desde texto

Convierte descripciones de variantes en una tabla con la misma
estructura que
[`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md),
lista para
[`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
e
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).
Sirve cuando no hay VCF, por ejemplo con un informe de otro centro de
secuenciacion o una lista de variantes escrita a mano. Las coordenadas
se obtienen del Variant Recoder de Ensembl.

## Usage

``` r
iei_variantes(
  x,
  muestra = "consulta",
  genotipo = NA,
  sexo = NULL,
  build = c("GRCh38", "GRCh37")
)
```

## Arguments

- x:

  Un vector con variantes (`"BTK c.1632-2A>G"`, `"rs699"`,
  `"NM_000061.3:c.1632-2A>G"`), un texto libre, o la ruta a un archivo
  de texto, PDF o Word.

- muestra:

  Codigo de la muestra.

- genotipo:

  Genotipo de las variantes (`"het"`, `"hom_alt"`, `"hemi"`). Puede ser
  un valor o uno por variante. Si es `NA` se usa el que diga el texto, y
  si no dice nada queda sin indicar.

- sexo:

  `"M"` o `"F"`. En un varon, una variante homocigota en el X fuera de
  la region pseudoautosomica pasa a ser hemicigota.

- build:

  `"GRCh38"` o `"GRCh37"`.

## Value

Una tabla como la de
[`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md)
con las columnas extra `gen` (si cae en un gen del panel), `entrada` y
`ambigua` (una descripcion de proteina que corresponde a varias
posiciones posibles).

## Examples

``` r
if (FALSE) { # \dontrun{
v <- iei_variantes(c("BTK c.1632-2A>G", "RAG1 p.Arg396Cys"), sexo = "M")
iei_priorizar(iei_anotar(v), fenotipo = "HP:0004313")

v <- iei_variantes("informe.pdf")
} # }
```
