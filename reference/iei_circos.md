# Mapa circular del panel en el genoma

Dibuja los cromosomas en circulo con los genes del panel coloreados por
categoria IUIS y, si se le pasa un resultado o una busqueda, marca
encima las variantes o los genes candidatos. Esta inspirado en los
graficos de circlize y sirve para ver de un vistazo donde estan los
genes del panel, en que categoria cae cada candidato y como se reparten
las variantes de un caso o de una cohorte.

## Usage

``` r
iei_circos(
  x = NULL,
  build = c("GRCh38", "GRCh37"),
  etiquetas = 15,
  titulo = NULL
)
```

## Arguments

- x:

  `NULL`, un resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  o una busqueda de
  [`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md).

- build:

  Version del genoma.

- etiquetas:

  Cuantos genes rotular.

- titulo:

  Titulo del grafico; `NULL` usa uno por defecto.

## Value

No devuelve nada; dibuja en el dispositivo grafico activo.

## Details

- Sin `x`, muestra el panel completo con la densidad de genes por
  region.

- Con un resultado de
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md),
  cada variante es un punto cuya distancia al centro es su puntaje y
  cuyo color es su prioridad. Se rotulan los genes de las mejores
  variantes.

- Con una busqueda de
  [`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md),
  cada gen candidato es una barra tan alta como su similitud, y los
  candidatos que comparten categoria IUIS se unen por dentro.

## See also

[`iei_cariograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cariograma.md)
para la misma informacion en cromosomas lineales.

## Examples

``` r
# \donttest{
if (requireNamespace("circlize", quietly = TRUE)) {
  iei_circos()
  b <- iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente"), n = 15)
  iei_circos(b)
}

#> ℹ Fenotipo interpretado como:
#> • agammaglobulinemia → HP:0004432 Agammaglobulinemia
#> • neumonia recurrente → HP:0006532 Recurrent pneumonia

# }
```
