# Ficha de un gen

Reune todo lo que se sabe de un gen del panel: herencia, evidencia en
PanelApp, categorias IUIS (con la fuente de cada una), mecanismo de
enfermedad, alias, ubicacion, fenotipos de PanelApp y terminos HPO. Si
hay conexion agrega las variantes conocidas en ClinVar y la proteina con
sus dominios Pfam, para ver donde se concentran las variantes
patogenicas.

## Usage

``` r
iei_gen(
  gen,
  clinvar = TRUE,
  proteina = TRUE,
  build = c("GRCh38", "GRCh37"),
  cache = tools::R_user_dir("ieiprio", "cache")
)
```

## Arguments

- gen:

  Simbolo del gen o un alias (`"BTK"`, `"TACI"`).

- clinvar:

  Si `TRUE`, consulta las variantes de ClinVar.

- proteina:

  Si `TRUE`, consulta la proteina y sus dominios en Ensembl.

- build:

  `"GRCh38"` o `"GRCh37"`.

- cache:

  Carpeta de la cache, o `FALSE` para no usarla.

## Value

Una lista de clase `iei_gen`. Se imprime como ficha y se grafica con
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) o
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html);
ver
[`autoplot.iei_gen()`](https://gentlenmoron.github.io/ieiprio/reference/autoplot.iei.md).

## Details

Las consultas en linea se guardan en la cache por dia, asi que repetir
la misma ficha no vuelve a consultar. Sin conexion la ficha sale igual
con la informacion incluida en el paquete.

## Examples

``` r
g <- iei_gen("BTK", clinvar = FALSE, proteina = FALSE)
g
#> 
#> ── BTK ──
#> 
#> • Herencia: XL ─ evidencia en PanelApp: verde
#> • Mecanismo: pérdida de función o no especificado
#> • Categoría IUIS: 3 Deficiencias predominantemente de anticuerpos (PanelApp,
#> tabla IUIS)
#> • Alias: ATK, PSCTK1, XLA
#> • Ubicación: X:101,349,447-101,390,796 (GRCh38)
#> • Fenotipos en PanelApp: Agammaglobulinemia and isolated hormone deficiency;
#> Agammaglobulinemia, X-linked; Agammaglobulinemia, X-linked 1 (y 6 mas)
#> • HPO: 77 términos; los más específicos: Complete or near-complete absence of
#> specific antibody response to protein-conjugated Haemophilus influenzae type b
#> vaccine, Enteroviral dermatomyositis syndrome, Enteroviral hepatitis,
#> Bronchiolitis obliterans, Lymph node hypoplasia, Recurrent enteroviral
#> infections
#> Grafico: plot(x) o autoplot(x); variantes de ClinVar en x$clinvar$variantes.
if (FALSE) { # \dontrun{
g <- iei_gen("BTK")
g$clinvar$variantes
plot(g)
} # }
```
