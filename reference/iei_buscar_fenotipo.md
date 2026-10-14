# Buscar genes candidatos por fenotipo

Ordena los genes del panel segun cuanto se parece el fenotipo del caso
al fenotipo conocido de cada gen. Sirve antes de tener variantes, para
saber que genes mirar primero, o para elegir el panel de una
secuenciacion dirigida.

## Usage

``` r
iei_buscar_fenotipo(
  fenotipo,
  n = 20,
  evidencia = "verde",
  herencia = NULL,
  similitud = c("enfermedad", "semantica", "exacta"),
  ausentes = NULL,
  p_valor = FALSE,
  permutaciones = 200
)
```

## Arguments

- fenotipo:

  Codigos HPO o texto en espanol o ingles, como en
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md).

- n:

  Numero de genes a devolver.

- evidencia:

  Niveles de evidencia de PanelApp a considerar.

- herencia:

  Opcional, limitar a uno o varios modos de herencia (`"AR"`, `"AD"`,
  `"XL"`, `"MT"`).

- similitud:

  `"enfermedad"`, `"semantica"` o `"exacta"`.

- ausentes:

  Signos que el caso no tiene y que ayudan a descartar, por ejemplo
  `"eccema"` cuando se piensa en Wiskott-Aldrich. Mismo formato que
  `fenotipo`.

- p_valor:

  Si `TRUE`, agrega un valor p empirico: la proporcion de fenotipos al
  azar del mismo tamano que dan una similitud igual o mayor con ese gen.
  Ayuda a ver si un gen sale primero por el caso o porque su fenotipo es
  tan amplio que encaja con casi todo.

- permutaciones:

  Fenotipos al azar para el valor p.

## Value

Un tibble de clase `iei_busqueda` con `gen`, `similitud` (0 a 1),
`especificidad` (0 a 1, cuanto del fenotipo conocido del gen se parece
al del caso; desempata), `exactos` (terminos identicos), `terminos_gen`,
`enfermedad` (la mas parecida, con su identificador OMIM u Orphanet),
`herencia`, `categorias_iuis`, `p_valor` si se pidio y `coincidencias`
(para cada termino del caso, el termino mas parecido de la enfermedad o
del gen). Se grafica con
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) o
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html).

## Details

Con `similitud = "enfermedad"` (por defecto) el caso se compara con cada
enfermedad anotada en HPO por separado, porque un mismo gen puede causar
cuadros muy distintos (STAT1 con ganancia o con perdida de funcion, por
ejemplo). Para cada signo del caso se busca el signo mas parecido de la
enfermedad usando la jerarquia de HPO, y ese parecido se pesa por la
frecuencia del signo en la enfermedad: un signo constante cuenta entero
y uno que aparece en el 5% de los pacientes cuenta la mitad. Si la
enfermedad excluye un signo que el caso tiene (anotacion NOT de HPO),
ese signo resta. Los `ausentes` restan segun lo frecuentes que son en la
enfermedad. El valor del gen es el de su enfermedad mas parecida, que se
muestra en la columna `enfermedad`. Los genes sin anotaciones por
enfermedad se comparan con sus terminos de gen.

Con `"semantica"` se compara con todos los terminos del gen juntos, como
en la version 0.3, y con `"exacta"` solo cuentan los terminos identicos.
[`iei_validar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_validar_fenotipo.md)
compara los tres metodos con casos simulados.

## Examples

``` r
iei_buscar_fenotipo(c("agammaglobulinemia", "neumonia recurrente"))
#> ℹ Fenotipo interpretado como:
#> • agammaglobulinemia → HP:0004432 Agammaglobulinemia
#> • neumonia recurrente → HP:0006532 Recurrent pneumonia
#> Genes candidatos por fenotipo (enfermedad)
#> Fenotipo: HP:0004432 Agammaglobulinemia; HP:0006532 Recurrent pneumonia
#> # A tibble: 20 × 6
#>    gen     similitud especificidad enfermedad           herencia categorias_iuis
#>    <chr>       <dbl>         <dbl> <chr>                <chr>    <chr>          
#>  1 CD79B       1             0.372 Agammaglobulinemia … AR       3              
#>  2 BLNK        1             0.359 Agammaglobulinemia … AR       3              
#>  3 IGLL1       1             0.349 Agammaglobulinemia … AR       3              
#>  4 BTK         0.948         0.308 X-linked agammaglob… XL       3              
#>  5 CD79A       0.921         0.333 Agammaglobulinemia … AR       3              
#>  6 CIITA       0.921         0.277 Bare lymphocyte syn… AR       1              
#>  7 IL2RG       0.89          0.316 Severe combined imm… XL       1              
#>  8 RAG1        0.89          0.247 Omenn syndrome (OMI… AR       1              
#>  9 RAG2        0.89          0.247 Omenn syndrome (OMI… AR       1              
#> 10 DCLRE1C     0.89          0.224 Omenn syndrome (OMI… AR       1              
#> 11 CD3E        0.797         0.438 Immunodeficiency 18… AR       1              
#> 12 CD19        0.797         0.431 Immunodeficiency, c… AR       3              
#> 13 ICOS        0.797         0.42  Immunodeficiency, c… AR       1              
#> 14 IL21R       0.797         0.402 Immunodeficiency 56… AR       NA             
#> 15 LAT         0.797         0.388 Immunodeficiency 52… AR       1              
#> 16 RFXAP       0.797         0.329 MHC class II defici… AR       1              
#> 17 MTHFD1      0.797         0.292 Combined immunodefi… AR       2              
#> 18 CR2         0.797         0.275 Immunodeficiency, c… AR       3              
#> 19 ADA         0.797         0.269 Severe combined imm… AR       1              
#> 20 PTCRA       0.797         0.247 Immunodeficiency 12… AR       NA             
#> Detalle de cada coincidencia en la columna coincidencias. Grafico: plot(x) o
#> autoplot(x).
iei_buscar_fenotipo("HP:0002728", herencia = "AD")
#> Genes candidatos por fenotipo (enfermedad)
#> Fenotipo: HP:0002728 Recurrent mucocutaneous candidiasis
#> # A tibble: 20 × 6
#>    gen     similitud especificidad enfermedad           herencia categorias_iuis
#>    <chr>       <dbl>         <dbl> <chr>                <chr>    <chr>          
#>  1 IKBKB       1             0.385 Immunodeficiency 15… AR/AD    1              
#>  2 STAT1       1             0.206 Immunodeficiency 31… AR/AD    6              
#>  3 AIRE        1             0.146 Autoimmune polyendo… AR/AD    NA             
#>  4 ITPR3       1             0.139 Immunodeficiency 13… AR/AD    1              
#>  5 IKZF2       1             0.134 Immunodysregulation… AD       1              
#>  6 GATA2       0.859         0.133 Immunodeficiency 21… AD       5              
#>  7 STAT3       0.859         0.113 Hyper-IgE recurrent… AD       2;4            
#>  8 CARD11      0.828         0.295 Card11 immunodefici… AR/AD    1;2;3          
#>  9 PIK3CD      0.828         0.173 Immunodeficiency 14… AR/AD    3              
#> 10 FOXN1       0.789         0.293 Immune defect due t… AR/AD    2              
#> 11 OTULIN      0.789         0.157 Immunodeficiency 10… AR/AD    6;7            
#> 12 TRAF3       0.789         0.14  Immunodeficiency 13… AD       4;7            
#> 13 SEC61A1     0.789         0.134 Tubulointerstitial … AD       3              
#> 14 TOP2B       0.789         0.106 B-cell immunodefici… AD       3              
#> 15 NFKBIA      0.772         0.226 Hypohidrotic ectode… AD       2              
#> 16 SPI1        0.748         0.241 Autosomal non-syndr… AD       3              
#> 17 TCF3        0.748         0.208 Autosomal non-syndr… AR/AD    3              
#> 18 PIK3R1      0.748         0.123 Autosomal non-syndr… AR/AD    3              
#> 19 NFE2L2      0.691         0.097 Immunodeficiency, d… AD       2              
#> 20 GFI1        0.61          0.209 Autosomal dominant … AD       5              
#> Detalle de cada coincidencia en la columna coincidencias. Grafico: plot(x) o
#> autoplot(x).
iei_buscar_fenotipo(c("trombocitopenia", "infecciones recurrentes"), ausentes = "eccema")
#> ℹ Fenotipo interpretado como:
#> • trombocitopenia → HP:0001873 Thrombocytopenia
#> • infecciones recurrentes → HP:0002719 Recurrent infections
#> ℹ Fenotipo interpretado como:
#> • eccema → HP:0000964 Eczematoid dermatitis
#> Genes candidatos por fenotipo (enfermedad)
#> Fenotipo: HP:0001873 Thrombocytopenia; HP:0002719 Recurrent infections
#> Ausentes: HP:0000964 Eczematoid dermatitis
#> # A tibble: 20 × 6
#>    gen     similitud especificidad enfermedad           herencia categorias_iuis
#>    <chr>       <dbl>         <dbl> <chr>                <chr>    <chr>          
#>  1 VPS45           1         0.325 Neutropenia, severe… AR       5              
#>  2 IRF8            1         0.322 Immunodeficiency 32… AR/AD    6              
#>  3 ITK             1         0.283 Lymphoproliferative… AR       4              
#>  4 LRBA            1         0.256 Immunodeficiency, c… AR       4              
#>  5 LCK             1         0.253 Immunodeficiency 22… AR       1              
#>  6 SH2D1A          1         0.253 Lymphoproliferative… XL       4              
#>  7 MAGT1           1         0.244 Immunodeficiency, X… XL       4              
#>  8 XIAP            1         0.231 Lymphoproliferative… XL       4              
#>  9 ADA             1         0.197 Severe combined imm… AR       1              
#> 10 CD81            1         0.193 Immunodeficiency, c… AR       3              
#> 11 CD40LG          1         0.188 Immunodeficiency, X… XL       1              
#> 12 PNP             1         0.186 Immunodeficiency du… AR       2              
#> 13 STXBP2          1         0.174 Hemophagocytic lymp… AR       4              
#> 14 G6PC3           1         0.163 Neutropenia, severe… AR       5              
#> 15 NHEJ1           1         0.163 Immunodeficiency 12… AR       1              
#> 16 SLC46A1         1         0.159 Folate malabsorptio… AR       2              
#> 17 SEC61A1         1         0.147 Neutropenia, severe… AD       3              
#> 18 STIM1           1         0.145 Immunodeficiency 10… AR       2              
#> 19 AP3B1           1         0.139 Hermansky-Pudlak sy… AR       4              
#> 20 LYST            1         0.136 Chediak-Higashi syn… AR       4              
#> Detalle de cada coincidencia en la columna coincidencias. Grafico: plot(x) o
#> autoplot(x).
```
