# Analizar un paciente con ieiprio

ieiprio toma un VCF de exoma o de panel y devuelve las variantes en
genes de errores innatos de la inmunidad (EII) ordenadas por prioridad,
con una explicación de cada puntaje y un reporte HTML por paciente.

> **Aviso.** Es una herramienta de apoyo e investigación. No es un
> informe diagnóstico ni reemplaza la clasificación ACMG/AMP hecha por
> un especialista.

``` r

library(ieiprio)
```

## Para qué sirve

Ante un paciente con sospecha de inmunodeficiencia primaria, el exoma o
el panel genético devuelve miles de variantes. La mayoría son normales,
comunes o están en genes que no tienen relación con el sistema inmune.
ieiprio reduce esa lista a unas pocas candidatas y las ordena con una
regla explícita.

1.  Se queda con los genes de errores innatos de la inmunidad del panel
    de PanelApp.
2.  Descarta lo que tiene mala calidad técnica.
3.  Pregunta a Ensembl qué le hace cada variante al gen y qué dicen
    gnomAD y ClinVar.
4.  Descarta lo que es demasiado frecuente para causar una enfermedad
    rara.
5.  Puntúa lo que queda según ClinVar, impacto, rareza, herencia y
    fenotipo, y escribe la razón de cada punto.
6.  Genera un reporte con control de calidad, gráficos y trazabilidad
    completa.

Cada parámetro de estos pasos se explica en detalle en la [guía de
parámetros](https://gentlenmoron.github.io/ieiprio/articles/parametros.md).

## El caso de ejemplo

El paquete trae un paciente ficticio, EJEMPLO01, varón con
hipogammaglobulinemia e infecciones recurrentes. Su VCF mezcla cinco
variantes reales tomadas de ClinVar con unas 4000 variantes de fondo
sintéticas, ubicadas fuera de los genes del panel, que solo sirven para
el control de calidad.

| Rol | Gen | Variante | ClinVar |
|----|----|----|----|
| Causal | BTK | c.1632-2A\>G | [265455](https://www.ncbi.nlm.nih.gov/clinvar/variation/265455/) |
| Portador | RAG1 | c.1186C\>T (p.Arg396Cys) | [13144](https://www.ncbi.nlm.nih.gov/clinvar/variation/13144/) |
| Incierta | STAT1 | c.\*1628A\>G | [333246](https://www.ncbi.nlm.nih.gov/clinvar/variation/333246/) |
| Benigna | BTK | c.\*334T\>G | [367692](https://www.ncbi.nlm.nih.gov/clinvar/variation/367692/) |
| Benigna | CYBB | c.1314+19C\>T | [35968](https://www.ncbi.nlm.nih.gov/clinvar/variation/35968/) |

## Todo en una línea

``` r

res <- iei_analizar(
  iei_ejemplo_vcf(),
  sexo = "M",
  fenotipo = c("HP:0004313", "HP:0002719"),
  carpeta = tempdir()
)
#> 
#> ── ieiprio ─────────────────────────────────────────────────────────────────────
#> ℹ Leyendo ejemplo.vcf.gz
#> ℹ Filtrando 4094 filas
#> ℹ Anotando 5 variantes con Ensembl VEP
#> ✔ Prioridad alta 1, media 0, baja 2
#> ✔ Reporte en /tmp/RtmpdebRA8/EJEMPLO01.html
res[, c("gen", "hgvs_c", "consecuencia", "puntaje", "categoria")]
#> # A tibble: 3 × 5
#>   gen   hgvs_c                        consecuencia            puntaje categoria
#>   <chr> <chr>                         <chr>                     <dbl> <chr>    
#> 1 BTK   ENST00000308731.8:c.1632-2A>G splice_acceptor_variant    11.5 alta     
#> 2 RAG1  ENST00000299440.6:c.1186C>T   missense_variant            4.5 baja     
#> 3 STAT1 NA                            downstream_gene_variant     3.5 baja
```

[`iei_analizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_analizar.md)
hace los mismos pasos que se explican abajo y deja un reporte HTML por
muestra en `carpeta`. El resto de esta guía los recorre uno por uno.

## Paso a paso

### 1. Leer el VCF

``` r

vcf <- iei_leer_vcf(iei_ejemplo_vcf(), sexo = "M")
vcf
#> # A tibble: 4,094 × 13
#>    muestra   chr       pos id    ref   alt   genotipo  qual    dp    gq    ab
#>    <chr>     <chr>   <int> <chr> <chr> <chr> <chr>    <dbl> <dbl> <dbl> <dbl>
#>  1 EJEMPLO01 1     1605175 NA    G     C     hom_alt    234    52    74 1    
#>  2 EJEMPLO01 1     1617029 NA    T     A     hom_alt    315    45    99 1    
#>  3 EJEMPLO01 1     2547757 NA    C     T     hom_alt    667    43    90 1    
#>  4 EJEMPLO01 1     3239502 NA    G     C     het        780    46    83 0.630
#>  5 EJEMPLO01 1     4235075 NA    A     C     het        245    56    99 0.625
#>  6 EJEMPLO01 1     4256441 NA    A     G     hom_alt    296    16    82 1    
#>  7 EJEMPLO01 1     6010805 NA    C     T     hom_alt    278    26    77 1    
#>  8 EJEMPLO01 1     6129328 NA    T     C     het        400     9    75 0.556
#>  9 EJEMPLO01 1     6429026 NA    A     G     het        695    68    74 0.5  
#> 10 EJEMPLO01 1     6891176 NA    A     G     het        891    36    83 0.611
#> # ℹ 4,084 more rows
#> # ℹ 2 more variables: filtro <chr>, build <chr>
```

El sexo importa. En un varón, una variante homocigota en el cromosoma X
fuera de las regiones pseudoautosómicas se marca como `hemi`, que es lo
que permite interpretar genes ligados al X como BTK, CYBB o WAS.

### 2. Verificar el sexo

Un sexo mal declarado cambia toda la interpretación de los genes ligados
al X. ieiprio lo estima con la proporción de heterocigotos en el X.

``` r

iei_verificar_sexo(vcf)
#> # A tibble: 1 × 6
#>   muestra   sexo_declarado sexo_inferido frac_het_x n_variantes_x concordante
#>   <chr>     <chr>          <chr>              <dbl>         <int> <lgl>      
#> 1 EJEMPLO01 M              M                      0           139 TRUE
```

### 3. Filtrar

``` r

filtrado <- iei_filtrar(vcf)
iei_bitacora(filtrado)
#> # A tibble: 4 × 3
#>   paso     filas descartadas
#>   <chr>    <int>       <int>
#> 1 Inicial   4094           0
#> 2 Genotipo  4094           0
#> 3 Panel        5        4089
#> 4 Calidad      5           0
filtrado[, c("gen", "chr", "pos", "ref", "alt", "genotipo", "dp")]
#> # A tibble: 5 × 7
#>   gen   chr         pos ref   alt   genotipo    dp
#>   <chr> <chr>     <int> <chr> <chr> <chr>    <dbl>
#> 1 STAT1 2     190969075 T     C     het         52
#> 2 RAG1  11     36574490 C     T     het         52
#> 3 CYBB  X      37805187 C     T     hemi        52
#> 4 BTK   X     101349551 A     C     hemi        52
#> 5 BTK   X     101353990 T     C     hemi        52
```

El panel por defecto son los genes con evidencia verde del panel 398 de
PanelApp. Puedes ver el panel completo con
[`iei_panel()`](https://gentlenmoron.github.io/ieiprio/reference/iei_panel.md)
o buscar un gen con `iei_buscar_gen("TACI")`.

### 4. Anotar con Ensembl VEP

``` r

anotado <- iei_anotar(filtrado)
anotado <- iei_filtrar_frecuencia(anotado)
anotado[, c("gen", "hgvs_c", "consecuencia", "af_gnomad", "clinvar")]
#> # A tibble: 3 × 5
#>   gen   hgvs_c                        consecuencia             af_gnomad clinvar
#>   <chr> <chr>                         <chr>                        <dbl> <chr>  
#> 1 STAT1 NA                            downstream_gene_variant  0.0000263 uncert…
#> 2 RAG1  ENST00000299440.6:c.1186C>T   missense_variant         0.0000376 pathog…
#> 3 BTK   ENST00000308731.8:c.1632-2A>G splice_acceptor_variant NA         pathog…
```

Las respuestas de VEP se guardan en una caché local, así que volver a
anotar al mismo paciente no consulta la API otra vez.

### 5. Priorizar

``` r

prior <- iei_priorizar(anotado, fenotipo = c("HP:0004313", "HP:0002719"))
prior[, c("gen", "puntaje", "categoria", "acmg")]
#> # A tibble: 3 × 4
#>   gen   puntaje categoria acmg                
#>   <chr>   <dbl> <chr>     <chr>               
#> 1 BTK      11.5 alta      PVS1, PM2_Supporting
#> 2 RAG1      4.5 baja      PM2_Supporting      
#> 3 STAT1     3.5 baja      PM2_Supporting
cat(prior$razones, sep = "\n\n")
#> ClinVar patogenica (+4); impacto alto (+3); ausente en gnomAD (+1); hemicigota en gen ligado al X (+2); fenotipo 1 de 2 terminos HPO (+1.5)
#> 
#> ClinVar patogenica (+4); impacto moderado (+1); heterocigoto unico en gen AR (-2); fenotipo 1 de 2 terminos HPO (+1.5)
#> 
#> compatible con herencia AD (+2); fenotipo 1 de 2 terminos HPO (+1.5)
```

Cada punto tiene una razón escrita. Para este paciente lo esperable es
esto.

- **BTK c.1632-2A\>G queda en prioridad alta.** Rompe el sitio aceptor
  de splicing, ClinVar la clasifica como patogénica, está ausente en
  gnomAD y el paciente es hemicigoto en un gen ligado al X. Es el cuadro
  de una agammaglobulinemia ligada al X.
- **RAG1 p.Arg396Cys queda en prioridad baja.** Es patogénica en
  ClinVar, pero es una sola variante en un gen recesivo. El reporte lo
  marca como portador y sugiere buscar una segunda variante si el
  fenotipo encaja.
- **STAT1 c.\*1628A\>G queda en prioridad baja.** Está en la región 3’
  no traducida y su significado es incierto.
- **Las dos benignas desaparecen** en el filtro de frecuencia o en el de
  ClinVar, y eso queda registrado en la bitácora.

Los puntos se pueden ajustar con
[`iei_pesos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_pesos.md).

### 6. Reporte

``` r

iei_reporte(prior, file.path(tempdir(), "EJEMPLO01.html"), paciente = "EJEMPLO01")
#> Reporte de "EJEMPLO01" guardado en
#> /tmp/RtmpdebRA8/EJEMPLO01.html.
```

El reporte es un único archivo HTML que se abre en cualquier navegador.
Incluye el control de calidad de la muestra, la tabla priorizada con
criterios ACMG sugeridos, fichas con el diagrama de cada proteína, la
coincidencia con el fenotipo, el embudo del filtrado y las versiones de
todas las fuentes.

## Reproducibilidad

Los datos incluidos están fijados a versiones concretas y cada archivo
tiene su huella md5.

``` r

iei_fuentes()
#> # A tibble: 3 × 6
#>   fuente        detalle                                version url   fecha md5  
#>   <chr>         <chr>                                  <chr>   <chr> <chr> <chr>
#> 1 PanelApp      Primary immunodeficiency or monogenic… 8.78    http… 2026… 6fe8…
#> 2 HPO           genes_to_phenotype.txt                 v2026-… http… 2026… 6dac…
#> 3 Curacion IUIS data-raw/iuis_curado.csv               git     http… 2026… d554…
```

Para usar una versión más nueva del panel sin tocar los datos del
paquete,
[`iei_actualizar_panel()`](https://gentlenmoron.github.io/ieiprio/reference/iei_actualizar_panel.md)
la descarga a tu caché.

## Limitaciones

- Solo analiza variantes puntuales e indels pequeños. No detecta CNV ni
  variantes estructurales.
- No evalúa cobertura, porque eso requiere el BAM.
- Los criterios ACMG son sugerencias automaticas y la fase de un posible
  heterocigoto compuesto debe confirmarse estudiando a los padres.
