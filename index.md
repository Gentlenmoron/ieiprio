# ieiprio

`ieiprio` encuentra, anota y prioriza variantes en genes de errores
innatos de la inmunidad a partir de secuencias Sanger, FASTA, VCF o el
texto de un informe.

![Cromatograma de BTK c.1632-2A\>G en una portadora, con el pico doble
en la posición de la variante” /\>\</p\> \</div\> \<div
id=](reference/figures/cromatograma.png)

## Instalación

``` r

install.packages("pak")
pak::pak("Gentlenmoron/ieiprio")
```

Con `remotes` también funciona.

``` r

remotes::install_github("Gentlenmoron/ieiprio")
```

## Ejemplo

``` r

library(ieiprio)
ab1 <- iei_ejemplo_sanger("het")
ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
v <- iei_sanger(ab1, gen = "BTK", sexo = "F", referencia = ref)
v[, c("gen", "chr", "pos", "ref", "alt", "genotipo", "qual", "ab", "lecturas")]
iei_cromatograma(v)
```

    #> BTK_c1632-2_F.ab1: hebra +, 461 bases alineadas en X:101353757-101354217, identidad 99.8%.
    #> BTK_c1632-2_R.ab1: hebra -, 461 bases alineadas en X:101353758-101354218, identidad 99.8%.
    #> # A tibble: 1 × 9
    #>   gen   chr         pos ref   alt   genotipo  qual    ab lecturas
    #>   <chr> <chr>     <int> <chr> <chr> <chr>    <dbl> <dbl>    <int>
    #> 1 BTK   X     101353990 T     C     het         45 0.445        2

Corre sin conexión.
[`iei_ejemplo_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_ejemplo_sanger.md)
escribe un par de cromatogramas `.ab1` de un amplicón de BTK con la
variante patogénica c.1632-2A\>G en una portadora.
[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
los alinea contra la referencia GRCh38 que trae el paquete y encuentra
la variante en las dos hebras, con 45 % de la señal en el alelo
alternativo. El cromatograma de arriba es la salida de
[`iei_cromatograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cromatograma.md).

Con conexión,
[`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
y
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
completan el análisis. La variante queda en prioridad alta con su razón
escrita punto por punto.

``` r

iei_priorizar(iei_anotar(v), fenotipo = c("agammaglobulinemia", "neumonía recurrente"))
```

## Entradas

| Entrada | Formato | Función |
|----|----|----|
| VCF | `.vcf`, `.vcf.gz`, GRCh37 o GRCh38 | [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md) |
| Cromatograma Sanger | `.ab1`, `.scf` | [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md) |
| FASTA | `.fa`, `.fasta` | [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md) |
| Secuencia escrita | Texto con las bases | [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md) |
| Informe | PDF, Word o texto | [`iei_variantes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variantes.md) |

## Galería

![Cromatograma de una variante
heterocigota](reference/figures/galeria-1-cromatograma.png)![Dot plot de
una lectura Sanger frente a la referencia de
BTK](reference/figures/galeria-2-dotplot.png)![Genes candidatos por
fenotipo sobre el genoma en
círculo](reference/figures/galeria-3-circos.png)![Diagrama de cuerdas
entre signos del caso y
genes](reference/figures/galeria-4-acordes.png)![Genes candidatos sobre
los cromosomas](reference/figures/galeria-5-cariograma.png)![Variantes
patogénicas de ClinVar sobre la proteína
BTK](reference/figures/galeria-6-clinvar.png)

## Funciones

| Función | Qué hace |
|----|----|
| [`iei_analizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_analizar.md) | Corre el flujo completo de un VCF y deja un reporte por muestra |
| [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md), [`iei_verificar_sexo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_verificar_sexo.md) | Leer un VCF y comparar el sexo declarado con el del cromosoma X |
| [`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md), [`iei_bitacora()`](https://gentlenmoron.github.io/ieiprio/reference/iei_bitacora.md) | Filtrar por panel y calidad, con el registro de cada paso |
| [`iei_cnv()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cnv.md) | CNV que trae el VCF cruzadas con el panel |
| [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md) | Llamar variantes desde `.ab1`, `.scf`, FASTA o una secuencia escrita |
| [`iei_cromatograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cromatograma.md) | Trazas alrededor de una variante, con la base leída y la de referencia |
| [`iei_dotplot()`](https://gentlenmoron.github.io/ieiprio/reference/iei_dotplot.md) | Dot plot de una secuencia frente a otra o frente a la referencia |
| [`iei_consenso()`](https://gentlenmoron.github.io/ieiprio/reference/iei_consenso.md) | Unir lecturas directa y reversa en una secuencia |
| [`iei_leer_ab1()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_ab1.md), [`iei_leer_scf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_scf.md) | Leer cromatogramas |
| [`iei_variantes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variantes.md), [`iei_variante()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variante.md) | Variantes desde PDF, Word, texto o HGVS |
| [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md), [`iei_gnomad()`](https://gentlenmoron.github.io/ieiprio/reference/iei_gnomad.md) | Anotar con Ensembl VEP, ClinVar y gnomAD |
| [`iei_filtrar_frecuencia()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar_frecuencia.md) | Descartar variantes comunes según la herencia del gen |
| [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md), [`iei_pesos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_pesos.md) | Puntaje explicado y criterios ACMG por puntos |
| [`iei_reporte()`](https://gentlenmoron.github.io/ieiprio/reference/iei_reporte.md) | Reporte HTML de un paciente |
| [`iei_exportar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_exportar.md), [`iei_phenopacket()`](https://gentlenmoron.github.io/ieiprio/reference/iei_phenopacket.md), [`iei_granges()`](https://gentlenmoron.github.io/ieiprio/reference/iei_granges.md) | Excel, VCF, BED, Phenopacket o `GRanges` |
| [`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md), [`iei_buscar_hpo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_hpo.md) | Genes del panel por parecido con el fenotipo |
| [`iei_validar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_validar_fenotipo.md) | Medir la búsqueda por fenotipo con casos simulados |
| [`iei_gen()`](https://gentlenmoron.github.io/ieiprio/reference/iei_gen.md) | Ficha de un gen con sus variantes de ClinVar |
| [`iei_circos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_circos.md), [`iei_cariograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cariograma.md) | Panel, candidatos o variantes sobre el genoma |
| [`iei_panel()`](https://gentlenmoron.github.io/ieiprio/reference/iei_panel.md), [`iei_buscar_gen()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_gen.md), [`iei_fuentes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_fuentes.md) | Panel de genes y versiones de los datos |
| [`iei_app()`](https://gentlenmoron.github.io/ieiprio/reference/iei_app.md) | Aplicación Shiny |

Las guías están en <https://gentlenmoron.github.io/ieiprio/>. [Primeros
pasos](https://gentlenmoron.github.io/ieiprio/articles/ieiprio.html)
recorre cada entrada, [Interpretación de
resultados](https://gentlenmoron.github.io/ieiprio/articles/interpretar.html)
explica el puntaje y la búsqueda por fenotipo, [Método y
parámetros](https://gentlenmoron.github.io/ieiprio/articles/metodo.html)
trae los cálculos y [Casos de
ejemplo](https://gentlenmoron.github.io/ieiprio/articles/casos.html)
tres análisis completos.

> Herramienta de apoyo e investigación. No es un informe diagnóstico ni
> reemplaza la clasificación ACMG/AMP hecha por un especialista.

## Cita

Rios Vasquez Y, Pavini N (2026). ieiprio: Priorización de variantes en
genes de errores innatos de la inmunidad. Paquete de R versión 1.0.0.
<https://github.com/Gentlenmoron/ieiprio>

`citation("ieiprio")` devuelve la referencia en BibTeX. El paquete se
distribuye con licencia MIT.
