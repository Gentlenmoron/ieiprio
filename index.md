# ieiprio

`ieiprio` encuentra, anota y prioriza variantes en genes de errores
innatos de la inmunidad a partir de secuencias Sanger, FASTA, VCF o el
texto de un informe. Trabaja solo con lo que sale de la secuencia y de
bases de datos curadas, para que el mismo archivo dé siempre el mismo
resultado.

## Características principales

- Lee variantes desde VCF, cromatogramas Sanger (`.ab1`, `.scf`), FASTA,
  una secuencia escrita o el texto de un informe en PDF o Word.
- Llama variantes Sanger con detección de heterocigotos por picos dobles
  y las muestra en el cromatograma, base por base.
- Anota con Ensembl VEP, ClinVar y gnomAD, con caché local y modo sin
  conexión.
- Prioriza con un puntaje explicado en texto y criterios ACMG por
  puntos.
- Compara el fenotipo del caso (términos HPO) con cada enfermedad del
  panel IUIS de PanelApp.
- Entrega un reporte HTML por paciente y exporta a Excel, VCF, BED o
  Phenopacket.
- Trae un banco de variantes reales en genes conocidos y menos
  conocidos, para practicar sin tener una secuencia propia.

## Instalación

``` r

install.packages("pak")
pak::pak("Gentlenmoron/ieiprio")
```

Para tener también las guías dentro de R, instala con `remotes` y
construye las viñetas.

``` r

install.packages("remotes")
remotes::install_github("Gentlenmoron/ieiprio", build_vignettes = TRUE)
browseVignettes("ieiprio")
```

## Entradas

| Entrada | Formato | Función |
|----|----|----|
| VCF | `.vcf`, `.vcf.gz`, GRCh37 o GRCh38 | [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md) |
| Cromatograma Sanger | `.ab1`, `.scf` | [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md) |
| FASTA | `.fa`, `.fasta` | [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md) |
| Secuencia escrita | Texto con las bases | [`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md) |
| Informe | PDF, Word o texto | [`iei_variantes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_variantes.md) |

## Flujo de trabajo

![Flujo de ieiprio. Las entradas VCF, Sanger o FASTA y texto pasan por
filtrar, anotar y priorizar hasta el reporte HTML y la exportación.
Debajo, seis gráficos: cromatograma, dot plot, variantes de ClinVar en
la proteína, genoma en círculo, diagrama de cuerdas y
cariograma.](reference/figures/flujo.svg)

*De la secuencia al reporte. En azul las entradas, en oscuro el proceso
y en verde las salidas. Debajo, el gráfico que acompaña cada paso. Los
ejemplos completos que los generan están en la guía [Casos de
ejemplo](https://gentlenmoron.github.io/ieiprio/articles/casos.html).*

## Funciones disponibles

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
| `variantes_practica` | Banco de doce variantes reales para practicar, con su enfermedad y signos HPO |
| [`iei_app()`](https://gentlenmoron.github.io/ieiprio/reference/iei_app.md) | Aplicación Shiny |

## Documentación

El sitio del paquete, <https://gentlenmoron.github.io/ieiprio/>, reúne
la referencia de cada función y cuatro guías.

| Guía | Contenido |
|----|----|
| [Primeros pasos](https://gentlenmoron.github.io/ieiprio/articles/ieiprio.html) | Instalación, fundamentos, cada tipo de entrada y datos para practicar |
| [Interpretación de resultados](https://gentlenmoron.github.io/ieiprio/articles/interpretar.html) | El puntaje, la clasificación ACMG y la búsqueda por fenotipo |
| [Método y parámetros](https://gentlenmoron.github.io/ieiprio/articles/metodo.html) | Fórmulas, umbrales, validación y referencias |
| [Casos de ejemplo](https://gentlenmoron.github.io/ieiprio/articles/casos.html) | Tres análisis completos con datos para copiar |

Las mismas guías se abren sin internet con `browseVignettes("ieiprio")`
si el paquete se instaló con `build_vignettes = TRUE`. Las insignias de
arriba muestran si las pruebas y la revisión del paquete pasan en macOS,
Windows y Linux.

> Herramienta de apoyo e investigación. No es un informe diagnóstico ni
> reemplaza la clasificación ACMG/AMP hecha por un especialista.

## Cita

Rios Vasquez Y, Pavini N (2026). ieiprio: Priorización de variantes en
genes de errores innatos de la inmunidad. Paquete de R versión 1.0.1.
<https://github.com/Gentlenmoron/ieiprio>

`citation("ieiprio")` devuelve la referencia en BibTeX. El paquete se
distribuye con licencia MIT.
