# ieiprio <img src="man/figures/logo.svg" align="right" height="110" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Gentlenmoron/ieiprio/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Gentlenmoron/ieiprio/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Gentlenmoron/ieiprio/actions/workflows/pkgdown.yaml/badge.svg)](https://gentlenmoron.github.io/ieiprio/)
<!-- badges: end -->

**Priorización de variantes en genes de errores innatos de la inmunidad.**

ieiprio toma un VCF de exoma o de panel y devuelve las variantes en genes de
errores innatos de la inmunidad (EII) ordenadas por prioridad, con una
explicación escrita de cada puntaje y un reporte HTML por paciente.

> **Aviso.** Herramienta de apoyo e investigación. No es un informe diagnóstico
> ni reemplaza la clasificación ACMG/AMP hecha por un especialista.

## Instalación

```r
# install.packages("pak")
pak::pak("Gentlenmoron/ieiprio")
```

## Uso

```r
library(ieiprio)

res <- iei_analizar(
  iei_ejemplo_vcf(),                      # o la ruta a tu VCF
  sexo = "M",
  fenotipo = c("HP:0004313", "HP:0002719"),
  carpeta = "reportes"
)
```

Eso lee el VCF, filtra por genes del panel y calidad, anota con Ensembl VEP,
filtra por frecuencia en gnomAD, prioriza y deja un reporte HTML por muestra.
La [guía de uso](https://gentlenmoron.github.io/ieiprio/articles/ieiprio.html)
recorre cada paso con un paciente de ejemplo.

## Qué hace

| Paso | Función | Qué aporta |
| --- | --- | --- |
| Leer | `iei_leer_vcf()` | VCF comprimido o no, multialélicos, GRCh37 o GRCh38, hemicigotos en el X |
| Verificar | `iei_verificar_sexo()` | Estima el sexo desde el X y lo compara con el declarado |
| Filtrar | `iei_filtrar()` | Genes del panel, profundidad, calidad y fracción alélica, con bitácora |
| Anotar | `iei_anotar()` | Consecuencia, HGVS, gnomAD y ClinVar vía Ensembl VEP, con caché |
| Priorizar | `iei_priorizar()` | Puntaje explicado, herencia, compuestos, fenotipo HPO y ACMG sugerido |
| Reportar | `iei_reporte()` | HTML con control de calidad, gráficos, diagramas de proteína y trazabilidad |

## Datos y reproducibilidad

El panel sale del [panel 398 de PanelApp](https://panelapp.genomicsengland.co.uk/panels/398/)
y los fenotipos de la [Human Phenotype Ontology](https://hpo.jax.org), ambos
fijados a una versión y con huella md5. `iei_fuentes()` muestra qué versión usa
el paquete instalado.

## Limitaciones

Solo analiza variantes puntuales e indels pequeños. No detecta CNV ni variantes
estructurales y no evalúa cobertura, porque eso requiere el BAM.
