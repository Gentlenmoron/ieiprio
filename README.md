# ieiprio

<!-- badges: start -->
[![R-CMD-check](https://github.com/TU_USUARIO/ieiprio/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TU_USUARIO/ieiprio/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Priorización de variantes en genes de errores innatos de la inmunidad (IEI).
ieiprio lee un VCF, se queda con las variantes en genes de la clasificación
IUIS, las anota con Ensembl VEP y las ordena según herencia y fenotipo clínico.

> **Aviso.** Herramienta de apoyo e investigación. No es un informe diagnóstico
> ni reemplaza la clasificación ACMG/AMP hecha por un especialista.

## Estado

En desarrollo. Por ahora incluye el panel de genes (versión semilla) y las
funciones para consultarlo.

## Instalación

```r
# install.packages("pak")
pak::pak("TU_USUARIO/ieiprio")
```

## Uso

```r
library(ieiprio)

iei_panel(categoria = 3)          # deficiencias de anticuerpos
iei_panel(herencia = "XL")        # genes ligados al X
iei_buscar_gen(c("TACI", "STAT1"))
```
