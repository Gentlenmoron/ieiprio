# Fenotipos por enfermedad de los genes del panel

Signos HPO de cada enfermedad (OMIM u Orphanet) asociada a un gen del
panel, con su frecuencia y los signos que la enfermedad excluye. Es la
base de la similitud por enfermedad de
[`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md)
e
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md),
y de los casos simulados de
[`iei_validar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_validar_fenotipo.md).

## Usage

``` r
hpo_enfermedades
```

## Format

Un tibble con columnas `gen`, `enfermedad_id`, `enfermedad`, `herencia`
(de la anotacion de la enfermedad; `NA` si HPO no la indica), `hpo_id`,
`frecuencia` (proporcion de pacientes con el signo, `NA` si no se
conoce) y `excluido` (`TRUE` cuando la enfermedad no presenta ese signo,
anotacion NOT).

## Source

`phenotype.hpoa` y `genes_to_disease.txt` de la Human Phenotype
Ontology, version indicada en
[fuentes_datos](https://gentlenmoron.github.io/ieiprio/reference/fuentes_datos.md).
<https://github.com/obophenotype/human-phenotype-ontology/releases>
