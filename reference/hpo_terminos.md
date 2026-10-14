# Ontologia HPO

Todos los terminos vigentes de la Human Phenotype Ontology con su
jerarquia, para la similitud semantica de
[`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
y
[`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md)
y la busqueda por texto de
[`iei_buscar_hpo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_hpo.md).

## Usage

``` r
hpo_terminos
```

## Format

Un tibble con una fila por termino:

- id:

  Codigo HPO.

- nombre:

  Nombre en ingles.

- padres:

  Lista con los codigos de sus padres (relacion is_a).

- sinonimos:

  Lista con sus sinonimos exactos.

- alternos:

  Lista con codigos antiguos que hoy apuntan a este termino.

- ic:

  Contenido de informacion: menos logaritmo de la fraccion de genes
  anotados con el termino o uno mas especifico. Mas alto es mas
  especifico.

## Source

Archivos hp.obo y genes_to_phenotype.txt del release indicado en
[`iei_fuentes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_fuentes.md),
<https://github.com/obophenotype/human-phenotype-ontology>.
