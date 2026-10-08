# Priorizar variantes

Asigna a cada variante un puntaje y una categoria (alta, media o baja),
con una columna `razones` que explica de donde sale cada punto. Las
variantes clasificadas como benignas en ClinVar se descartan y quedan en
la bitacora.

## Usage

``` r
iei_priorizar(
  x,
  fenotipo = NULL,
  categoria = NULL,
  pesos = iei_pesos(),
  panel = ieiprio::iuis_panel,
  hpo = ieiprio::hpo_genes,
  categorias = ieiprio::iuis_categorias
)
```

## Arguments

- x:

  Tabla devuelta por
  [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
  o
  [`iei_filtrar_frecuencia()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar_frecuencia.md).

- fenotipo:

  Terminos HPO del paciente, por ejemplo
  `c("HP:0004313", "HP:0002719")`.

- categoria:

  Categorias IUIS sospechadas (1 a 10), alternativa a `fenotipo`.

- pesos:

  Lista de pesos, ver
  [`iei_pesos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_pesos.md).

- panel, hpo, categorias:

  Tablas de referencia. Por defecto los datos del paquete.

## Value

La tabla con el desglose del puntaje (`p_clinvar`, `p_impacto`,
`p_frecuencia`, `p_herencia`, `p_fenotipo`), el `puntaje` total, la
`categoria`, las `razones`, una `nota` clinica y `acmg`, con los
criterios ACMG/AMP sugeridos (PVS1, PM2_Supporting, PM3, PP4, BA1). Los
criterios son orientativos y deben ser revisados por un especialista.
Ordenada por muestra y puntaje.

## Details

Los componentes del puntaje son estos.

- **ClinVar.** Patogenica o probablemente patogenica suma 4. Con
  interpretaciones en conflicto suma 1.

- **Impacto.** Alto (stop, frameshift, splicing canonico) suma 3.
  Moderado (missense, inframe) suma 1. En genes con enfermedad por
  ganancia de funcion una missense cuenta como impacto alto.

- **Frecuencia.** Ausente de gnomAD suma 1.

- **Herencia.** Genotipo compatible con el modo del gen suma 2. En genes
  autosomicos recesivos, dos heterocigotos en el mismo gen y la misma
  muestra se marcan como posible heterocigoto compuesto y suman 2; un
  heterocigoto solo resta 2, se marca como portador y nunca queda en
  prioridad alta, aunque el fenotipo encaje.

- **Fenotipo.** Proporcion de los terminos HPO del paciente presentes en
  el gen, escalada de 0 a 3. O 3 puntos si el gen esta en la `categoria`
  IUIS indicada.

## Examples

``` r
if (FALSE) { # \dontrun{
prior <- iei_priorizar(anot, fenotipo = c("HP:0004313", "HP:0002719"))
} # }
```
