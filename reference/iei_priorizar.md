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
  similitud = c("enfermedad", "semantica", "exacta"),
  fenotipo_ausente = NULL,
  panel = ieiprio::iuis_panel,
  hpo = ieiprio::hpo_genes,
  categorias = ieiprio::iuis_categorias,
  enfermedades = ieiprio::hpo_enfermedades
)
```

## Arguments

- x:

  Tabla devuelta por
  [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
  o
  [`iei_filtrar_frecuencia()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar_frecuencia.md).

- fenotipo:

  Fenotipo del caso, como codigos HPO (`c("HP:0004313", "HP:0002719")`)
  o como texto en espanol o ingles
  (`c("agammaglobulinemia", "neumonia recurrente")`). El texto se
  traduce a HPO y se informa la traduccion; ver
  [`iei_buscar_hpo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_hpo.md).

- categoria:

  Categorias IUIS sospechadas (1 a 10), alternativa a `fenotipo`.

- pesos:

  Lista de pesos, ver
  [`iei_pesos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_pesos.md).

- similitud:

  `"enfermedad"`, `"semantica"` o `"exacta"`, como se compara el
  fenotipo con el del gen.

- fenotipo_ausente:

  Signos que el caso no tiene, en el mismo formato que `fenotipo`.
  Restan a los genes cuyas enfermedades suelen tenerlos.

- panel, hpo, categorias, enfermedades:

  Tablas de referencia. Por defecto los datos del paquete.

## Value

La tabla con el desglose del puntaje (`p_clinvar`, `p_impacto`,
`p_frecuencia`, `p_predictores`, `p_herencia`, `p_fenotipo`), el
`puntaje` total, la `categoria`, las `razones`, una `nota` clinica,
`acmg` con los criterios ACMG/AMP sugeridos (PVS1 graduado, PS1, PM1,
PM2, PM3, PM4, PM5, PP2, PP3, PP4, BA1, BS1, BS2, BP4 y BP7, segun los
datos disponibles), `acmg_puntos` con su suma en el sistema de puntos de
Tavtigian y colaboradores (2020) y `acmg_clase` con la clase que
corresponde. Los criterios son orientativos y deben ser revisados por un
especialista. Ordenada por muestra y puntaje. Es un objeto
`iei_resultado`, que se comporta como tibble y tiene metodos
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) y
[`plot()`](https://rdrr.io/r/graphics/plot.default.html); ver tambien
[`iei_tabla()`](https://gentlenmoron.github.io/ieiprio/reference/iei_tabla.md)
y
[`iei_metadatos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_metadatos.md).

## Details

Los componentes del puntaje son estos.

- **ClinVar.** Patogenica o probablemente patogenica suma 4 con dos o
  mas estrellas de revision (o si no se pudieron verificar), 3 con una
  estrella y 1 sin estrellas. Con interpretaciones en conflicto suma 1.

- **Impacto.** Alto (stop, frameshift, splicing canonico) suma 3.
  Moderado (missense, inframe) suma 1. En genes con enfermedad por
  ganancia de funcion una missense cuenta como impacto alto.

- **Frecuencia.** Ausente de gnomAD suma 1.

- **Predictores.** REVEL, AlphaMissense o SpliceAI por encima de sus
  umbrales calibrados suman 1 (PP3). Si todos los disponibles indican
  efecto benigno restan 1 (BP4).

- **Herencia.** Genotipo compatible con el modo del gen suma 2. En genes
  autosomicos recesivos, dos heterocigotos en el mismo gen y la misma
  muestra se marcan como posible heterocigoto compuesto y suman 2; un
  heterocigoto solo resta 2, se marca como portador y nunca queda en
  prioridad alta, aunque el fenotipo encaje.

- **Fenotipo.** Parecido entre el fenotipo del caso y el del gen,
  escalado de 0 a 3. Cada termino del caso se compara con el mas
  parecido usando la jerarquia de HPO: un termino igual o mas especifico
  vale 1 y uno emparentado vale menos segun cuanta informacion
  comparten, asi "neumonia recurrente" suma para un gen anotado con
  "infecciones respiratorias recurrentes". Con
  `similitud = "enfermedad"` (por defecto) la comparacion se hace con
  cada enfermedad del gen por separado, pesando cada signo por su
  frecuencia y restando los signos que la enfermedad excluye y los
  `fenotipo_ausente`; las razones nombran la enfermedad mas parecida.
  Con `"semantica"` se usan todos los terminos del gen juntos y con
  `"exacta"` solo los identicos (ver
  [`iei_buscar_fenotipo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_buscar_fenotipo.md)).
  Un gen en la `categoria` IUIS indicada recibe 3 puntos.

## Examples

``` r
if (FALSE) { # \dontrun{
prior <- iei_priorizar(anot, fenotipo = c("HP:0004313", "HP:0002719"))
} # }
```
