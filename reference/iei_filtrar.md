# Filtrar variantes por genotipo, panel y calidad

Reduce la tabla de
[`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md)
a las variantes presentes en genes del panel que pasan los filtros de
calidad. Cada paso queda registrado en una bitacora que se consulta con
[`iei_bitacora()`](https://gentlenmoron.github.io/ieiprio/reference/iei_bitacora.md).

## Usage

``` r
iei_filtrar(
  vcf,
  genes = NULL,
  evidencia = "verde",
  min_dp = 10,
  min_gq = 20,
  min_ab = 0.2,
  solo_pass = TRUE,
  margen = 20,
  coordenadas = ieiprio::genes_coordenadas
)
```

## Arguments

- vcf:

  Tabla devuelta por
  [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md).

- genes:

  Genes a conservar. Por defecto, los de
  [`iei_panel()`](https://gentlenmoron.github.io/ieiprio/reference/iei_panel.md)
  con el nivel de `evidencia` indicado.

- evidencia:

  Niveles de evidencia de PanelApp para el panel por defecto. Se ignora
  si das `genes`.

- min_dp:

  Profundidad minima.

- min_gq:

  Calidad de genotipo minima.

- min_ab:

  Fraccion alelica minima para heterocigotos. El maximo es `1 - min_ab`,
  y homocigotos y hemicigotos deben tener al menos `1 - min_ab`.

- solo_pass:

  Si `TRUE`, solo variantes con FILTER igual a `PASS` o `.`.

- margen:

  Pares de bases que se agregan a cada lado del gen.

- coordenadas:

  Tabla de coordenadas de genes. Por defecto
  [genes_coordenadas](https://gentlenmoron.github.io/ieiprio/reference/genes_coordenadas.md).

## Value

La tabla filtrada con la columna `gen` agregada despues de `muestra`, y
el atributo `ieiprio_bitacora`.

## Details

Los pasos, en orden, son estos.

1.  **Genotipo.** Quita `hom_ref` y genotipos sin llamar.

2.  **Panel.** Deja las variantes dentro de un gen del panel, con un
    margen a cada lado para no perder variantes de splicing. Agrega la
    columna `gen`. Una variante en dos genes solapados aparece una vez
    por gen.

3.  **Calidad.** Exige profundidad, calidad de genotipo, FILTER igual a
    PASS y una fraccion alelica coherente con el genotipo. Un valor
    faltante no descarta la variante.

El filtro por frecuencia poblacional se aplica despues de anotar.

## Examples

``` r
if (FALSE) { # \dontrun{
vcf <- iei_leer_vcf("paciente01.vcf.gz", sexo = "M")
filt <- iei_filtrar(vcf)
iei_bitacora(filt)
} # }
```
