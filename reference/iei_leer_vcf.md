# Leer un archivo VCF

Lee un VCF (comprimido o no) y lo convierte en una tabla larga con una
fila por variante, alelo alternativo y muestra. Separa los sitios
multialelicos, normaliza los nombres de cromosoma, detecta el build del
genoma y marca como hemicigotas las variantes del X en varones.

## Usage

``` r
iei_leer_vcf(
  ruta,
  muestras = NULL,
  build = c("auto", "GRCh37", "GRCh38"),
  sexo = NULL
)
```

## Arguments

- ruta:

  Ruta al archivo `.vcf` o `.vcf.gz`.

- muestras:

  Nombres de las muestras a leer. `NULL` lee todas.

- build:

  `"auto"` detecta el build desde el encabezado. Tambien se puede
  indicar `"GRCh37"` o `"GRCh38"`.

- sexo:

  Sexo de las muestras, `"M"` o `"F"`. Puede ser un solo valor para
  todas o un vector con nombres de muestra, por ejemplo
  `c(P01 = "M", P02 = "F")`. `NULL` si no se conoce; en ese caso no se
  marcan hemicigotos salvo que el genotipo ya sea haploide.

## Value

Un tibble con las columnas `muestra`, `chr`, `pos`, `id`, `ref`, `alt`,
`genotipo` (`het`, `hom_alt`, `hemi`, `hom_ref` o `NA` si no se llamo),
`qual`, `dp`, `gq`, `ab` (fraccion de lecturas con el alelo
alternativo), `filtro` y `build`. Lleva los atributos `build` y
`archivo`.

## Examples

``` r
if (FALSE) { # \dontrun{
vcf <- iei_leer_vcf("paciente01.vcf.gz", sexo = "M")
} # }
```
