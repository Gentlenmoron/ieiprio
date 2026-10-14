# Consultar gnomAD para las variantes candidatas

Agrega datos de gnomAD que Ensembl VEP no entrega: cuantas personas
sanas son homocigotas o hemicigotas para la variante (la prueba mas
fuerte contra una variante en un gen recesivo o ligado al X, criterio
BS2), la frecuencia filtrada al 95 % que recomienda ClinGen para BA1 y
BS1, y la restriccion del gen (intolerancia a perdida de funcion y a
missense, criterio PP2).

## Usage

``` r
iei_gnomad(
  x,
  max_variantes = 50,
  cache = tools::R_user_dir("ieiprio", "cache")
)
```

## Arguments

- x:

  Tabla anotada con
  [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md).

- max_variantes:

  Maximo de variantes a consultar; las demas quedan en NA.

- cache:

  Carpeta de la cache, o `FALSE`.

## Value

`x` con `gnomad_hom`, `gnomad_hemi`, `faf95`, `grupo_faf95`, `loeuf`
(limite superior de observado/esperado para perdida de funcion; menor a
0.6 indica intolerancia) y `mis_z` (mayor a 3.09 indica intolerancia a
missense).

## Details

gnomAD limita las consultas a unas diez por minuto, asi que conviene
aplicarlo solo a las candidatas: despues de
[`iei_filtrar_frecuencia()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar_frecuencia.md)
o a las primeras filas de un resultado, y luego volver a priorizar.

## Examples

``` r
if (FALSE) { # \dontrun{
anot <- iei_gnomad(iei_filtrar_frecuencia(iei_anotar(filtrado)))
} # }
```
