# Leer un cromatograma Sanger (.ab1)

Lee un archivo ABIF de un secuenciador Applied Biosystems y devuelve las
bases llamadas, su calidad, la posicion de cada pico y las cuatro
trazas. Esta escrito en R base, sin dependencias de Bioconductor.

## Usage

``` r
iei_leer_ab1(archivo)

# S3 method for class 'iei_ab1'
print(x, ...)
```

## Arguments

- archivo:

  Ruta al archivo `.ab1` o `.abi`.

- x:

  Un objeto `iei_ab1`.

- ...:

  No se usa.

## Value

Una lista de clase `iei_ab1` con `bases` (texto), `calidad` (Phred por
base, `NA` si el archivo no la trae), `picos` (posicion de cada base en
la traza), `trazas` (tibble con las columnas `A`, `C`, `G`, `T`),
`muestra` y `archivo`.

## Examples

``` r
if (FALSE) { # \dontrun{
ab1 <- iei_leer_ab1("BTK_exon15_F.ab1")
ab1$bases
} # }
```
