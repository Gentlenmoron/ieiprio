# Categorias IUIS de los genes del panel

Categoria de la clasificacion IUIS de errores innatos de la inmunidad de
cada gen. Un gen puede estar en mas de una. La mayoria sale de los
fenotipos de PanelApp, que nombran la tabla IUIS del gen; cuando no la
nombran, de palabras clave de esos fenotipos (por ejemplo SCID o
agammaglobulinemia); y unas pocas de una curacion manual. La columna
`fuente` dice de donde salio cada una.

## Usage

``` r
iuis_categorias
```

## Format

Un tibble con columnas `gen`, `categoria_iuis` (1 a 10),
`categoria_nombre` y `fuente` (`"PanelApp, tabla IUIS"`,
`"PanelApp, fenotipo"` o `"curado"`).

## Source

PanelApp panel 398 y clasificacion IUIS (Poli y colaboradores).
