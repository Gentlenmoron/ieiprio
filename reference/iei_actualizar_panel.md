# Descargar una version actualizada del panel desde PanelApp

El panel incluido en el paquete esta fijado a una version de PanelApp
para que los resultados sean reproducibles (ver
[`iei_fuentes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_fuentes.md)).
Esta funcion baja otra version, o la mas reciente, sin modificar los
datos del paquete, y la guarda en la cache del usuario.

## Usage

``` r
iei_actualizar_panel(version = NULL, panel_id = 398, guardar = TRUE)
```

## Arguments

- version:

  Version del panel en PanelApp, por ejemplo `"8.78"`. `NULL` baja la
  mas reciente.

- panel_id:

  Identificador del panel en PanelApp. Por defecto 398,
  inmunodeficiencias primarias.

- guardar:

  Si `TRUE`, guarda el resultado como RDS en la cache.

## Value

Un tibble con la misma estructura que
[iuis_panel](https://gentlenmoron.github.io/ieiprio/reference/iuis_panel.md),
con los atributos `version` y `panel_id`.

## Examples

``` r
if (FALSE) { # \dontrun{
nuevo <- iei_actualizar_panel()
attr(nuevo, "version")
} # }
```
