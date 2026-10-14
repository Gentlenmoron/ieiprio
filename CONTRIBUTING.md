# Cómo contribuir a ieiprio

Gracias por querer mejorar ieiprio. Toda ayuda sirve, desde avisar de un
error hasta proponer un gen que falta en la curación.

## Avisar de un error o proponer una mejora

Abre un [issue](https://github.com/Gentlenmoron/ieiprio/issues) con la
plantilla que corresponda. Incluye lo que corriste, el mensaje completo
y la salida de
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

**Nunca adjuntes datos de pacientes.** Ni VCF, ni cromatogramas, ni
informes, ni códigos que permitan identificar a alguien. Si hace falta
un ejemplo, escribe la variante a mano (`"BTK c.1632-2A>G"`) o usa los
datos de ejemplo del paquete.

## Proponer cambios en el código

1.  Haz un *fork* del repositorio y crea una rama con un nombre
    descriptivo.
2.  Instala las dependencias de desarrollo con
    `pak::pak(c("devtools", "roxygen2", "testthat"))`.
3.  Haz el cambio y agrega pruebas en `tests/testthat/`. Las pruebas no
    deben depender de internet; usa `local_mocked_bindings()` para
    simular los servicios, como hacen las pruebas existentes.
4.  Corre `devtools::document()`, `devtools::test()` y
    `devtools::check()`.
5.  Agrega una línea en `NEWS.md` y abre un *pull request* explicando
    qué cambia y por qué.

## Estilo

- Funciones públicas con el prefijo `iei_` y nombres en español;
  funciones internas con un punto al inicio (`.normalizar_variantes()`).
- El código fuente va en ASCII. Las tildes se escriben como `\u00e1`
  dentro de las cadenas y se omiten en los comentarios.
- Los mensajes al usuario con `cli`, en español claro y con tildes.
- Cada decisión que cambia un puntaje debe quedar escrita en `razones` o
  en `nota`.

## Datos curados

Las tablas de `data-raw/` (`iuis_curado.csv`, `iuis_mecanismo.csv`,
`genes_dificiles.csv`) se pueden corregir o ampliar con un *pull
request*. Cada fila nueva necesita la referencia que la justifica.
Después de editar, corre `data-raw/construir_datos.R` para regenerar los
datos.

## Autoría

Quien haga un aporte sustancial se agrega como autor en `DESCRIPTION`.
