# Changelog

## ieiprio 1.0.2

- `variantes_practica` suma las columnas `transcrito` (RefSeq MANE
  Select) y `clinvar_id` (identificador VCV de ClinVar), verificadas
  contra ClinVar el 2026-10-20. Cada fila se puede rastrear hasta su
  registro en ClinVar.
- La documentación del banco aclara que RAG1 en heterocigosis y
  TNFRSF13B c.310T\>C son contraejemplos a propósito, y que el banco
  sirve para practicar el flujo y no para medir rendimiento.
- La prueba en vivo de BTK se salta cuando Ensembl no devuelve la
  proteína, en vez de fallar.
- `bslib (>= 0.9.0)` en Suggests, la versión mínima que usa la
  aplicación.

## ieiprio 1.0.1

- Nuevo dato `variantes_practica`, un banco de doce variantes reales de
  ClinVar en genes conocidos y menos conocidos del panel (BTK, RAG1,
  ADA, STAT3, PIK3CD, TNFRSF13B, FOXP3, AIRE, CYBB, NCF1, STAT1 y MEFV).
  Cada una trae su genotipo típico, la enfermedad, tres signos HPO de
  esa enfermedad y qué aspecto del análisis muestra.
- La guía Primeros pasos suma la sección “Datos para practicar”, con el
  banco, cómo explorar cualquier gen del panel, cómo simular pacientes
  desde `hpo_enfermedades`, cómo armar secuencias propias desde la
  referencia y las fuentes públicas que entran en el paquete.
- Portada del sitio al estilo de catGenes, con los gráficos dentro del
  diagrama de flujo.

## ieiprio 1.0.0

La versión 1.0.0 vuelve al análisis de secuencias. El paquete trabaja
solo con lo que sale de la secuencia y de bases de datos curadas, y deja
fuera la información que no se puede verificar en ella.

### Qué salió

Se eliminó `iei_laboratorio()` con sus rangos por edad
(`referencias_laboratorio`), su gráfico y la pestaña de la aplicación.
También se eliminó `iei_familia()`, el argumento `ped` de
[`iei_analizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_analizar.md)
y todo lo que leía tríos o archivos PED. La priorización ya no usa el
origen materno o paterno, y dejan de sugerirse PS2, PM6, PP1 y BS4. La
fase del propio VCF se conserva, así que un compuesto en trans sigue
sumando PM3 y un par en cis se descarta. Las guías ya no anuncian CNV
por cobertura, ARN ni lecturas largas.

### Qué volvió

[`iei_cromatograma()`](https://gentlenmoron.github.io/ieiprio/reference/iei_cromatograma.md)
recupera el aspecto de su primera versión, con fondo limpio, trazas
finas y la base leída arriba de la de referencia. El nuevo argumento
`estilo = "iei"` dibuja con el tema del resto de los gráficos. El
cromatograma de una portadora de BTK c.1632-2A\>G es la imagen de la
portada.

### Qué se agregó

[`iei_dotplot()`](https://gentlenmoron.github.io/ieiprio/reference/iei_dotplot.md)
compara una secuencia con otra o con la referencia del gen por el método
clásico de ventanas, escrito en R base. Marca la variante y muestra los
cortes por indels como diagonales desplazadas.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) sobre un
resultado de
[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
dibuja el cromatograma o el dot plot.

[`iei_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_sanger.md)
acepta una `referencia` propia, como FASTA con la región en la cabecera
o como lista. El paquete trae la referencia GRCh38 del amplicón de BTK y
su anotación de Ensembl VEP en `inst/extdata`, así que
[`iei_ejemplo_sanger()`](https://gentlenmoron.github.io/ieiprio/reference/iei_ejemplo_sanger.md)
y los casos de las guías corren sin conexión.

[`iei_reporte()`](https://gentlenmoron.github.io/ieiprio/reference/iei_reporte.md)
se adapta a la entrada. El control de calidad del VCF aparece solo con
un VCF, y con Sanger el reporte trae el cromatograma de cada variante
como imagen dentro del HTML.

La galería tiene seis imágenes en 4:3, a 1200 por 900 píxeles, generadas
por `data-raw/figuras_readme.R`.
[`iei_circos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_circos.md)
y el diagrama de cuerdas ponen la leyenda al costado en un lienzo
apaisado. Hay una cuarta guía, Casos de ejemplo, con tres análisis
completos, y las otras tres abren con un ejemplo mínimo.

### Correcciones

El mínimo pasa a R 4.2, porque GlobalOptions, que usa circlize, no
compila en R 4.1. La integración continua fija `ubuntu-24.04` en todos
los flujos y quita la entrada de R 4.1. La aplicación usa
`navbar_options()` de bslib, el gráfico de genes ya no usa el argumento
`label.size` y la conversión a plotly no avisa por las etiquetas. Las
piruletas de ClinVar rotulan con el cambio de aminoácido y la leyenda
del cariograma ya no se corta.

### Versiones anteriores

La 0.1 leía, filtraba, anotaba con Ensembl VEP y priorizaba variantes de
un VCF con un reporte HTML. La 0.2 sumó predictores, variantes desde
texto, PDF y Word, y secuencias Sanger. La 0.3 agregó la exploración de
genes y fenotipos con similitud semántica de HPO y las categorías IUIS.
La 0.4 comparó el fenotipo con cada enfermedad, sugirió criterios ACMG
por puntos, leyó CNV del VCF y mejoró el análisis Sanger y la
aplicación. La 0.5 sumó los gráficos del genoma y la teoría de las
guías.
