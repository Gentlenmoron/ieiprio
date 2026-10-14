# ieiprio 1.0.0

La versión 1.0.0 vuelve al análisis de secuencias. El paquete trabaja solo con
lo que sale de la secuencia y de bases de datos curadas, y deja fuera la
información que no se puede verificar en ella.

## Qué salió

Se eliminó `iei_laboratorio()` con sus rangos por edad (`referencias_laboratorio`),
su gráfico y la pestaña de la aplicación. También se eliminó `iei_familia()`,
el argumento `ped` de `iei_analizar()` y todo lo que leía tríos o archivos PED.
La priorización ya no usa el origen materno o paterno, y dejan de sugerirse
PS2, PM6, PP1 y BS4. La fase del propio VCF se conserva, así que un compuesto en
trans sigue sumando PM3 y un par en cis se descarta. Las guías ya no anuncian
CNV por cobertura, ARN ni lecturas largas.

## Qué volvió

`iei_cromatograma()` recupera el aspecto de su primera versión, con fondo
limpio, trazas finas y la base leída arriba de la de referencia. El nuevo
argumento `estilo = "iei"` dibuja con el tema del resto de los gráficos. El
cromatograma de una portadora de BTK c.1632-2A>G es la imagen de la portada.

## Qué se agregó

`iei_dotplot()` compara una secuencia con otra o con la referencia del gen por
el método clásico de ventanas, escrito en R base. Marca la variante y muestra
los cortes por indels como diagonales desplazadas. `plot()` sobre un resultado
de `iei_sanger()` dibuja el cromatograma o el dot plot.

`iei_sanger()` acepta una `referencia` propia, como FASTA con la región en la
cabecera o como lista. El paquete trae la referencia GRCh38 del amplicón de BTK
y su anotación de Ensembl VEP en `inst/extdata`, así que
`iei_ejemplo_sanger()` y los casos de las guías corren sin conexión.

`iei_reporte()` se adapta a la entrada. El control de calidad del VCF aparece
solo con un VCF, y con Sanger el reporte trae el cromatograma de cada variante
como imagen dentro del HTML.

La galería tiene seis imágenes en 4:3, a 1200 por 900 píxeles, generadas por
`data-raw/figuras_readme.R`. `iei_circos()` y el diagrama de cuerdas ponen la
leyenda al costado en un lienzo apaisado. Hay una cuarta guía, Casos de
ejemplo, con tres análisis completos, y las otras tres abren con un ejemplo
mínimo.

## Correcciones

El mínimo pasa a R 4.2, porque GlobalOptions, que usa circlize, no compila en
R 4.1. La integración continua fija `ubuntu-24.04` en todos los flujos y
quita la entrada de R 4.1. La aplicación usa `navbar_options()` de bslib, el
gráfico de genes ya no usa el argumento `label.size` y la conversión a plotly
no avisa por las etiquetas. Las piruletas de ClinVar rotulan con el cambio de
aminoácido y la leyenda del cariograma ya no se corta.

## Versiones anteriores

La 0.1 leía, filtraba, anotaba con Ensembl VEP y priorizaba variantes de un VCF
con un reporte HTML. La 0.2 sumó predictores, variantes desde texto, PDF y Word,
y secuencias Sanger. La 0.3 agregó la exploración de genes y fenotipos con
similitud semántica de HPO y las categorías IUIS. La 0.4 comparó el fenotipo con
cada enfermedad, sugirió criterios ACMG por puntos, leyó CNV del VCF y mejoró
el análisis Sanger y la aplicación. La 0.5 sumó los gráficos del genoma y la
teoría de las guías.
