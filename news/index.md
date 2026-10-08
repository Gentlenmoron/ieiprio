# Changelog

## ieiprio 0.1.0

- Primera versión pública. Flujo completo desde el VCF hasta un reporte
  HTML por paciente, con datos reproducibles de PanelApp, HPO, Ensembl
  VEP y ClinVar.

## ieiprio 0.0.7.9000

- Nueva funcion
  [`iei_analizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_analizar.md)
  que corre todo el flujo en una linea y genera un reporte por muestra.
- Nueva funcion
  [`iei_ejemplo_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_ejemplo_vcf.md)
  con un paciente ficticio que combina cinco variantes reales de ClinVar
  en BTK, RAG1, STAT1 y CYBB con variantes de fondo sinteticas. Se
  genera con `data-raw/ejemplo_vcf.R` y la seleccion queda fijada en
  `data-raw/ejemplo_clinvar.csv`.

## ieiprio 0.0.6.9000

- Nueva funcion
  [`iei_reporte()`](https://gentlenmoron.github.io/ieiprio/reference/iei_reporte.md)
  que genera un HTML autocontenido por paciente, sin Quarto ni otras
  dependencias. Incluye resumen, control de calidad de la muestra
  (histogramas de DP, GQ y fraccion alelica, het/hom, Ti/Tv),
  verificacion de sexo, tabla priorizada con criterios ACMG sugeridos,
  grafico del desglose del puntaje, fichas con diagrama de la proteina y
  sus dominios Pfam, matriz de coincidencia con el fenotipo, embudo del
  filtrado, parametros del analisis con la huella md5 del VCF, fuentes
  con versiones y limitaciones.
- Nueva funcion
  [`iei_verificar_sexo()`](https://gentlenmoron.github.io/ieiprio/reference/iei_verificar_sexo.md)
  que estima el sexo desde el cromosoma X y lo compara con el declarado.
- [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  agrega el desglose del puntaje (`p_clinvar`, `p_impacto`,
  `p_frecuencia`, `p_herencia`, `p_fenotipo`) y la columna `acmg` con
  criterios sugeridos (PVS1, PM2_Supporting, PM3, PP4, BA1). Un portador
  en gen recesivo nunca queda en prioridad alta.
- [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
  agrega `proteina_id` y `pos_proteina`. La cache cambia de version y la
  anterior se ignora.
- [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md)
  guarda el sexo declarado y la huella md5 del VCF, y
  [`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md)
  el control de calidad y los parametros usados. Todo viaja hasta el
  reporte.

## ieiprio 0.0.5.9000

- Nueva funcion
  [`iei_priorizar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_priorizar.md)
  que asigna puntaje, categoria (alta, media, baja), razones en texto y
  notas clinicas. Combina ClinVar, impacto, ganancia de funcion,
  frecuencia, compatibilidad con la herencia (incluye posibles
  heterocigotos compuestos y portadores) y fenotipo HPO o categoria
  IUIS. Descarta benignas y lo registra en la bitacora.
- Nueva funcion
  [`iei_pesos()`](https://gentlenmoron.github.io/ieiprio/reference/iei_pesos.md)
  para ajustar los puntos.
- Los tests de anotacion ya no imprimen mensajes.

## ieiprio 0.0.4.9000

- Nueva funcion
  [`iei_anotar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_anotar.md)
  que consulta la API REST de Ensembl VEP (GRCh37 o GRCh38 segun el VCF)
  en lotes de hasta 200, elige un transcrito por variante (mismo gen,
  MANE Select, canonico, mayor impacto) y agrega consecuencia, HGVS,
  frecuencia maxima en gnomAD y ClinVar. Guarda una cache local y
  registra el release de Ensembl.
- Nueva funcion
  [`iei_filtrar_frecuencia()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar_frecuencia.md)
  con umbrales por modo de herencia, que suma el paso a la bitacora.
- Nueva funcion
  [`iei_limpiar_cache()`](https://gentlenmoron.github.io/ieiprio/reference/iei_limpiar_cache.md).

## ieiprio 0.0.3.9000

- Nueva funcion
  [`iei_filtrar()`](https://gentlenmoron.github.io/ieiprio/reference/iei_filtrar.md)
  que filtra por genotipo, por genes del panel (con margen para
  splicing) y por calidad, dejando una bitacora que se consulta con
  [`iei_bitacora()`](https://gentlenmoron.github.io/ieiprio/reference/iei_bitacora.md).
- Nuevo dataset `genes_coordenadas` con la ubicacion de cada gen en
  GRCh37 y GRCh38, tomada de PanelApp.

## ieiprio 0.0.2.9000

- Nueva funcion
  [`iei_leer_vcf()`](https://gentlenmoron.github.io/ieiprio/reference/iei_leer_vcf.md)
  que lee VCF comprimidos o no, separa sitios multialelicos, normaliza
  cromosomas, detecta el build (GRCh37 o GRCh38) y marca hemicigotos en
  el X de varones fuera de las regiones pseudoautosomicas. Escrita en R
  base, sin dependencias nuevas.

## ieiprio 0.0.1.9000

- Los datos ahora se descargan por API con versiones fijas (PanelApp
  panel 398 y HPO), definidas en `data-raw/config.R`.
- Nuevos datasets `iuis_categorias`, `hpo_genes` y `fuentes_datos`.
- Nuevas funciones
  [`iei_fuentes()`](https://gentlenmoron.github.io/ieiprio/reference/iei_fuentes.md)
  e
  [`iei_actualizar_panel()`](https://gentlenmoron.github.io/ieiprio/reference/iei_actualizar_panel.md).
- [`iei_panel()`](https://gentlenmoron.github.io/ieiprio/reference/iei_panel.md)
  filtra por nivel de evidencia de PanelApp.

## ieiprio 0.0.0.9000

- Esqueleto del paquete con panel semilla.
