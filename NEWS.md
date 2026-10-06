# ieiprio 0.0.5.9000

* Nueva funcion `iei_priorizar()` que asigna puntaje, categoria (alta, media, baja), razones en texto y notas clinicas. Combina ClinVar, impacto, ganancia de funcion, frecuencia, compatibilidad con la herencia (incluye posibles heterocigotos compuestos y portadores) y fenotipo HPO o categoria IUIS. Descarta benignas y lo registra en la bitacora.
* Nueva funcion `iei_pesos()` para ajustar los puntos.
* Los tests de anotacion ya no imprimen mensajes.

# ieiprio 0.0.4.9000

* Nueva funcion `iei_anotar()` que consulta la API REST de Ensembl VEP (GRCh37 o GRCh38 segun el VCF) en lotes de hasta 200, elige un transcrito por variante (mismo gen, MANE Select, canonico, mayor impacto) y agrega consecuencia, HGVS, frecuencia maxima en gnomAD y ClinVar. Guarda una cache local y registra el release de Ensembl.
* Nueva funcion `iei_filtrar_frecuencia()` con umbrales por modo de herencia, que suma el paso a la bitacora.
* Nueva funcion `iei_limpiar_cache()`.

# ieiprio 0.0.3.9000

* Nueva funcion `iei_filtrar()` que filtra por genotipo, por genes del panel (con margen para splicing) y por calidad, dejando una bitacora que se consulta con `iei_bitacora()`.
* Nuevo dataset `genes_coordenadas` con la ubicacion de cada gen en GRCh37 y GRCh38, tomada de PanelApp.

# ieiprio 0.0.2.9000

* Nueva funcion `iei_leer_vcf()` que lee VCF comprimidos o no, separa sitios multialelicos, normaliza cromosomas, detecta el build (GRCh37 o GRCh38) y marca hemicigotos en el X de varones fuera de las regiones pseudoautosomicas. Escrita en R base, sin dependencias nuevas.

# ieiprio 0.0.1.9000

* Los datos ahora se descargan por API con versiones fijas (PanelApp panel 398 y HPO), definidas en `data-raw/config.R`.
* Nuevos datasets `iuis_categorias`, `hpo_genes` y `fuentes_datos`.
* Nuevas funciones `iei_fuentes()` e `iei_actualizar_panel()`.
* `iei_panel()` filtra por nivel de evidencia de PanelApp.

# ieiprio 0.0.0.9000

* Esqueleto del paquete con panel semilla.
