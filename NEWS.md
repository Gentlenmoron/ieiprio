# ieiprio 0.0.2.9000

* Nueva funcion `iei_leer_vcf()` que lee VCF comprimidos o no, separa sitios multialelicos, normaliza cromosomas, detecta el build (GRCh37 o GRCh38) y marca hemicigotos en el X de varones fuera de las regiones pseudoautosomicas. Escrita en R base, sin dependencias nuevas.

# ieiprio 0.0.1.9000

* Los datos ahora se descargan por API con versiones fijas (PanelApp panel 398 y HPO), definidas en `data-raw/config.R`.
* Nuevos datasets `iuis_categorias`, `hpo_genes` y `fuentes_datos`.
* Nuevas funciones `iei_fuentes()` e `iei_actualizar_panel()`.
* `iei_panel()` filtra por nivel de evidencia de PanelApp.

# ieiprio 0.0.0.9000

* Esqueleto del paquete con panel semilla.
