# Versiones fijas de cada fuente de datos.
# Cambiar una version aqui y volver a correr construir_datos.R es la UNICA
# forma de actualizar los datos del paquete. Asi cualquiera reproduce lo mismo.
#
# md5: dejar en NULL la primera vez. El script imprime la huella de cada
# archivo descargado; copiala aqui para que las siguientes corridas verifiquen
# que se bajo exactamente el mismo archivo.

config <- list(
  panelapp = list(
    panel_id = 398,          # Primary immunodeficiency or monogenic IBD
    version  = "8.78",
    md5      = NULL
  ),
  hpo = list(
    tag = "v2026-09-01",     # release de github.com/obophenotype/human-phenotype-ontology
    md5 = NULL
  )
)
