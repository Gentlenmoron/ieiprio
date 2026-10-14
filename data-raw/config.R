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
    md5 = "6dac74395b2a9192dc0e71594093594f",      # genes_to_phenotype.txt
    md5_obo = "97620d1eae29ff3665ee644a7315b0f1",  # hp.obo
    md5_hpoa = "b4125a21c9ba916d10ce7551e55ead53", # phenotype.hpoa
    md5_g2d = "337ee42c2c17f2d89ec5d21ba7515c67"   # genes_to_disease.txt
  )
)
