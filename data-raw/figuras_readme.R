# Figuras del README y de la galeria del sitio. Todas salen en 4:3, 1200 x 900
# px (8 x 6 pulgadas a 150 dpi) con ragg. Corre de punta a punta; la ficha de
# BTK necesita conexion y, si no la hay, la ultima figura pasa a ser el
# ranking por fenotipo.
devtools::load_all()
library(ggplot2)
options(ieiprio.silenciar_privacidad = TRUE)
dir <- "man/figures"
ancho <- 8
alto <- 6
dpi <- 150
guardar <- function(g, f) {
  ggsave(file.path(dir, f), g, width = ancho, height = alto, dpi = dpi, device = ragg::agg_png, bg = "white")
}
base_png <- function(f) ragg::agg_png(file.path(dir, f), width = ancho, height = alto, units = "in", res = dpi)

unlink(Sys.glob(file.path(dir, "galeria-*.png")))

ref <- system.file("extdata", "btk_referencia.fa", package = "ieiprio")
het <- suppressMessages(iei_sanger(iei_ejemplo_sanger("het"), gen = "BTK", sexo = "F", referencia = ref))
cromatograma <- iei_cromatograma(het, ventana = 9)
guardar(cromatograma, "cromatograma.png")
guardar(cromatograma, "galeria-1-cromatograma.png")

hemi <- suppressMessages(iei_sanger(iei_ejemplo_sanger("hemi"), gen = "BTK", sexo = "M", referencia = ref))
guardar(iei_dotplot(hemi), "galeria-2-dotplot.png")

options(ieiprio.sin_red = TRUE)
b <- iei_buscar_fenotipo(c("agammaglobulinemia", "neumonía recurrente", "diarrea crónica"), n = 20)
base_png("galeria-3-circos.png"); iei_circos(b); grDevices::dev.off()
base_png("galeria-4-acordes.png"); plot(b, tipo = "acordes"); grDevices::dev.off()
guardar(iei_cariograma(b), "galeria-5-cariograma.png")

options(ieiprio.sin_red = FALSE)
ficha <- tryCatch(iei_gen("BTK"), error = function(e) NULL)
if (!is.null(ficha) && !is.null(ficha$clinvar)) {
  guardar(plot(ficha), "galeria-6-clinvar.png")
} else {
  message("Sin conexion para la ficha de BTK; se usa el ranking por fenotipo.")
  guardar(plot(b), "galeria-6-ranking.png")
}
