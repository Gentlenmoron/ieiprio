# Paleta y tema comunes a todos los graficos del paquete.
#
# Los colores categoricos siguen un orden fijo validado para daltonismo en
# pares adyacentes (azul, naranja, aguamarina, amarillo, magenta, verde). Las
# magnitudes usan una sola rampa azul y la patogenicidad una rampa divergente
# rojo a azul con gris en el centro. El azul marino es el color de la marca.

.marca <- "#0b2545"
.marca_2 <- "#1d6fa5"
.tinta <- c(primaria = "#1b1b1a", secundaria = "#52514e", tenue = "#8a8984",
            rejilla = "#ecebe7", fondo = "#fcfcfb", eje = "#d3d2cd")
.categoricos <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
.rampa_azul <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
.divergente <- c(Patogenica = "#a32626", `Probablemente patogenica` = "#e66767",
                 `Significado incierto` = "#c9c8c2", `Probablemente benigna` = "#6da7ec",
                 Benigna = "#184f95")
.formas_tipo <- c("missense" = 21, "sin sentido" = 24, "frameshift o indel" = 22, "splicing" = 23,
                  "sinonima" = 25, "otra" = 21)

.tema_iei <- function(base_size = 12, leyenda = "top") {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.3, colour = .marca,
                                         margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(colour = .tinta[["secundaria"]], size = base_size * 0.95,
                                            margin = ggplot2::margin(b = 10)),
      plot.caption = ggplot2::element_text(colour = .tinta[["tenue"]], size = base_size * 0.75, hjust = 0,
                                           margin = ggplot2::margin(t = 10)),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      axis.title = ggplot2::element_text(colour = .tinta[["secundaria"]], size = base_size * 0.85),
      axis.text = ggplot2::element_text(colour = .tinta[["secundaria"]]),
      panel.grid.major = ggplot2::element_line(colour = .tinta[["rejilla"]], linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = leyenda,
      legend.justification = "left",
      legend.title = ggplot2::element_text(colour = .tinta[["secundaria"]], size = base_size * 0.85),
      legend.text = ggplot2::element_text(colour = .tinta[["primaria"]], size = base_size * 0.85),
      legend.key.size = ggplot2::unit(0.9, "lines"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(14, 18, 10, 14)
    )
}

.colores_tipo <- stats::setNames(c(.categoricos[1:4], "#a3a29b", "#c9c8c2"),
                                 c("missense", "sin sentido", "frameshift o indel", "splicing", "sinonima", "otra"))

# Colores por herencia en orden fijo, con gris para lo desconocido.
.colores_herencia <- c(AR = .categoricos[1], AD = .categoricos[2], XL = .categoricos[3],
                       `AR/AD` = .categoricos[4], `AD/AR` = .categoricos[4], MT = .categoricos[5],
                       `sin dato` = "#a3a29b")

# Texto claro sobre fondos oscuros y oscuro sobre fondos claros.
.texto_sobre <- function(valor, corte = 0.55) ifelse(valor >= corte, "white", .tinta[["primaria"]])

.acortar <- function(x, n = 42) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)

.pct <- function(x) paste0(round(100 * x), " %")
