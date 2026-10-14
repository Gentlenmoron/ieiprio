# Aplicacion Shiny de ieiprio. Se abre con ieiprio::iei_app().
library(shiny)
library(ieiprio)

`%||%` <- function(a, b) if (is.null(a)) b else a

tabla <- function(x) {
  x <- as.data.frame(tibble::as_tibble(unclass(x)[names(x)]))
  if (requireNamespace("DT", quietly = TRUE)) {
    DT::datatable(x, rownames = FALSE, class = "compact stripe hover",
                  options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  } else {
    x
  }
}
salida_tabla <- function(id) {
  if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput(id) else tableOutput(id)
}
# env se evalua aqui y no dentro del render, que corre en otro marco.
render_tabla <- function(expr, env = parent.frame()) {
  force(env)
  f <- substitute(expr)
  if (requireNamespace("DT", quietly = TRUE)) {
    DT::renderDT(tabla(eval(f, env)))
  } else {
    renderTable(as.data.frame(eval(f, env)))
  }
}
lineas <- function(txt) {
  v <- trimws(unlist(strsplit(txt %||% "", "[\n;,]")))
  v[nzchar(v)]
}
ayuda <- function(...) tags$p(class = "text-muted small mb-3", ...)

tema <- bslib::bs_theme(
  version = 5, primary = "#0b2545", secondary = "#475467", success = "#12b76a",
  "navbar-bg" = "#0b2545", "border-radius" = "0.5rem",
  base_font = bslib::font_collection("-apple-system", "Segoe UI", "Roboto", "Helvetica Neue", "Arial", "sans-serif")
)

barra <- function(...) bslib::sidebar(width = 320, bg = "#f8fafc", ...)

ui <- bslib::page_navbar(
  title = tags$strong("ieiprio"),
  theme = tema,
  window_title = "ieiprio",
  navbar_options = bslib::navbar_options(bg = "#0b2545", theme = "dark"),
  fillable = FALSE,
  header = tagList(
    tags$style(HTML(paste(
      ".navbar .form-switch .form-check-input{background-color:rgba(255,255,255,.2);border-color:rgba(255,255,255,.7);",
      "--bs-form-switch-bg:url(\"data:image/svg+xml,%3csvg xmlns='http://www.w3.org/2000/svg' viewBox='-4 -4 8 8'%3e%3ccircle r='3' fill='%23fff'/%3e%3c/svg%3e\")}",
      ".navbar .form-switch .form-check-input:checked{background-color:#12b76a;border-color:#12b76a}",
      ".navbar .form-switch label{color:#fff}"
    ))),
    tags$div(
      class = "alert alert-warning py-2 px-3 mb-0 rounded-0 small",
      "Herramienta de apoyo e investigaci\u00f3n, no de diagn\u00f3stico. ",
      "Las consultas en l\u00ednea env\u00edan solo coordenadas de variantes, nunca datos del paciente."
    )
  ),

  bslib::nav_panel("Fenotipo", icon = icon("stethoscope"),
    bslib::layout_sidebar(
      sidebar = barra(
        textAreaInput("fen", "Signos del caso", "agammaglobulinemia\nneumon\u00eda recurrente", rows = 5),
        ayuda("Uno por l\u00ednea, en espa\u00f1ol, ingl\u00e9s o c\u00f3digo HPO."),
        textAreaInput("aus", "Signos ausentes", "", rows = 2),
        selectInput("sim", "M\u00e9todo de similitud",
                    c("Por enfermedad" = "enfermedad", "Por gen" = "semantica", "Exacta" = "exacta")),
        sliderInput("n", "Genes a mostrar", 5, 50, 20, step = 5),
        radioButtons("vista_fen", "Gr\u00e1fico", inline = FALSE,
                     c("Ranking" = "ranking", "Mapa de similitud" = "mapa", "Cuerdas" = "acordes",
                       "Genoma en c\u00edrculo" = "circos", "Cromosomas" = "cariotipo")),
        actionButton("buscar", "Buscar genes", class = "btn-primary w-100", icon = icon("magnifying-glass"))
      ),
      bslib::card(
        bslib::card_header("Genes del panel ordenados por parecido con el caso"),
        ayuda("La similitud compara los signos del caso con la enfermedad m\u00e1s parecida de cada gen, de 0 a 1. ",
              "Las cuerdas unen cada signo con los genes que explica, y el genoma en c\u00edrculo ubica a los candidatos ",
              "con el color de su herencia y une los que comparten categor\u00eda IUIS."),
        plotOutput("plot_fen", height = "620px")
      ),
      bslib::card(bslib::card_header("Detalle"), salida_tabla("tab_fen"))
    )),

  bslib::nav_panel("Gen", icon = icon("dna"),
    bslib::layout_sidebar(
      sidebar = barra(
        textInput("gen", "Gen o alias", "BTK"),
        bslib::input_switch("gen_clinvar", "Consultar ClinVar y la prote\u00edna", TRUE),
        actionButton("ver_gen", "Ver ficha", class = "btn-primary w-100", icon = icon("file-lines"))
      ),
      bslib::card(bslib::card_header("Ficha del gen"), verbatimTextOutput("txt_gen")),
      bslib::card(
        bslib::card_header("Variantes de ClinVar sobre la prote\u00edna"),
        plotOutput("plot_gen", height = "380px")
      )
    )),

  bslib::nav_panel("Variante", icon = icon("pen"),
    bslib::layout_sidebar(
      sidebar = barra(
        textInput("var", "Variante", "BTK c.1632-2A>G"),
        ayuda("Gen y cambio en HGVS, o cromosoma, posici\u00f3n y alelos."),
        selectInput("var_gt", "Genotipo", c("No indicado" = "NA", "Heterocigota" = "het",
                                             "Homocigota" = "hom_alt", "Hemicigota" = "hemi")),
        textAreaInput("var_fen", "Fenotipo", "", rows = 3),
        actionButton("ver_var", "Interpretar", class = "btn-primary w-100", icon = icon("play"))
      ),
      bslib::card(bslib::card_header("Resultado"), salida_tabla("tab_var"))
    )),

  bslib::nav_panel("Archivo", icon = icon("folder-open"),
    bslib::layout_sidebar(
      sidebar = barra(
        fileInput("arch", "VCF, cromatogramas, FASTA, PDF o Word", multiple = TRUE,
                  buttonLabel = "Elegir", placeholder = "Ning\u00fan archivo"),
        textInput("arch_gen", "Gen secuenciado", ""),
        ayuda("Solo para cromatogramas Sanger."),
        selectInput("arch_sexo", "Sexo", c("No indicado" = "", "Masculino" = "M", "Femenino" = "F")),
        textAreaInput("arch_fen", "Fenotipo", "", rows = 3),
        actionButton("analizar", "Analizar", class = "btn-primary w-100", icon = icon("play")),
        tags$hr(),
        downloadButton("bajar_html", "Reporte HTML", class = "btn-outline-secondary w-100 mb-2"),
        downloadButton("bajar_xlsx", "Tabla en Excel", class = "btn-outline-secondary w-100")
      ),
      bslib::card(bslib::card_header("Resumen"), verbatimTextOutput("txt_arch")),
      bslib::card(bslib::card_header("Variantes priorizadas"), salida_tabla("tab_arch"))
    )),

  bslib::nav_panel("Validaci\u00f3n", icon = icon("chart-line"),
    bslib::layout_sidebar(
      sidebar = barra(
        sliderInput("val_n", "Casos simulados", 20, 300, 100, step = 20),
        sliderInput("val_imp", "Signos imprecisos", 0, 0.9, 0.3, step = 0.1),
        ayuda("Proporci\u00f3n de signos reemplazados por un t\u00e9rmino m\u00e1s general."),
        sliderInput("val_ruido", "Signos al azar por caso", 0, 5, 1),
        actionButton("validar", "Simular", class = "btn-primary w-100", icon = icon("play"))
      ),
      bslib::card(
        bslib::card_header("Qu\u00e9 tan bien encuentra el gen correcto cada m\u00e9todo"),
        ayuda("Se simulan pacientes a partir de enfermedades conocidas del panel y se mide en qu\u00e9 puesto ",
              "queda el gen verdadero. Una curva m\u00e1s alta y m\u00e1s a la izquierda es un m\u00e9todo mejor. ",
              "Sirve para comparar m\u00e9todos y ajustar par\u00e1metros, no para estimar el rendimiento con pacientes reales."),
        plotOutput("plot_val", height = "380px")
      ),
      bslib::card(bslib::card_header("Resumen"), tableOutput("tab_val"))
    )),

  bslib::nav_spacer(),
  bslib::nav_item(
    tags$div(class = "pt-2 text-white",
             bslib::input_switch("sin_red", "Sin conexi\u00f3n", isTRUE(getOption("ieiprio.sin_red"))))
  )
)

server <- function(input, output, session) {
  # Las opciones de la sesion de R se restauran al cerrar la aplicacion.
  previas <- options("ieiprio.sin_red", "ieiprio.silenciar_privacidad")
  session$onSessionEnded(function() options(previas))
  observe(options(ieiprio.sin_red = isTRUE(input$sin_red), ieiprio.silenciar_privacidad = TRUE))

  busqueda <- eventReactive(input$buscar, {
    aus <- lineas(input$aus)
    iei_buscar_fenotipo(lineas(input$fen), n = input$n, similitud = input$sim,
                        ausentes = if (length(aus)) aus)
  })
  output$plot_fen <- renderPlot({
    g <- plot(busqueda(), tipo = input$vista_fen, n = input$n)
    if (inherits(g, "ggplot")) print(g)
  }, res = 96)
  output$tab_fen <- render_tabla(busqueda())

  ficha <- eventReactive(input$ver_gen, {
    iei_gen(input$gen, clinvar = input$gen_clinvar, proteina = input$gen_clinvar)
  })
  output$txt_gen <- renderPrint(print(ficha()))
  output$plot_gen <- renderPlot({
    g <- ficha()
    validate(need(!is.null(g$clinvar), "Activa la consulta a ClinVar para ver las variantes sobre la prote\u00edna."))
    plot(g)
  }, res = 96)

  variante <- eventReactive(input$ver_var, {
    fen <- lineas(input$var_fen)
    gt <- if (identical(input$var_gt, "NA")) NA else input$var_gt
    iei_priorizar(iei_anotar(iei_variantes(input$var, genotipo = gt)), fenotipo = if (length(fen)) fen)
  })
  output$tab_var <- render_tabla(iei_tabla(variante()))

  resultado <- eventReactive(input$analizar, {
    req(input$arch)
    rutas <- file.path(dirname(input$arch$datapath[1]), input$arch$name)
    file.rename(input$arch$datapath, rutas)
    ext <- tolower(tools::file_ext(rutas[1]))
    sexo <- if (nzchar(input$arch_sexo)) input$arch_sexo
    fen <- lineas(input$arch_fen)
    v <- if (ext %in% c("vcf", "gz")) {
      iei_filtrar(iei_leer_vcf(rutas[1], sexo = sexo))
    } else if (ext %in% c("ab1", "abi", "scf", "fa", "fasta", "seq")) {
      iei_sanger(rutas, gen = if (nzchar(input$arch_gen)) input$arch_gen, sexo = sexo)
    } else {
      iei_variantes(rutas[1], sexo = sexo)
    }
    iei_priorizar(iei_anotar(v), fenotipo = if (length(fen)) fen)
  })
  output$txt_arch <- renderPrint(summary(resultado()))
  output$tab_arch <- render_tabla(iei_tabla(resultado()))
  output$bajar_html <- downloadHandler(
    filename = function() "ieiprio_reporte.html",
    content = function(file) {
      r <- resultado()
      iei_reporte(r, file, paciente = unique(r$muestra)[1])
    })
  output$bajar_xlsx <- downloadHandler(
    filename = function() "ieiprio_resultado.xlsx",
    content = function(file) iei_exportar(resultado(), file, formato = "xlsx"))

  validacion <- eventReactive(input$validar, {
    iei_validar_fenotipo(n = input$val_n, imprecision = input$val_imp, ruido = input$val_ruido)
  })
  output$plot_val <- renderPlot(plot(validacion()), res = 96)
  output$tab_val <- renderTable({
    r <- ieiprio:::.resumen_validacion(validacion())
    r$metodo <- unname(c(enfermedad = "Por enfermedad", semantica = "Por gen", exacta = "Exacta")[r$metodo])
    names(r) <- c("M\u00e9todo", "Primer puesto (%)", "Entre los 5 primeros (%)",
                  "Entre los 10 primeros (%)", "Puesto mediano")
    r
  }, digits = 0)
}

shinyApp(ui, server)
