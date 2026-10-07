#' Generar el reporte HTML de un paciente
#'
#' Crea un archivo HTML autocontenido, sin dependencias externas, con las
#' variantes priorizadas de una muestra, una ficha por cada variante de
#' prioridad alta, la bitacora del filtrado, las versiones de todas las
#' fuentes de datos y las limitaciones del analisis. Se puede abrir en
#' cualquier navegador y enviar por correo.
#'
#' El reporte usa solo el codigo de la muestra. No incluyas nombres ni
#' datos que identifiquen al paciente.
#'
#' @param x Tabla devuelta por [iei_priorizar()].
#' @param archivo Ruta del HTML a crear.
#' @param paciente Codigo de la muestra. Obligatorio si `x` tiene varias.
#' @param responsable Nombre de quien firma el analisis (opcional).
#' @param panel Tabla de genes, por defecto [iuis_panel].
#' @param dominios Si `TRUE`, consulta en Ensembl el largo y los dominios
#'   de cada proteina para dibujar el diagrama de la variante. Usa cache.
#' @param cache Carpeta de la cache, o `FALSE` para no usarla.
#'
#' @return La ruta del archivo, de forma invisible.
#' @export
#' @examples
#' \dontrun{
#' iei_reporte(prior, "P01.html", paciente = "P01")
#' }
iei_reporte <- function(x, archivo, paciente = NULL, responsable = NULL,
                        panel = ieiprio::iuis_panel, dominios = TRUE,
                        cache = tools::R_user_dir("ieiprio", "cache")) {
  necesarias <- c("muestra", "gen", "puntaje", "categoria", "razones")
  faltan <- setdiff(necesarias, names(x))
  if (length(faltan) > 0) {
    cli::cli_abort("{.arg x} debe venir de {.fn iei_priorizar}. Faltan: {.val {faltan}}.")
  }
  muestras <- unique(x$muestra)
  if (is.null(paciente)) {
    if (length(muestras) > 1) {
      cli::cli_abort(c("{.arg x} tiene {length(muestras)} muestras; indica {.arg paciente}.",
                       "i" = "Disponibles: {.val {muestras}}."))
    }
    paciente <- muestras %||% "sin_muestra"
    if (length(paciente) == 0) paciente <- "sin_muestra"
  } else if (!paciente %in% muestras && nrow(x) > 0) {
    cli::cli_abort("La muestra {.val {paciente}} no esta en {.arg x}.")
  }

  datos <- x[x$muestra == paciente, ]
  a <- .attrs_ieiprio(x)
  build <- a$build %||% NA_character_
  fichas <- datos[datos$categoria %in% c("alta", "media"), ]
  proteinas <- if (dominios && "proteina_id" %in% names(fichas) && nrow(fichas) > 0) {
    .proteinas_info(utils::head(fichas$proteina_id, 12), build, cache)
  } else list()
  html <- .html_reporte(datos, paciente, responsable, panel, a, proteinas)
  dir.create(dirname(archivo), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, archivo, useBytes = TRUE)
  cli::cli_inform("Reporte de {.val {paciente}} guardado en {.path {archivo}}.")
  invisible(archivo)
}

.esc <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

.version_paquete <- function() {
  tryCatch(as.character(utils::packageVersion("ieiprio")), error = function(e) "desarrollo")
}

.enlace <- function(url, texto) {
  ifelse(is.na(url) | url == "", "",
         sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', .esc(url), .esc(texto)))
}

.enlaces_variante <- function(d, build) {
  dataset <- if (identical(build, "GRCh37")) "gnomad_r2_1" else "gnomad_r4"
  gnomad <- if (all(c("chr", "pos", "ref", "alt") %in% names(d))) {
    sprintf("https://gnomad.broadinstitute.org/variant/%s-%s-%s-%s?dataset=%s",
            d$chr, d$pos, d$ref, d$alt, dataset)
  } else rep(NA_character_, nrow(d))
  rs <- if ("id" %in% names(d)) d$id else rep(NA_character_, nrow(d))
  clinvar <- ifelse(!is.na(rs) & startsWith(rs, "rs"),
                    paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=", rs), NA_character_)
  panelapp <- paste0("https://panelapp.genomicsengland.co.uk/panels/398/", d$gen, "/")
  paste(.enlace(panelapp, "PanelApp"), .enlace(gnomad, "gnomAD"), .enlace(clinvar, "ClinVar"),
        sep = " ")
}

.celda <- function(d, col) if (col %in% names(d)) d[[col]] else rep(NA, nrow(d))

.formato_af <- function(af) ifelse(is.na(af), "ausente", formatC(af, format = "e", digits = 1))

.tabla_html <- function(encabezados, filas, clase = NULL, num = integer(0)) {
  clases <- ifelse(seq_along(encabezados) %in% num, ' class="num"', "")
  th <- paste0("<th", clases, ">", encabezados, "</th>", collapse = "")
  tr <- if (length(filas) == 0) {
    sprintf('<tr><td colspan="%d" class="vacio">Sin variantes</td></tr>', length(encabezados))
  } else filas
  sprintf('<div class="tabla"><table%s><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>',
          if (is.null(clase)) "" else sprintf(' class="%s"', clase), th, paste(tr, collapse = ""))
}

.filas_variantes <- function(d) {
  if (nrow(d) == 0) return(character(0))
  variante <- ifelse(is.na(.celda(d, "hgvs_c")),
                     sprintf("%s:%s %s&gt;%s", .esc(.celda(d, "chr")), .esc(.celda(d, "pos")),
                             .esc(.celda(d, "ref")), .esc(.celda(d, "alt"))),
                     .esc(sub("^[^:]+:", "", .celda(d, "hgvs_c"))))
  proteina <- .esc(sub("^[^:]+:", "", .celda(d, "hgvs_p")))
  sprintf(paste0('<tr class="%s"><td><b>%s</b></td><td>%s<br><span class="sub">%s</span></td>',
                 '<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class="acmg">%s</td><td class="num">%s</td>',
                 '<td><span class="etiqueta %s">%s</span></td><td>%s</td></tr>'),
          .esc(d$categoria), .esc(d$gen), variante, proteina,
          .esc(.celda(d, "genotipo")), .esc(gsub("_", " ", .celda(d, "consecuencia"))),
          .esc(gsub("_", " ", .celda(d, "clinvar"))), .formato_af(.celda(d, "af_gnomad")),
          .esc(.celda(d, "acmg")), .esc(d$puntaje), .esc(d$categoria), .esc(d$categoria), .esc(d$razones))
}

.ficha <- function(r, info, build, svg) {
  enlaces <- .enlaces_variante(r, build)
  sprintf(paste0('<div class="ficha %s"><h3>%s <span class="sub">%s %s</span> ',
                 '<span class="etiqueta %s">%s &middot; %s</span></h3>%s',
                 '<p><b>Genotipo</b> %s &middot; <b>Herencia del gen</b> %s &middot; ',
                 '<b>Consecuencia</b> %s &middot; <b>gnomAD</b> %s</p>',
                 '<p><b>ACMG sugerido</b> %s</p>',
                 '<p><b>Razones</b> %s</p>%s',
                 '<p><b>Fenotipos asociados al gen (PanelApp)</b> %s</p>',
                 '<p class="enlaces">%s</p></div>'),
          .esc(r$categoria), .esc(r$gen), .esc(sub("^[^:]+:", "", .celda(r, "hgvs_c"))),
          .esc(sub("^[^:]+:", "", .celda(r, "hgvs_p"))),
          .esc(r$categoria), .esc(r$categoria), .esc(r$puntaje), svg,
          .esc(.celda(r, "genotipo")), .esc(info$herencia),
          .esc(gsub("_", " ", .celda(r, "consecuencia"))), .formato_af(.celda(r, "af_gnomad")),
          ifelse(is.na(.celda(r, "acmg")), "ninguno", .esc(.celda(r, "acmg"))),
          .esc(r$razones),
          if (is.na(.celda(r, "nota"))) "" else sprintf('<p class="nota">%s</p>', .esc(r$nota)),
          .esc(info$fenotipos), enlaces)
}

.caja_sexo <- function(qc) {
  if (is.null(qc)) return("")
  dec <- qc$sexo_declarado; inf <- qc$sexo_inferido
  txt_inf <- if (is.na(inf)) sprintf("no estimable (%d variantes en X)", qc$n_x) else
    sprintf("%s (%.0f%% heterocigotas en %d variantes del X)", inf, 100 * qc$frac_het_x, qc$n_x)
  txt_dec <- if (is.na(dec)) "no indicado" else dec
  clase <- if (is.na(dec) || is.na(inf) || inf == "indeterminado") "neutro" else if (dec == inf) "bien" else "mal"
  mensaje <- switch(clase,
    bien = "El sexo declarado coincide con el estimado desde el cromosoma X.",
    mal = "El sexo declarado NO coincide con el estimado. Revisa la identidad de la muestra antes de interpretar genes ligados al X.",
    neutro = if (is.na(dec)) "Sin sexo declarado, las variantes del X no se marcan como hemicigotas. Indica el sexo en iei_leer_vcf()." else
      "No hay suficientes variantes en el X para confirmar el sexo declarado.")
  sprintf('<div class="sexo %s"><b>Verificaci&oacute;n de sexo</b> Declarado %s &middot; Estimado %s. %s</div>',
          clase, .esc(txt_dec), .esc(txt_inf), .esc(mensaje))
}

.bloque_qc <- function(qc) {
  if (is.null(qc)) return("<p class=\"sub\">Sin datos de control de calidad (el reporte no viene de iei_filtrar()).</p>")
  ratio <- if (qc$n_hom > 0) qc$n_het / qc$n_hom else NA
  metricas <- .tabla_html(
    c("Variantes llamadas", "Sin llamar", "Het / hom", "Ti/Tv", "DP mediana", "GQ mediana"),
    sprintf('<tr><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td></tr>',
            format(qc$n_variantes, big.mark = " "), format(qc$sin_llamar, big.mark = " "),
            ifelse(is.na(ratio), "-", sprintf("%.2f", ratio)),
            ifelse(is.na(qc$titv), "-", sprintf("%.2f", qc$titv)),
            ifelse(is.na(qc$mediana_dp), "-", .num(qc$mediana_dp)),
            ifelse(is.na(qc$mediana_gq), "-", .num(qc$mediana_gq))),
    num = 1:6)
  paste0(.caja_sexo(qc), metricas,
         '<p class="sub">Referencias habituales: Ti/Tv cerca de 2.0 en genomas y de 2.8 a 3.0 en exomas; relaci&oacute;n het/hom entre 1.5 y 2.0. Valores muy distintos sugieren problemas de calidad o de llamado.</p>',
         '<div class="grafs">',
         .svg_histograma(qc$hist_dp, "Profundidad (DP)", "lecturas (100 o m&aacute;s al final)", "#0ba5ec"),
         .svg_histograma(qc$hist_gq, "Calidad del genotipo (GQ)", "GQ", "#7a5af8"),
         .svg_histograma(qc$hist_ab, "Fracci&oacute;n al&eacute;lica en heterocigotos", "fracci&oacute;n del alelo alternativo", "#12b76a"),
         "</div>")
}

.bloque_parametros <- function(a, release) {
  vcf <- a$ieiprio_vcf; p <- a$ieiprio_parametros; pesos <- a$ieiprio_pesos
  maxaf <- a$ieiprio_max_af
  fila <- function(k, v) sprintf("<tr><td>%s</td><td>%s</td></tr>", k, v)
  filas <- c(
    fila("Archivo VCF", .esc(vcf$archivo %||% "-")),
    fila("Huella md5 del VCF", sprintf("<code>%s</code>", .esc(vcf$md5 %||% "-"))),
    fila("Genoma", .esc(a$build %||% "-")),
    fila("Sexo declarado", .esc(paste(sprintf("%s %s", names(a$ieiprio_sexo), ifelse(is.na(a$ieiprio_sexo), "no indicado", a$ieiprio_sexo)), collapse = "; "))),
    if (!is.null(p)) c(
      fila("Genes analizados", .esc(p$genes)),
      fila("Filtros de calidad", sprintf("DP &ge; %s, GQ &ge; %s, fracci&oacute;n al&eacute;lica %s a %s en heterocigotos, %s, margen %s pb",
                                        p$min_dp, p$min_gq, p$min_ab, 1 - p$min_ab,
                                        if (isTRUE(p$solo_pass)) "solo FILTER PASS" else "cualquier FILTER", p$margen))),
    if (!is.null(maxaf)) fila("Frecuencia m&aacute;xima en gnomAD",
                              .esc(paste(sprintf("%s %s", names(maxaf), maxaf), collapse = "; "))),
    if (!is.null(pesos)) fila("Pesos y cortes", .esc(paste(sprintf("%s=%s", names(pesos), unlist(pesos)), collapse = ", "))),
    fila("Ensembl VEP", sprintf("release %s", .esc(release)))
  )
  .tabla_html(c("Par&aacute;metro", "Valor"), filas, clase = "param")
}

.html_reporte <- function(d, paciente, responsable, panel, a, proteinas) {
  build <- a$build %||% NA_character_
  release <- a$vep_release %||% NA_character_
  fenotipo <- a$ieiprio_fenotipo; categoria <- a$ieiprio_categoria
  qc <- a$ieiprio_qc[[paciente]]
  fuentes <- ieiprio::fuentes_datos
  n_cat <- table(factor(d$categoria, levels = c("alta", "media", "baja")))
  principales <- d[d$categoria %in% c("alta", "media"), ]
  bajas <- d[d$categoria == "baja", ]

  resumen <- sprintf(paste0('<div class="resumen"><div class="caja alta"><span>%d</span>alta</div>',
                            '<div class="caja media"><span>%d</span>media</div>',
                            '<div class="caja baja"><span>%d</span>baja</div></div>'),
                     n_cat[["alta"]], n_cat[["media"]], n_cat[["baja"]])
  destacado <- if (nrow(d) > 0) {
    sprintf('<p>Variante con mayor puntaje: <b>%s</b> %s (%s, puntaje %s).</p>',
            .esc(d$gen[1]), .esc(sub("^[^:]+:", "", .celda(d, "hgvs_p")[1])),
            .esc(d$categoria[1]), .esc(d$puntaje[1]))
  } else "<p>No quedaron variantes despu&eacute;s del filtrado.</p>"
  alerta_sexo <- if (!is.null(qc) && !is.na(qc$sexo_declarado) && !is.na(qc$sexo_inferido) &&
                     qc$sexo_inferido != "indeterminado" && qc$sexo_declarado != qc$sexo_inferido) .caja_sexo(qc) else ""

  encab <- c("Gen", "Variante", "Genotipo", "Consecuencia", "ClinVar", "gnomAD",
             "ACMG sugerido", "Puntaje", "Prioridad", "Razones")
  tabla <- .tabla_html(encab, .filas_variantes(principales), num = 8)
  tabla_baja <- sprintf('<details><summary>Variantes de prioridad baja (%d)</summary>%s</details>',
                        nrow(bajas), .tabla_html(encab, .filas_variantes(bajas), num = 8))

  fichas <- if (nrow(principales) == 0) "<p>No hay variantes de prioridad alta ni media.</p>" else {
    paste(vapply(seq_len(min(nrow(principales), 12)), function(i) {
      r <- principales[i, ]
      info <- if ("proteina_id" %in% names(r) && !is.na(r$proteina_id)) proteinas[[r$proteina_id]] else NULL
      marcas <- if (!is.null(info)) {
        mismas <- d[!is.na(.celda(d, "proteina_id")) & d$proteina_id == r$proteina_id, ]
        data.frame(pos = mismas$pos_proteina,
                   etiqueta = sub("^[^:]+:p\\.", "", .celda(mismas, "hgvs_p")),
                   categoria = mismas$categoria, stringsAsFactors = FALSE)
      } else NULL
      .ficha(r, panel[match(r$gen, panel$gen), ], build,
             if (is.null(info)) "" else .svg_proteina(info, marcas))
    }, character(1)), collapse = "")
  }

  notas <- d[!is.na(.celda(d, "nota")) & !d$categoria %in% c("alta", "media"), ]
  bloque_notas <- if (nrow(notas) == 0) "" else sprintf(
    "<h2>Notas de variantes de prioridad baja</h2><ul>%s</ul>",
    paste(sprintf("<li><b>%s</b>%s %s</li>", .esc(notas$gen),
                  ifelse(is.na(.celda(notas, "hgvs_p")), "",
                         paste0(" ", .esc(sub("^[^:]+:", "", .celda(notas, "hgvs_p"))), ".")),
                  .esc(notas$nota)), collapse = ""))

  bitacora <- a$ieiprio_bitacora
  bloque_bitacora <- if (is.null(bitacora)) "<p>Sin bit&aacute;cora.</p>" else paste0(
    .svg_embudo(bitacora),
    "<details><summary>Ver tabla</summary>",
    .tabla_html(c("Paso", "Filas", "Descartadas"),
                sprintf('<tr><td>%s</td><td class="num">%s</td><td class="num">%s</td></tr>',
                        .esc(bitacora$paso), .esc(bitacora$filas), .esc(bitacora$descartadas)),
                num = 2:3), "</details>",
    '<p class="sub">Cada fila es una variante en una muestra del VCF, antes de separar por paciente.</p>')

  bloque_fuentes <- .tabla_html(
    c("Fuente", "Versi&oacute;n", "Fecha", "md5"),
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td><code>%s</code></td></tr>",
            .enlace(fuentes$url, fuentes$fuente), .esc(fuentes$version),
            .esc(fuentes$fecha), .esc(fuentes$md5)))

  fenotipo_txt <- if (!is.null(fenotipo)) {
    nombres <- ieiprio::hpo_genes$hpo_nombre[match(fenotipo, ieiprio::hpo_genes$hpo_id)]
    paste(sprintf("%s%s", .esc(fenotipo), ifelse(is.na(nombres), "", paste0(" ", .esc(nombres)))),
          collapse = "; ")
  } else if (!is.null(categoria)) {
    paste("Categor&iacute;a IUIS", .esc(paste(categoria, collapse = ", ")))
  } else "No indicado"
  bloque_fenotipo <- if (is.null(fenotipo) || nrow(principales) == 0) "" else paste0(
    "<h2>Coincidencia con el fenotipo</h2>",
    '<p class="sub">Genes de prioridad alta y media frente a los t&eacute;rminos HPO del paciente. Un punto indica que el gen tiene ese t&eacute;rmino anotado en HPO.</p>',
    .tabla_fenotipo(principales, fenotipo, ieiprio::hpo_genes))

  paste0(
    '<!DOCTYPE html><html lang="es"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    sprintf("<title>ieiprio &middot; %s</title>", .esc(paciente)),
    "<style>", .css_reporte(), "</style></head><body><main>",
    '<header><p class="marca">ieiprio</p>',
    sprintf("<h1>Reporte de variantes &middot; %s</h1>", .esc(paciente)),
    sprintf('<p class="sub">Generado el %s &middot; ieiprio %s &middot; genoma %s &middot; VCF %s%s</p>',
            format(Sys.Date()), .version_paquete(), .esc(build),
            .esc(a$ieiprio_vcf$archivo %||% "-"),
            if (is.null(responsable)) "" else paste0(" &middot; ", .esc(responsable))),
    sprintf("<p><b>Fenotipo usado en la priorizaci&oacute;n</b> %s</p></header>", fenotipo_txt),
    '<div class="aviso">Herramienta de apoyo e investigaci&oacute;n. No es un informe diagn&oacute;stico ',
    "ni reemplaza la clasificaci&oacute;n ACMG/AMP realizada por un especialista. ",
    "Los criterios ACMG mostrados son sugerencias autom&aacute;ticas. ",
    "Toda variante relevante debe confirmarse por un m&eacute;todo independiente.</div>",
    "<h2>Resumen</h2>", alerta_sexo, resumen, destacado,
    "<h2>Control de calidad de la muestra</h2>", .bloque_qc(qc),
    "<h2>Variantes priorizadas</h2>", tabla, .svg_puntaje(principales), tabla_baja,
    "<h2>Fichas de variantes</h2>", fichas,
    bloque_fenotipo, bloque_notas,
    "<h2>Filtrado</h2>", bloque_bitacora,
    "<h2>Par&aacute;metros y trazabilidad</h2>", .bloque_parametros(a, release),
    "<h3>Fuentes de datos</h3>", bloque_fuentes,
    "<h2>Limitaciones</h2><ul>",
    "<li>Solo se analizan variantes puntuales e indels peque&ntilde;os del VCF. No se detectan CNV, variantes estructurales ni expansiones de repetidos.</li>",
    "<li>No se eval&uacute;a la cobertura de los genes del panel, porque requiere el archivo BAM. Una regi&oacute;n sin cobertura no aparece como variante, pero tampoco queda descartada.</li>",
    "<li>Los criterios ACMG sugeridos se infieren solo de consecuencia, frecuencia, fenotipo y genotipo. PVS1 requiere confirmar que la p&eacute;rdida de funci&oacute;n es el mecanismo de enfermedad y que el transcrito es relevante.</li>",
    "<li>La coincidencia de fenotipo compara t&eacute;rminos HPO exactos, sin usar la jerarqu&iacute;a de la ontolog&iacute;a.</li>",
    "<li>Los heterocigotos compuestos son posibles mientras no se confirme la fase estudiando a los padres.</li>",
    "<li>La verificaci&oacute;n de sexo necesita al menos 20 variantes en el cromosoma X; en paneles peque&ntilde;os puede no ser estimable.</li>",
    "<li>El panel de genes y las anotaciones dependen de las versiones listadas arriba.</li></ul>",
    "</main></body></html>"
  )
}

.css_reporte <- function() {
  paste0(
    ":root{--tinta:#1f2328;--suave:#59636e;--borde:#d8dee4;--fondo:#ffffff;--panel:#f6f8fa;",
    "--alta:#b42318;--alta-f:#fef3f2;--media:#b54708;--media-f:#fffaeb;--baja:#475467;--baja-f:#f2f4f7;--enlace:#175cd3}",
    "*{box-sizing:border-box}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;",
    "color:var(--tinta);background:var(--fondo);line-height:1.5;font-size:15px}",
    "main{max-width:1100px;margin:0 auto;padding:32px 20px 64px}",
    ".marca{font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--suave);font-size:12px;margin:0}",
    "h1{font-size:26px;margin:4px 0 4px}h2{font-size:18px;margin:36px 0 12px;padding-bottom:6px;border-bottom:1px solid var(--borde)}",
    "h3{font-size:16px;margin:0 0 8px}.sub{color:var(--suave);font-size:13px;font-weight:400}",
    ".aviso{background:var(--panel);border-left:4px solid var(--suave);padding:12px 16px;margin:20px 0;font-size:14px}",
    ".resumen{display:flex;gap:12px;margin:8px 0 12px}.caja{flex:1;padding:12px 16px;border-radius:8px;font-size:13px;text-transform:uppercase;letter-spacing:.05em}",
    ".caja span{display:block;font-size:28px;font-weight:700;letter-spacing:0}",
    ".caja.alta{background:var(--alta-f);color:var(--alta)}.caja.media{background:var(--media-f);color:var(--media)}.caja.baja{background:var(--baja-f);color:var(--baja)}",
    ".tabla{overflow-x:auto;border:1px solid var(--borde);border-radius:8px;margin:8px 0}",
    "table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px 10px;text-align:left;vertical-align:top;border-bottom:1px solid var(--borde)}",
    "th{background:var(--panel);font-weight:600;white-space:nowrap}tr:last-child td{border-bottom:0}td.num{text-align:right;white-space:nowrap}",
    "td.vacio{color:var(--suave);text-align:center}",
    ".etiqueta{display:inline-block;padding:1px 8px;border-radius:99px;font-size:12px;font-weight:600}",
    ".etiqueta.alta{background:var(--alta-f);color:var(--alta)}.etiqueta.media{background:var(--media-f);color:var(--media)}.etiqueta.baja{background:var(--baja-f);color:var(--baja)}",
    ".ficha{border:1px solid var(--borde);border-left:4px solid var(--alta);border-radius:8px;padding:14px 18px;margin:12px 0}",
    ".ficha p{margin:6px 0}.nota{background:var(--media-f);padding:6px 10px;border-radius:6px}",
    ".enlaces a{margin-right:12px}a{color:var(--enlace)}details{margin:8px 0}summary{cursor:pointer;color:var(--suave)}",
    "code{font-size:12px}",
    ".ficha.media{border-left-color:var(--media)}.ficha h3 .etiqueta{margin-left:6px;vertical-align:2px}",
    ".sexo{padding:10px 14px;border-radius:8px;margin:10px 0;font-size:14px}.sexo.bien{background:#ecfdf3;color:#067647}",
    ".sexo.mal{background:var(--alta-f);color:var(--alta);font-weight:600}.sexo.neutro{background:var(--panel);color:var(--suave)}",
    ".grafs{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:12px;margin:12px 0}",
    ".graf{margin:0}.graf figcaption{font-size:13px;font-weight:600;margin-bottom:4px}.graf.ancho figcaption{font-weight:400;color:var(--suave);font-size:12px}",
    ".svg{width:100%;height:auto;font-size:11px;fill:var(--tinta)}.svg .eje{stroke:var(--borde)}.svg .cero{stroke:var(--tinta);stroke-dasharray:3 3}",
    ".svg .eti{fill:var(--suave)}.svg .gen{font-size:12px}.svg .total{font-size:12px;font-weight:600}.svg .sub{fill:var(--suave);font-weight:400}",
    ".svg .dom{font-size:10px;fill:#344054}.svg .pin{font-size:11px;font-weight:600}.prot{margin:4px 0 6px}",
    "th.num{text-align:right}h3{margin-top:24px}td.acmg{font-size:12px;white-space:nowrap}.matriz td.si{color:var(--alta);text-align:center;font-size:16px}.matriz td.no{color:var(--borde);text-align:center}",
    ".param td:first-child{white-space:nowrap;color:var(--suave);width:220px}",
    "@media print{.tabla{overflow:visible}details{display:block}details>summary{display:none}.ficha{break-inside:avoid}}"
  )
}
