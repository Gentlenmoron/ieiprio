.hpo_url <- function(tag) {
  sprintf(paste0("https://github.com/obophenotype/human-phenotype-ontology/",
                 "releases/download/%s/genes_to_phenotype.txt"), tag)
}

.hpo_descargar <- function(tag, destino) {
  httr2::request(.hpo_url(tag)) |>
    .agente() |>
    httr2::req_timeout(300) |>
    httr2::req_retry(max_tries = 4) |>
    httr2::req_perform(path = destino)
  invisible(destino)
}

# Lee genes_to_phenotype.txt y se queda con los genes indicados.
.hpo_procesar <- function(ruta, genes) {
  g2p <- utils::read.delim(ruta, quote = "", stringsAsFactors = FALSE,
                           colClasses = "character")
  necesarias <- c("gene_symbol", "hpo_id", "hpo_name")
  if (!all(necesarias %in% names(g2p))) {
    cli::cli_abort("El archivo de HPO no tiene las columnas {.val {necesarias}}.")
  }
  g2p <- g2p[g2p$gene_symbol %in% genes, necesarias]
  g2p <- unique(g2p)
  g2p <- g2p[order(g2p$gene_symbol, g2p$hpo_id), ]
  tibble::tibble(gen = g2p$gene_symbol, hpo_id = g2p$hpo_id,
                 hpo_nombre = g2p$hpo_name)
}
