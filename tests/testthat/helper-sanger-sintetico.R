# Cromatogramas sinteticos: picos gaussianos para uno o dos alelos.
cromatograma_sintetico <- function(alelo1, alelo2 = alelo1, alto2 = 0.9, paso = 10L, archivo = "sintetico") {
  n <- min(nchar(alelo1), nchar(alelo2))
  b1 <- strsplit(substr(alelo1, 1, n), "")[[1]]
  b2 <- strsplit(substr(alelo2, 1, n), "")[[1]]
  largo <- (n + 1L) * paso
  tr <- matrix(0, largo, 4, dimnames = list(NULL, c("A", "C", "G", "T")))
  picos <- paso * seq_len(n)
  forma <- exp(-((-6:6)^2) / 8)
  for (i in seq_len(n)) {
    r <- picos[i] + (-6:6)
    tr[r, b1[i]] <- tr[r, b1[i]] + 1000 * forma
    tr[r, b2[i]] <- tr[r, b2[i]] + 1000 * alto2 * forma
  }
  structure(list(bases = paste(b1, collapse = ""), calidad = rep(40L, n), picos = as.integer(picos),
                 trazas = tibble::as_tibble(as.data.frame(round(tr))), muestra = NA_character_,
                 archivo = archivo), class = "iei_ab1")
}

ref_sintetica <- function(n = 400, semilla = 7) {
  withr::with_seed(semilla, paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = ""))
}
