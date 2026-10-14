# Genera un cromatograma ABIF sintetico para las pruebas, sin depender de
# archivos binarios en el repositorio. Cada base es un pico gaussiano; en las
# posiciones de `het` se agrega un segundo pico de la base indicada.
escribir_ab1 <- function(...) .escribir_ab1(...)

# Referencia sintetica reproducible y una funcion falsa para Ensembl.
referencia_prueba <- function(n = 3000, semilla = 11) {
  withr::with_seed(semilla, paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = ""))
}

revcomp_prueba <- function(x) {
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  paste(rev(comp[strsplit(x, "")[[1]]]), collapse = "")
}
