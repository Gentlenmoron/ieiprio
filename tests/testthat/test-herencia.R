test_that(".mapear_herencia traduce los textos de PanelApp", {
  entrada <- c(
    "MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted",
    "BIALLELIC, autosomal or pseudoautosomal",
    "BOTH monoallelic and biallelic (but BIALLELIC mutations cause a more SEVERE disease form), autosomal or pseudoautosomal",
    "X-LINKED: hemizygous mutation in males, biallelic mutations in females",
    "MITOCHONDRIAL", "Unknown", "", NA
  )
  expect_equal(.mapear_herencia(entrada),
               c("AD", "AR", "AR/AD", "XL", "MT", NA, NA, NA))
})

test_that(".mapear_evidencia usa los niveles de PanelApp", {
  expect_equal(.mapear_evidencia(c("3", "2", "1", "0", NA)),
               c("verde", "ambar", "rojo", "rojo", NA))
})
