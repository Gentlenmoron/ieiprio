.o <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.chr1 <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  x <- x[[1]]
  if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
}

.agente <- function(req) {
  httr2::req_user_agent(req, "ieiprio (https://github.com/Gentlenmoron/ieiprio)")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Atributos que el paquete arrastra de una etapa a otra.
.attrs_ieiprio <- function(x) {
  a <- attributes(x)
  a[names(a) %in% c("build", "vep_release", "archivo") | startsWith(names(a), "ieiprio_")]
}

.poner_attrs <- function(x, a) {
  for (n in names(a)) attr(x, n) <- a[[n]]
  x
}
