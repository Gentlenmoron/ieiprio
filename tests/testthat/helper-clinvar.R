# Respuestas simuladas de ClinVar (E-utilities) y Ensembl para iei_gen().
clinvar_falso <- function(endpoint, ...) {
  q <- list(...)
  if (endpoint == "esearch.fcgi") {
    if (identical(q$retmax, 0)) {
      n <- c(clinsig_pathogenic = 120, clinsig_likely_pathogenic = 40, clinsig_vus = 300,
             clinsig_likely_benign = 90, clinsig_benign = 60)
      cl <- sub(".*AND ([a-z_]+)\\[Properties\\].*", "\\1", q$term)
      return(list(esearchresult = list(count = as.character(n[[cl]]))))
    }
    return(list(esearchresult = list(idlist = list("101", "102", "103", "104", "105"))))
  }
  v <- function(uid, titulo, clas, rev, cons, cambio, gen = "BTK") {
    list(uid = uid, title = titulo, protein_change = cambio,
         molecular_consequence_list = as.list(cons),
         germline_classification = list(description = clas, review_status = rev),
         genes = list(list(symbol = gen)))
  }
  list(result = list(
    uids = list("101", "102", "103", "104", "105"),
    `101` = v("101", "NM_000061.3(BTK):c.1630A>G (p.Arg544Gly)", "Pathogenic",
              "criteria provided, multiple submitters, no conflicts", "missense variant", "R544G"),
    `102` = v("102", "NM_000061.3(BTK):c.1632-2A>G", "Pathogenic", "criteria provided, single submitter",
              "splice acceptor variant", ""),
    `103` = v("103", "NM_000061.3(BTK):c.862C>T (p.Arg288Trp)", "Likely pathogenic",
              "reviewed by expert panel", "missense variant", "R288W"),
    `104` = v("104", "NM_000061.3(BTK):c.1581C>T (p.Gln527Ter)", "Pathogenic",
              "criteria provided, single submitter", "nonsense", "Q527*"),
    `105` = v("105", "NM_000000.1(OTRO):c.1A>G (p.Met1Val)", "Pathogenic",
              "no assertion criteria provided", "missense variant", "M1V", gen = "OTRO")
  ))
}

proteina_falsa <- function(id, build) {
  list(largo = 659, dominios = data.frame(inicio = c(6, 280, 402, 400 + 0), fin = c(133, 330, 500, 400),
                                          nombre = c("PH", "SH3", "SH2", "x"))[1:3, ])
}
