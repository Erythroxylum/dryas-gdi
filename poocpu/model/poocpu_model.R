# Exact PoOcPu1 submodel used for the two Eurasian gdi comparisons.
# Topology supplied by Dawson White from the H4d-derived model diagram:
# (Po,((octo_EU,octo_Carp_MK)Octo,(Pu,RU_SJ)Punc)AOcPu)R;

population_order <- c("Po", "octo_EU", "octo_Carp_MK", "Pu", "RU_SJ")

comparison_specs <- tibble::tribble(
  ~comparison, ~config, ~pop_a, ~pop_b, ~focal_population, ~tau_node,
  "octo_EU_vs_MK", "aab", "octo_EU", "octo_Carp_MK", "octo_EU",      "Octo",
  "octo_EU_vs_MK", "abb", "octo_EU", "octo_Carp_MK", "octo_Carp_MK", "Octo",
  "Pu_vs_RUSJ",    "aab", "Pu",      "RU_SJ",         "Pu",           "Punc",
  "Pu_vs_RUSJ",    "abb", "Pu",      "RU_SJ",         "RU_SJ",        "Punc"
)

# Migration edges retained from PoOcPu1 / H4d.
# BPP direction follows the fitted-model specification supplied with the analysis.
migration_pairs <- list(
  c("octo_EU", "RU_SJ"),
  c("RU_SJ", "octo_EU"),
  c("Po", "Punc")
)

migration_parameter_names <- c(
  "W_octo_EU_to_RU_SJ" = "W_EU_to_RUSJ",
  "W_RU_SJ_to_octo_EU" = "W_RUSJ_to_EU",
  "W_Po_to_Punc"        = "W_Po_to_Punc"
)

make_species_tree <- function(row) {
  sprintf(
    "(Po #%g, ((octo_EU #%g, octo_Carp_MK #%g)Octo:%g #%g, (Pu #%g, RU_SJ #%g)Punc:%g #%g)AOcPu:%g #%g)R:%g #%g;",
    row$theta_Po,
    row$theta_octo_EU,
    row$theta_octo_Carp_MK,
    row$tau_Octo,
    row$theta_Octo,
    row$theta_Pu,
    row$theta_RU_SJ,
    row$tau_Punc,
    row$theta_Punc,
    row$tau_AOcPu,
    row$theta_AOcPu,
    row$tau_R,
    row$theta_R
  )
}

required_parameter_columns <- c(
  "chromosome",
  "tau_Octo", "tau_Punc", "tau_AOcPu", "tau_R",
  "theta_octo_EU", "theta_octo_Carp_MK", "theta_Pu", "theta_RU_SJ", "theta_Po",
  "theta_Octo", "theta_Punc", "theta_AOcPu", "theta_R",
  "W_EU_to_RUSJ", "W_RUSJ_to_EU", "W_Po_to_Punc"
)
