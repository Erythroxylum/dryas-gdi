# Four-population ajanensis/alaskensis demographic model used by Yuttapong.
#
# Population tree:
#   ((ajan_Interior, ajan_Seward)J,
#    (alas_Interior, alas_Seward)L)R;
#
# The fitted model includes four migration edges between geographically
# corresponding Interior and Seward populations across the two species.

# Keep the exact population order used in the original BPP analysis.
population_order <- c(
  "ajan_Interior",
  "ajan_Seward",
  "alas_Interior",
  "alas_Seward"
)

# Test the two within-species population pairs separately using three-sequence
# aab/abb simulations. This reproduces Yuttapong's biological hypotheses while
# using the same Eq. 13 parser as the other analyses in this repository.
comparison_specs <- tibble::tribble(
  ~comparison,        ~config, ~pop_a,          ~pop_b,          ~focal_population, ~tau_node,
  "ajan_I_vs_S",      "aab",   "ajan_Interior", "ajan_Seward",   "ajan_Interior",   "J",
  "ajan_I_vs_S",      "abb",   "ajan_Interior", "ajan_Seward",   "ajan_Seward",     "J",
  "alas_I_vs_S",      "aab",   "alas_Interior", "alas_Seward",   "alas_Interior",   "L",
  "alas_I_vs_S",      "abb",   "alas_Interior", "alas_Seward",   "alas_Seward",     "L"
)

# Retain all four migration edges from the fitted four-population MSC-M model.
migration_pairs <- list(
  c("alas_Seward",   "ajan_Seward"),
  c("ajan_Seward",   "alas_Seward"),
  c("alas_Interior", "ajan_Interior"),
  c("ajan_Interior", "alas_Interior")
)

# Map each migration edge to the compact column name used in our extracted
# chromosome-specific parameter table.
migration_parameter_names <- c(
  "W_alas_Seward_to_ajan_Seward" = "W_alas_Seward_to_ajan_Seward",
  "W_ajan_Seward_to_alas_Seward" = "W_ajan_Seward_to_alas_Seward",
  "W_alas_Interior_to_ajan_Interior" = "W_alas_Interior_to_ajan_Interior",
  "W_ajan_Interior_to_alas_Interior" = "W_ajan_Interior_to_alas_Interior"
)

# Convert one chromosome's posterior means into the numeric BPP simulation tree.
make_species_tree <- function(row) {
  sprintf(
    "((ajan_Interior #%g, ajan_Seward #%g)J:%g #%g, (alas_Interior #%g, alas_Seward #%g)L:%g #%g)R:%g #%g;",
    row$theta_ajan_Interior,
    row$theta_ajan_Seward,
    row$tau_J,
    row$theta_J,
    row$theta_alas_Interior,
    row$theta_alas_Seward,
    row$tau_L,
    row$theta_L,
    row$tau_R,
    row$theta_R
  )
}

# List every parameter column required by the control generator so incomplete
# input tables fail early rather than producing malformed BPP controls.
required_parameter_columns <- c(
  "chromosome",
  "theta_ajan_Interior",
  "theta_ajan_Seward",
  "theta_alas_Interior",
  "theta_alas_Seward",
  "theta_R",
  "theta_J",
  "theta_L",
  "tau_R",
  "tau_J",
  "tau_L",
  "W_alas_Seward_to_ajan_Seward",
  "W_ajan_Seward_to_alas_Seward",
  "W_alas_Interior_to_ajan_Interior",
  "W_ajan_Interior_to_alas_Interior"
)
