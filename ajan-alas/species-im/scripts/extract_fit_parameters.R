# Extract posterior means from the nine completed species-level IM fits.
# Each MCMC file contains chromosome-specific estimates of theta, tau,
# and directional migration rates for ajanensis and alaskensis.

library(readr)
library(dplyr)
library(purrr)

chromosomes <- paste0("ch", 1:9)

params <- map_dfr(chromosomes, function(chrom) {

  f <- file.path(
    "ajan-alas/species-im/fit/output",
    paste0(chrom, ".mcmc.txt")
  )

  x <- read_tsv(f, show_col_types = FALSE)

  tibble(
    chromosome = chrom,
    theta_ajan = mean(x$`theta:1:ajan`),
    theta_alas = mean(x$`theta:2:alas`),
    theta_R    = mean(x$`theta:3:R`),
    tau_R      = mean(x$`tau:3:R`),
    W_ajan_to_alas = mean(x$`W:1->2:ajan->alas`),
    W_alas_to_ajan = mean(x$`W:2->1:alas->ajan`)
  )
})

write_csv(
  params,
  "ajan-alas/species-im/fit/species_im_means.csv"
)

print(params)
