# Species-level ajan vs alas gdi simulation: ch1 aab
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch1_aab.tree.txt
Imapfile = imap/ajan_alas_aab.imap.txt
species&tree = 2 ajan alas
                 2 1
(ajan #0.00775191808, alas #0.00909450528)R:0.00050758784 #0.00312194378;
loci&length = 1000000 50
migration = 2
            ajan alas 354.2521204
            alas ajan 704.0553088
