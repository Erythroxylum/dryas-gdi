# Species-level ajan vs alas gdi simulation: ch7 aab
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch7_aab.tree.txt
Imapfile = imap/ajan_alas_aab.imap.txt
species&tree = 2 ajan alas
                 2 1
(ajan #0.00746452176, alas #0.0073036575)R:0.00062057358 #0.0034379244;
loci&length = 1000000 50
migration = 2
            ajan alas 596.9625713
            alas ajan 217.5499949
