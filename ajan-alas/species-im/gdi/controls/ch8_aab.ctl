# Species-level ajan vs alas gdi simulation: ch8 aab
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch8_aab.tree.txt
Imapfile = imap/ajan_alas_aab.imap.txt
species&tree = 2 ajan alas
                 2 1
(ajan #0.00811689452, alas #0.01002898252)R:0.0004889371 #0.00364058468;
loci&length = 1000000 50
migration = 2
            ajan alas 848.7356586
            alas ajan 786.9869437
