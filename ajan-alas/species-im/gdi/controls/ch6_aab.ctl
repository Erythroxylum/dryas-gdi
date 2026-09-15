# Species-level ajan vs alas gdi simulation: ch6 aab
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch6_aab.tree.txt
Imapfile = imap/ajan_alas_aab.imap.txt
species&tree = 2 ajan alas
                 2 1
(ajan #0.00718033758, alas #0.0108134556)R:0.00051674482 #0.003373695;
loci&length = 1000000 50
migration = 2
            ajan alas 1373.909377
            alas ajan 1447.889097
