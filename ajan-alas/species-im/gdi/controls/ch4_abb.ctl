# Species-level ajan vs alas gdi simulation: ch4 abb
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch4_abb.tree.txt
Imapfile = imap/ajan_alas_abb.imap.txt
species&tree = 2 ajan alas
                 1 2
(ajan #0.00685977052, alas #0.0101965533)R:0.00056236228 #0.00352509714;
loci&length = 1000000 50
migration = 2
            ajan alas 417.138773
            alas ajan 1192.76206
