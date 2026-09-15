# Species-level ajan vs alas gdi simulation: ch5 abb
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch5_abb.tree.txt
Imapfile = imap/ajan_alas_abb.imap.txt
species&tree = 2 ajan alas
                 1 2
(ajan #0.00637325904, alas #0.00811242926)R:0.0005729294 #0.00329948964;
loci&length = 1000000 50
migration = 2
            ajan alas 728.15609
            alas ajan 771.8630417
