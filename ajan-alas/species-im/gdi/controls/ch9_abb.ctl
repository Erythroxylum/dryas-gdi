# Species-level ajan vs alas gdi simulation: ch9 abb
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch9_abb.tree.txt
Imapfile = imap/ajan_alas_abb.imap.txt
species&tree = 2 ajan alas
                 1 2
(ajan #0.00569416528, alas #0.00845283732)R:0.00060249078 #0.00311534774;
loci&length = 1000000 50
migration = 2
            ajan alas 539.7557026
            alas ajan 582.9772252
