# Species-level ajan vs alas gdi simulation: ch3 abb
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch3_abb.tree.txt
Imapfile = imap/ajan_alas_abb.imap.txt
species&tree = 2 ajan alas
                 1 2
(ajan #0.0070065882, alas #0.00806023984)R:0.00063712126 #0.00326632736;
loci&length = 1000000 50
migration = 2
            ajan alas 690.4813885
            alas ajan 658.370632
