# Species-level ajan vs alas gdi simulation: ch2 aab
# Parameters are chromosome-specific posterior means from the empirical IM fit.
# aab estimates gdi_ajan; abb estimates gdi_alas.
seed = -1
treefile = output/trees/ch2_aab.tree.txt
Imapfile = imap/ajan_alas_aab.imap.txt
species&tree = 2 ajan alas
                 2 1
(ajan #0.00514907714, alas #0.00690321738)R:0.00075964238 #0.00324898302;
loci&length = 1000000 50
migration = 2
            ajan alas 394.2850715
            alas ajan 653.2638739
