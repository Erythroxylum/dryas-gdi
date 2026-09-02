# Species-level ajan vs alas IM fit for ch8
# Empirical sequence file: /Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/0Dryas/panDryas/BPP-yuttapong/bpp_may12/data/ajan-alas-s20/multilocus/ch8.txt
# Source sample map: /Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/0Dryas/panDryas/BPP-yuttapong/bpp_may12/imap/ajan-alas-s20-p4.imap.txt
# Source BPP template: /Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/0Dryas/panDryas/BPP-yuttapong/bpp_may12/ctl/bpp-a00-ajan-alas-s20-p4.ctl
# Interior and Seward are collapsed within species; migration is fitted in both directions.
# Phase = 1 1 because the collapsed model has two unphased species.
# Migration prior: W ~ Gamma(2, 0.01), matching m3-prior2.
          seed = -1

seqfile = /Users/dawsonwhite/Library/Mobile Documents/com~apple~CloudDocs/0Dryas/panDryas/BPP-yuttapong/bpp_may12/data/ajan-alas-s20/multilocus/ch8.txt
Imapfile = /Users/dawsonwhite/Desktop/dryas-gdi/ajan-alas/species-im/fit/imap/ajan_alas.imap.txt
jobname = fit/output/ch8

  speciesdelimitation = 0                * fixed species tree
          speciestree = 0                * speciestree fixed

species&tree = 2 ajan alas
                 9 11
(ajan, alas)R;
phase = 1 1

       usedata = 1   * 0: no data(prior); 1:seq Like
         nloci = 2000  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 400  # gamma(a, b) for theta
      tauprior = gamma 2 400  # gamma(a, b) for root tau & Dirichlet(a) for other tau's
wprior = 2 0.01
migration = 2
            ajan alas
            alas ajan

      finetune = 1

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 10000
      sampfreq = 20
       nsample = 50000
