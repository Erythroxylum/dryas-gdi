          seed = -1

       seqfile = loci.txt
      Imapfile = Imap.txt
       jobname = out

  speciesdelimitation = 0                * fixed species tree
          speciestree = 0                * speciestree fixed

 species&tree = 4  ajan_Interior ajan_Seward alas_Interior alas_Seward
                               5           4             6           5
                   ((ajan_Interior, alas_Interior)I, (ajan_Seward, alas_Seward)S)R;
         phase = 1 1 1 1  * 0: phased sequences, 1: unphased diploid sequences

       usedata = 1   * 0: no data(prior); 1:seq Like
         nloci = 2000  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 400  # gamma(a, b) for theta
      tauprior = gamma 2 400  # gamma(a, b) for root tau & Dirichlet(a) for other tau's
        wprior = 2 1  # gamma(a, b)
     migration = 14
                 ajan_Interior alas_Interior
                 ajan_Interior ajan_Seward
                 ajan_Interior alas_Seward
                 alas_Interior ajan_Interior
                 alas_Interior ajan_Seward
                 alas_Interior alas_Seward
                 ajan_Seward ajan_Interior
                 ajan_Seward alas_Interior
                 ajan_Seward alas_Seward
                 alas_Seward ajan_Interior
                 alas_Seward alas_Interior
                 alas_Seward ajan_Seward
                 I S
                 S I

      finetune = 1

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 10000
      sampfreq = 20
       nsample = 50000
