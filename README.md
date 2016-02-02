# tgirtERCC

These are the scripts for generating plots for the paper: [RNA-seq of human reference RNA samples using a thermostable group II intron reverse transcriptase](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3683930/).

In this study, we compare TGIRT-seq with TruSeq v2 and v3 using well-characterized MAQC samples. We find that TGIRT-seq recapitulates the relative abundance of human transcripts and RNA spike-ins in ribo-depleted, fragmented RNA samples comparably to non-strand-specific TruSeq v2 and better than strand-specific TruSeq v3. Moreover, TGIRT-seq is more strand specific than TruSeq v3 and eliminates sampling biases from random hexamer priming, which are inherent to TruSeq. The TGIRT-seq data sets also show more uniform 5' to 3' gene coverage and identify more splice junctions, particularly near the 5' ends of mRNAs, than do the TruSeq data sets. Finally, TGIRT-seq enables the simultaneous profiling of mRNAs and lncRNAs in the same RNA-seq experiment as structured small ncRNAs, including tRNAs, which are essentially absent with TruSeq.

Some **Rscrpts** in this repository depends on the *R* package: [tgirtABRF](https://github.com/wckdouglas/tgirtABRF)
