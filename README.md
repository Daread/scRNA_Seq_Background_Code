# scRNA_Seq_Background_Code
Methods for background correction in single cell RNA-Seq data

Of the background methods we use, the SoupX method and Packer
method aren't built into Monocle readily.

**SoupX**
See https://doi.org/10.1093/gigascience/giaa151, by Young and Behjati

For SoupX, you can do one version that subtracts signal as if all
samples are exposed to a uniform supernatant (if most leakage is
during the library prep itself). That is under
"correctCDSbySoupX_all_at_once". Another version does correction
separately based on cells/background from each sample (if most
background is coming from initial dissociation and leakage), under
"correctCDSbySoupX_individually". I might lean towards the "individually" version, especially if you're correcting a CDS with a heterogenous mix of tissues. When I ran this on a dataset that was all very similar kidneys, it seemed to not really matter, but there wasn't a ton of inter-sample diversity.

Some other parameters need to be given (like SoupX needs an
estimate for background fraction, where I often use something like 40
as an estimate of 40% background). Note: The SoupX publication itself
recommends empirically determining the background percentage by doing something 
like counting globin mRNA captures by nuclei in snRNA-seq data, and applying that
going forward. I have not implemented a function to encompass that approach,
though it could be useful in the future.


**The "Packer" or "Jonathan" Method**
See PMID: 31488706, Packer et al. 2019

For the Packer correction method, run findBackground_Packer_Method.

**General**
For both, the expected input is a Monocle CDS object for cells, and a
second CDS for background ("cells" with low UMI, typically <15 in my
work).  Also, make sure that there is a
column named "sample" in the background cds "coldata" dataframe when
running the Packer method or the by-sample SoupX method or it'll
crash. 
