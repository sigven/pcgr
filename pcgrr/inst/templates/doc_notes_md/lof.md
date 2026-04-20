For SNVs and InDels, PCGR and CPSR consider the following consequence types as
candidates for loss-of-function variation:

* Stop gains - [SO:0001587]({seqontology_url}SO:0001587)
* Frameshift variants - [SO:0001589]({seqontology_url}SO:0001589)
* Splice-site disruptions (2bp donor/acceptor site)
    - [SO:0001574]({seqontology_url}SO:0001574)
    - [SO:0001575]({seqontology_url}SO:0001575)
* Splice donor 5th base disruptions - [SO:0001787]({seqontology_url}SO:0001787)
* Start losses - [SO:0002012]({seqontology_url}SO:0002012)

If variants of other consequence types (e.g. synonymous variants, or splice
site variants beyond the canonical 2bp site/5th donor base) are found to
affect splicing, specifically through records in
[MutSpliceDB](https://brb.nci.nih.gov/splicing/), these are also marked as
loss-of-function candidates.

A collection of filters is next applied, which potentially drops the
loss-of-function annotation for candidates identified above:

* Frameshifts/stop gains within the last 5% of the CDS
* Splice site variants that are not predicted to affect a donor site (GC -> GT)
* Variants where [MaxEntScan](https://pubmed.ncbi.nlm.nih.gov/15285897/)
  does not predict an effect on splicing
    - A drop in MaxEntScan score (between reference and alternative allele)
      of less than 60%. This is only considered where the minimum MES score
      for the reference site is at least 4.0 and 5.0 for donor and
      acceptor sites, respectively.
    - Annotations in MutSpliceDB will have precedence if any conflicting
      evidence with MaxEntScan output is found
      
If a variant is filtered as non-LoF through any of these criteria, this will
be evident from the `LOF_FILTER` variable (found in the interactive tables
of the HTML report as well as the TSV/Excel output).
