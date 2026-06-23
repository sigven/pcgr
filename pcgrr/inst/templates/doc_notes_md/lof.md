For SNVs and InDels, PCGR and CPSR consider the following consequence types as
candidates for loss-of-function variation:

* Stop gains - [SO:0001587]({seqontology_url}SO:0001587)
* Frameshift variants - [SO:0001589]({seqontology_url}SO:0001589)
* Splice-site disruptions (2bp donor/acceptor site)
    - [SO:0001574]({seqontology_url}SO:0001574)
    - [SO:0001575]({seqontology_url}SO:0001575)
* Start losses - [SO:0002012]({seqontology_url}SO:0002012)

A collection of filters is next applied, which potentially drops the
loss-of-function annotation for candidates identified above:

* Frameshifts/stop gains within the last 2.5% of the CDS
* Canonical splice site variants that are not predicted to affect a donor site (GC -> GT)
      
If a variant is filtered as non-LoF through any of these criteria, this will
be evident from the `LOF_FILTER` variable (found in the interactive tables
of the HTML report as well as the TSV/Excel output).

Variants at positions immediately flanking the canonical splice sites may
also qualify as loss-of-function candidates when sufficient evidence of
splicing disruption is present. PCGR/CPSR evaluates two classes of
extended splice positions using [MaxEntScan](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html)
delta-loss scores, stratified by position:

* **Donor side** (positions +3 to +6 from the exon/intron boundary):
  strata *Donor_+3*, *Donor_+4*, *Donor_+5*, *Donor_+6* — positions +3/+4/+5
  carry a high prior probability for functional impact; position +6 is medium.
* **Acceptor side** (positions -3 to -20 from the intron/exon boundary):
  stratum *Acceptor_-3* (high prior probability), *Acceptor_PPT_Close*
  (-4 to -8, polypyrimidine tract) and *Acceptor_PPT_Distant* (-9 to -20),
  both with a low prior probability.

For each stratum, position-calibrated thresholds determine an evidence tier
(*supporting*, *moderate*, or *strong*). Only variants reaching the *strong*
tier are marked as loss-of-function candidates; *moderate* and *supporting*
tiers are recorded for reference but do not trigger a loss-of-function call.
The evidence tier is stored in the `MAXENTSCAN` variable.

In addition, variants of other consequence types (e.g. synonymous variants)
found to affect splicing through records in
[MutSpliceDB](https://brb.nci.nih.gov/splicing/) are also marked as
loss-of-function candidates.


