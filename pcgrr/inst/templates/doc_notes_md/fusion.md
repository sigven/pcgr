RNA fusion analysis in PCGR processes user-provided fusion calls — e.g. output from tools
such as [Arriba](https://github.com/suhrig/arriba) or [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion-Tutorial/wiki) — and
evaluates them for clinical significance through biomarker matching, oncogenicity 
assessment (OncoKB), and druggability annotation. 
The module also incorporates frequency statistics of the gene-level fusion 
events in the [Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer]({mitelman_url}).

Each user-provided fusion call is expected to specify at least one gene partner (5' or 3'), and
the chromosomal breakpoint coordinates. PCGR resolves Ensembl transcript annotations for each
breakpoint, confirming and determining the 5' and 3' gene partners where possible. Fusions where only one
partner is named in the input (e.g. `::ALK` or `FGFR2::`) are supported; the known gene
drives biomarker matching.

Fusion events are matched against clinical biomarker databases (CIViC, CGI, optionally OncoKB) using 
gene-level identifiers. When the matching biomarker is defined without
both partners specified (e.g. _ALK fusions_, _FGFR2 fusions_), only the known gene is
reflected in the `SAMPLE_ALTERATION` column — showing the resolved fusion partner would imply
a partner-specific match that did not occur. Tier classification of gene fusions 
follows the AMP/ASCO/CAP framework [@Li2017-ew] described under _Variant actionability classification_. 

PCGR further annotates whether the 5' or 3' partner genes are targets of cancer drugs, 
drawing on the same druggability resources used for SNV/InDel 
and CNA analysis.
