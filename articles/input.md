# Input files

## Input files

The PCGR workflow accepts the following main input files:

- An unannotated, single-sample [VCF
  file](https://github.com/samtools/hts-specs#variant-calling-data-files)
  (\>= v4.2) with called somatic variants (SNVs/InDels)
- A file with allele-specific copy number segments (tab-separated
  values - TSV)
- A file with RNA fusion transcripts detected in the tumor
  (tab-separated values - TSV)
- A file with transcript/gene expression levels (tab-separated values -
  TSV)

At least one of these input files is required. PCGR can be run with any
combination of the above, and will only produce output sections relevant
to the data provided. The following arguments to the `pcgr` command are
used for input files:

- `--input_vcf`
- `--input_cna`
- `--input_rna_fusion`
- `--input_rna_expression`

In addition to these main input files, the user can also opt to provide
a [panel-of-normals VCF](#panel-of-normals-pon-vcf) file (`--pon_vcf`)
for tumor-only variant filtering, as well as an input file with
[CPSR-classified germline variants (TSV)](#germline-variants)
(`--input_cpsr`).

### VCF

- We **strongly** recommend that the input VCF is compressed and indexed
  using [bgzip](http://www.htslib.org/doc/bgzip.md) and
  [tabix](http://www.htslib.org/doc/tabix.md)
- If the input VCF contains [multi-allelic
  sites](https://glow.readthedocs.io/en/latest/etl/variant-splitter.html),
  these will be subject to
  [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose). Either
  way, we encourage that users prepare the input VCF *without the
  presence of multi-allelic sites*.
- Variants used for reporting should be designated as `PASS` in the VCF
  FILTER column. Variants denoted with e.g. `Reject` as a FILTER value
  will not be subject to analysis in PCGR. For records with undefined
  values in the FILTER column (`'.'`), these will be considered as
  `PASS` variants.

#### Formatting of variant sequencing depth/allelic support (DP/AF)

The representation of variant genotype data (allelic depth and support
in tumor vs. control sample) is usually formatted in the genotype fields
of a VCF file, on a per-sample basis. However, considering the VCF
output for the [numerous somatic SNV/InDel
callers](https://www.biostars.org/p/19104/) that are in use, we have
experienced a general lack of uniformity for how this information is
encoded in the genotype fields. In order for PCGR to recognize this type
of information robustly, we currently require that you as a user encode
this unambiguously in the INFO field of your VCF file.

Shown below is how the VCF header for these entries should look like in
your input VCF, and how the corresponding variant data is encoded per
record:

**VCF header**

    ##INFO=<ID=TVAF,Number=.,Type=Float,Description="Allelic fraction of alternative allele in tumor">
    ##INFO=<ID=TDP,Number=.,Type=Integer,Description="Read depth across variant site in tumor">

**VCF records**

    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO
    1       23519660        .       G       C       .       PASS    TVAF=0.0980;TDP=51

When you run PCGR, make sure to feed the name of the INFO tags to the
tool, i.e. using the `tumor_dp_tag`and `tumor_af_tag` options:

    pcgr --tumor_dp_tag TDP --tumor_af_tag TVAF

As an effort to support the users with this procedure, we are hoping to
establish a *VCF pre-processing script*, that accepts any VCF (from any
caller), and that will ensure that the file adheres to the formatting
requirements of PCGR. See also [this
issue](https://github.com/sigven/pcgr/issues/136).

#### Other notes regarding input VCF

- PCGR generates a number of VCF INFO annotation tags that are appended
  to the query VCF. We therefore encourage the users to submit query VCF
  files that *have not* been subject to annotations by other means, but
  rather a VCF file that comes directly from variant calling. If not,
  there are likely to be INFO tags in the query VCF file that coincide
  with those produced by PCGR.
- Note that you can preserve particular tags in the INFO field towards
  the TSV output of PCGR. Sometimes, it can be convenient to investigate
  particular properties of the variants (encoded in the VCF) against
  functional annotations (as provided by PCGR). To achieve this, use the
  option `--retained_info_tags <TAG1>,<TAG2>` etc.

### Panel-of-normals VCF

The user can submit a panel-of-normals (PoN) VCF file to PCGR,
indicating variants that should be excluded/filtered from the input VCF
in a **tumor-only setting** (i.e. option `--tumor_only`).

- `--pon_vcf <VCF_FILE>`

The PoN VCF file needs to contain the following.

1.  The following description needs to be present in the VCF header:

&nbsp;

    ##INFO=<ID=PANEL_OF_NORMALS,Number=0,Type=Flag,Description="Overlap with germline call among panel of normals">

2.  Each record in the PoN VCF file needs to contain the
    `PANEL_OF_NORMALS` flag in the INFO field

### Copy number segments - allele-specific

A tab-separated values file with allele-specific copy number aberrations
**MUST** contain the following four columns:

- `Chromosome`
- `Start`
- `End`
- `nMajor`
- `nMinor`

Here, *Chromosome*, *Start*, and *End* denote the chromosomal segment,
while *nMajor* and *nMinor* denote major and minor copy number for
particular segment, a standard output of somatic copy number alteration
callers. Note that coordinates must be **one-based** (i.e. chromosomes
start at 1, not 0), and that copy numbers should be formatted as whole
integers, not floating values. Below shows the initial part of a copy
number segment file that is formatted correctly according to PCGR’s
requirements:

      Chromosome    Start   End nMajor nMinor
      1 3218329 3550598 1 0
      1 3552451 4593614 1 1
      1 4593663 6433129 2 1

Importantly, you can configure predefined thresholds for segments that
are considered to be high-level copy-number amplifications through the
following option:

- `--n_copy_gain`: minimum (total) copy number for a segment to be
  considered an amplification, *default 6*

### Germline variants

The user can submit a file with germline variants processed and
classified with [CPSR](https://sigven.github.io/cpsr), which caters for
an integration of germline and somatic findings in the output report.
This file corresponds to output file
`<sample_id>.cpsr.<genome_assembly>.classification.tsv.gz` from the CPSR
pipeline. Make sure the genome assembly is the same as the one used for
the somatic variant input files.

### RNA fusions

The user can submit a file with RNA fusion transcripts detected in the
tumor sample. PCGR will cross-reference detected fusions against curated
cancer gene databases and known fusion events, and incorporate them into
the output report.

The tab-separated values file with RNA fusion calls **MUST** contain the
following four columns:

| Column | Type | Description |
|----|----|----|
| `FusionGene` | string | Fusion gene pair, two gene symbols separated by `--` or `::` (e.g. `EML4--ALK` or `TMPRSS2::ERG`) |
| `LeftBreakpoint` | string | Genomic position of the left/5’ breakpoint, formatted as `chrom:position` (e.g. `2:42522694`) |
| `RightBreakpoint` | string | Genomic position of the right/3’ breakpoint, formatted as `chrom:position` (e.g. `2:29446394`) |
| `SplitReads` | integer | Number of split reads supporting the fusion (non-negative integer) |

An optional column `Score` (float) may also be included, representing a
caller-specific confidence score for the fusion event.

#### Formatting notes

- The chromosome field in breakpoints must be one of `1`–`22`, `X`, or
  `Y` — with or without a `chr` prefix
- Each gene symbol in `FusionGene` must be a valid alphanumeric
  identifier (letters, digits, `.`, `-`, `_` allowed); purely numeric
  gene names are not accepted
- Fusions with an empty or missing partner (e.g. intragenic or
  read-through events) are allowed but will generate a warning during
  validation

Below shows the first few records of a correctly formatted RNA fusion
file:

    FusionGene       LeftBreakpoint    RightBreakpoint   SplitReads
    EML4--ALK        2:42522694        2:29446394        100
    BCR--ABL1        22:23582899       9:13358994        50
    EWSR1--FLI1      22:29683136       11:128675244      75
    TMPRSS2--ERG     21:42836478       21:39737961       60

**Quality control of RNA fusion events**: We strongly advice users to to
perform quality control of fusion events upstream of PCGR analysis - the
focus in PCGR is strictly on functional annotation and potential
clinical impact.

- We recommend to carefully check e.g. the number of supporting reads
  for each fusion event, and to visually inspect fusion events in a
  genome browser (e.g. IGV) or other third-party tools
  (e.g. [FusViz](https://github.com/senzhaocode/FuSViz)) to conduct
  frame inference, investigate domain retention, and rule out potential
  false positives.
- Importantly, while most RNA fusion detection tools will attempt to
  filter out known normal/artefactual events, PCGR makes no attempt to
  further filter fusion events based on e.g. presence in normal tissue
  databases, and thus relies on the user to perform this step upstream
  of PCGR analysis.
- If multiple input entries share the same gene pair
  (e.g. *GENE_A–GENE_B*) but differ only in their exact breakpoint
  coordinates, each will appear as a separate entry in the PCGR report.
  This can make the fusion section appear repetitive. PCGR does not
  attempt to collapse these, as there is no principled way to select a
  canonical breakpoint without the full caller context. We therefore
  recommend that users deduplicate fusion calls to one entry per gene
  pair — retaining the call with the highest read support or confidence
  score — before submitting to PCGR.

### Gene expression

The user can submit a file with bulk gene/transcript expression data to
PCGR, indicating the relative expression levels of genes in the query
sample. PCGR may conduct an expression outlier analysis (compared
against other cohorts of samples), and also perform an RNA expression
similarity analysis, i.e. correlation of expression profile with other
tumor samples.

The tab-separated values file with gene expression estimates **MUST**
contain the following two columns:

- `TargetID`
- `TPM`

The `TPM` column indicates the expression level of a given target, and
must be provided as *Transcripts per million* (not log-transformed), a
relatively common output measure from bulk RNA-seq pipelines. The
`TargetID` should list transcript or gene identifiers, preferably at the
transcript level. Permitted identifier types including Ensembl/RefSeq
identifiers for transcript-level expression estimates, and gene symbols
and Ensembl gene identifiers for gene-level expression estimates.
