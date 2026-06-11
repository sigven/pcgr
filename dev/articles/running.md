# Running

The PCGR software comes with a range of options to configure the
analyses performed and to customize the output report.

## Key settings

Below, we will outline key settings that are important with respect to
the sequencing assay that has been conducted.

### Whole-genome, whole-exome or targeted sequencing assay

By default, PCGR expects that the input VCF comes from whole-exome
sequencing. This can be adjusted with the `assay` option, which takes
three alternative values:

- `--assay <WGS|WES|TARGETED>`

If the input VCF comes from a targeted sequencing assay/panel
(i.e. `--assay TARGETED`), the target size of the assay/panel should be
adjusted accordingly, which is critical for tumor mutational burden
estimation. Assay target size, reflecting the *protein-coding portion of
the assayed genomic region* in megabases, can be configured through the
following option:

- `--effective_target_size_mb <value>`

Ideally, this should reflect the *callable* target size of the assay,
i.e. the size of the target regions for which coding variants is
sufficiently covered (i.e. with a certain sequencing depth) by the
sequencing assay. This can be estimated by the user based on the design
of the sequencing assay, or by using tools such as
[mosdepth](https://github.com/brentp/mosdepth).

### Tumor-control vs. tumor-only

By default, PCGR expects that the input VCF contains somatic variants
identified from a tumor-control sequencing setup. This implies that the
VCF contains information with respect to variant allelic depth/support
both for the tumor sample and the corresponding control sample.

If the input VCF comes from a **tumor-only** assay, turn on the
`--tumor_only` option. In this mode, PCGR conducts a set of successive
filtering steps on the raw input set of variants, aiming to exclude the
majority of germline variants from the tumor-only input set. In addition
to default filtering applied against variants found in the gnomAD
database (population-specific minor allele frequency thresholds can be
configured, see below), additional filtering procedures can be
explicitly set, i.e.:

- `--exclude_dbsnp_nonsomatic`
- `--exclude_likely_het_germline`
- `--exclude_likely_hom_germline`
- `--exclude_clinvar_germline`
- `--exclude_nonexonic`

Users should note that these filters may occasionally remove true
somatic variants, and should thus be used with caution.

If the user has provided a [panel-of-normals VCF
file](https://sigven.github.io/pcgr/dev/articles/input.html#panel-of-normals-vcf)
as input(`--pon_vcf`), one may also opt to remove variants found in that
variant collection using a dedicated option:

- `--exclude_pon`

------------------------------------------------------------------------

**IMPORTANT NOTE 1:** A number of the analyses available in PCGR are not
particularly well-suited for tumor-only input, e.g. mutational signature
analysis and MSI status prediction.

**IMPORTANT NOTE 2:** If you run PCGR on tumor-only WGS assays, we
strongly recommend that you pre-filter your VCF against the exome, since
we have frequently encountered that the large unfiltered variant set
from such assays will cause a crash in PCGR. See also [this
issue](https://github.com/sigven/pcgr/issues/178).

### Sample properties

A few selected properties of the tumor sample that should be determined
prior to PCGR analysis can be provided

- `--tumor_ploidy <value>` - ploidy estimate
- `--tumor_purity <value>` - purity estimate

If not provided, PCGR will attempt to estimate ploidy from the copy
number segments (if provided) using a weighted median approach. If copy
number segments are not provided, these properties will be set to `NA`
in the report.

#### Tumor site

The primary tissue/site of the tumor sample that is analyzed can and
should be be provided, through the following option:

- `--tumor_site <value>`

This option takes a value between 0 and 30, reflecting the main primary
sites/tissues for human cancers. This information is used e.g. to
discriminate between on and off-label targeted therapies from the
actionable biomarkers detected in the sample. Note that the default
value is **0** (implies **Any** tissue/site) - only pan-cancer
biomarkers will be considered for Tier I assignment, and all
tumor-specific biomarkers will be considered as non-matching site
evidence for Tier II/III assignment.

The following lists the possible encodings, and the corresponding
primary tissues/sites:

    0 = Any
    1 = Adrenal Gland
    2 = Ampulla of Vater
    3 = Biliary Tract
    4 = Bladder/Urinary Tract
    5 = Bone
    6 = Breast
    7 = Cervix
    8 = CNS/Brain
    9 = Colon/Rectum
    10 = Esophagus/Stomach
    11 = Eye
    12 = Head and Neck
    13 = Kidney
    14 = Liver
    15 = Lung
    16 = Lymphoid
    17 = Myeloid
    18 = Ovary/Fallopian Tube
    19 = Pancreas
    20 = Peripheral Nervous System
    21 = Peritoneum
    22 = Pleura
    23 = Prostate
    24 = Skin
    25 = Soft Tissue
    26 = Testis
    27 = Thymus
    28 = Thyroid
    29 = Uterus
    30 = Vulva/Vagina
    (default: 0 - any tumor type)

If you want to check the nature of cancer subtypes attributed to each of
the primary tissue/sites, take a look at the [table with phenotype terms
and associated primary
sites](https://sigven.github.io/pcgr/dev/articles/primary_tumor_sites.md).

### Allelic support

Proper designation of INFO tags in the input VCF that denote variant
depth/allelic fraction is critical to make the report as comprehensive
as possible. See the [input
section](https://sigven.github.io/pcgr/dev/articles/input.md) on how to
accomplish this.

Several options are used to inform PCGR about tags in the VCF INFO field
that encodes this information:

- `--tumor_dp_tag <value>`
- `--tumor_af_tag <value>`
- `--control_dp_tag <value>`
- `--control_af_tag <value>`

If these tags are set correctly, one may set thresholds (sequencing
depth and/or allelic depth/fraction) on the variants to be included in
the report, i.e. through:

- `--tumor_dp_min <value>`
- `--tumor_af_min <value>`
- `--tumor_ad_min <value>`
- `--control_dp_min <value>`
- `--control_af_max <value>`
- `--control_ad_max <value>`

The *allelic depths tresholds* refer to the minimum number of reads
supporting the variant allele in the tumor sample (`--tumor_ad_min`),
and the maximum number of reads supporting the variant allele in the
control sample (`--control_ad_max`).

### Tumor mutational burden (TMB)

If tags for allelic support is provided in the VCF and configured by the
user, users can configure the TMB calculation by setting minimum
requirements for sequencing coverage and allelic fraction, i.e.:

- `--tmb_dp_min <value>`
- `--tmb_af_min <value>`

### Large input sets (VCF)

For input VCF files with \> 500,000 variants, note that these will be
subject to filtering prior to the final step in PCGR (step 4:
reporting). For such scenarios, raw input sets will be filtered
according to consequence type, in the sense that
intergenic/intronic/upstream/downstream gene variants will be excluded
prior to reporting. This is a necessary step to avoid memory issues etc.
during processing with the PCGR R package.

If you have a large input VCF, and have sufficient memory capacity on
your compute platform, we also recommend to increase the [VEP buffer
size](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#cacheopt)
(option `--vep_buffer_size`), as this will speed up the VEP processing
significantly.

### Somatic copy number variant class thresholds

When copy number segments are provided (`--input_cna`), PCGR classifies
each segment as an amplification, gain, heterozygous deletion, or
homozygous deletion based on configurable thresholds. The thresholding
strategy is controlled by:

- `--cna_threshold_mode <absolute|relative|combined>` (default:
  `absolute`)

Three modes are supported:

- **absolute** — thresholds are expressed as total copy numbers (integer
  values). A segment is classified based on its total copy number
  (`nMajor + nMinor`), irrespective of tumor ploidy.
- **relative** — thresholds are expressed as fold-change over the
  estimated tumor ploidy. If `--tumor_ploidy` is not provided explicitly
  by the user, PCGR will auto-estimate the ploidy from the copy number
  segments using a genome-wide weighted median approach.
- **combined** — a segment must satisfy *both* the absolute and relative
  threshold criteria to be assigned a given class. This is the most
  conservative mode.

The per-class thresholds that can be configured are:

| Class | Absolute option | Default | Relative option | Default |
|----|----|----|----|----|
| Amplification | `--cna_amp_threshold_absolute` | Total CN ≥ 5 | `--cna_amp_threshold_relative` | ≥ 2.5× ploidy |
| Gain | `--cna_gain_threshold_absolute` | Total CN ≥ 3 | `--cna_gain_threshold_relative` | ≥ 1.5× ploidy |
| Heterozygous deletion | `--cna_del_threshold_absolute` | Total CN ≤ 1 | `--cna_del_threshold_relative` | ≤ 0.5× ploidy |

Segments with total copy number of zero are always classified as
homozygous deletions, regardless of the threshold mode.

An additional option controls which genes are considered affected by a
CNA event:

- `--cna_transcript_overlap_pct <value>` — minimum mean percent overlap
  between the copy number segment and gene transcripts required to
  report a gene as affected (default: 50)

### OncoKB biomarker integration

PCGR can optionally query the [OncoKB](https://www.oncokb.org) precision
oncology knowledge base to annotate SNVs/InDels, copy number alterations
and fusions with functional annotations and clinical actionability
evidence. To enable this, provide a valid OncoKB API token:

- `--oncokb_api_token <ONCOKB_API_TOKEN>`

To specify the tumor type for OncoKB queries, provide the relevant
[OncoTree](https://oncotree.mskcc.org) code:

- `--oncokb_oncotree_code <CODE>` (e.g. `COAD` for colon adenocarcinoma)

If no OncoTree code is provided, PCGR will map the primary tumor site to
its corresponding OncoTree code based on a predefined mapping. If the
tumor site is set to “Any” (code 0), OncoKB queries will be performed
without specifying a tumor type.

By default, PCGR integrates biomarker evidence from multiple knowledge
sources (CIViC, CGI, OncoKB). If you want to limit biomarker reporting
to OncoKB only and skip CIViC and CGI, use:

- `--oncokb_exclusive`

**Licensing notice — OncoKB terms of use:** OncoKB is freely accessible
for non-commercial research use. However, use of OncoKB data for
**commercial purposes or for patient care in a hospital or clinical
setting** — including the generation of PCGR reports used to inform
diagnostic or therapeutic decisions — requires a separate commercial
license from Memorial Sloan Kettering Cancer Center. Users are solely
responsible for ensuring that their use of OncoKB through PCGR complies
with the [OncoKB Terms of Use](https://www.oncokb.org/terms). If in
doubt, contact the OncoKB team at <contact@oncokb.org> before use.

### RNA fusions

If RNA fusion data is available for the tumor sample, it can be provided
via:

- `--input_rna_fusion <FILE>`

PCGR will cross-reference the detected fusions against curated cancer
gene databases and known fusion events. A minimum split-read support
threshold can be configured to filter low-confidence fusions:

- `--fusion_min_split_reads <value>` (default: 3, minimum: 2)

See the [input
section](https://sigven.github.io/pcgr/dev/articles/input.html#rna-fusions)
for details on the required file format.

### Configuration of output files

If you only want the Excel/TSV output (with all variant classifications
and auxiliary analyses (e.g. MSI classifications, TMB estimates,
mutational signatures predictions) of PCGR, you may turn off the HTML
output by using the `--no_html` option. This will speed up the analysis,
and could be a useful option if you foresee that the sample input
datasets are too large for the HTML report generation to work properly.

If you only want the annotated VCF (also converted to TSV), use the
option `--no_reporting`. This will skip the final steps of the PCGR
workflow, and only generate the annotated VCF and TSV files.

## All options

A tumor sample report is generated by running the **pcgr** command,
which takes the following arguments and options:

``` text
usage:
  pcgr -h [options]
  [--input_vcf <INPUT_VCF>]
  [--input_cna <INPUT_CNA>]
  [--input_rna_fusion <INPUT_RNA_FUSION>]
  [--input_rna_expression <INPUT_RNA_EXPRESSION>]
  [--vep_dir <VEP_DIR>]
  --refdata_dir <REFDATA_DIR>
  --output_dir <OUTPUT_DIR>
  --genome_assembly <GENOME_ASSEMBLY>
  --sample_id <SAMPLE_ID>

Personal Cancer Genome Reporter (PCGR) workflow for clinical translation of tumor omics data
(SNVs/InDels, CNA, RNA expression, RNA fusions)

Required arguments:
  --refdata_dir REFDATA_DIR
                        Directory where PCGR reference data bundle was downloaded and unpacked
  --output_dir OUTPUT_DIR
                        Output directory
  --genome_assembly {grch37,grch38}
                        Human genome assembly build: grch37 or grch38
  --sample_id SAMPLE_ID
                        Tumor sample/cancer genome identifier - prefix for output files

Input file options:
  --input_vcf INPUT_VCF
                        VCF input file with somatic variants in tumor sample, SNVs/InDels
  --input_cna INPUT_CNA
                        Somatic copy number alteration segments (tab-separated values)
  --input_rna_fusion INPUT_RNA_FUSION
                        File with RNA fusion transcripts detected in tumor (tab-separated values)
  --input_rna_expression INPUT_RNA_EXP
                        File with bulk RNA expression counts (TPM) of transcripts in tumor (tab-separated values)

SNV/InDel analysis options:
  --vep_dir VEP_DIR     Directory of VEP cache, e.g. $HOME/.vep (required when --input_vcf is provided)
  --assay {WGS,WES,TARGETED}
                        Type of DNA sequencing assay performed for input data (VCF), default: WES
  --effective_target_size_mb EFFECTIVE_TARGET_SIZE_MB
                        Effective target size in Mb (potentially limited by read depth) of sequencing
                        assay (for TMB analysis) (default: 34 (WES/WGS))
  --tumor_only          Input VCF comes from tumor-only sequencing, calls will be filtered for variants
                        of germline origin (default: False)
  --vcf2maf             Generate a MAF file for input VCF using https://github.com/mskcc/vcf2maf
                        (default: False)
  --vcfanno_n_proc VCFANNO_N_PROC
                        Number of vcfanno processes (option '-p' in vcfanno), default: 4
  --retained_info_tags RETAINED_INFO_TAGS
                        Comma-separated string of VCF INFO tags from query VCF to keep in PCGR output TSV
  --ignore_noncoding    Ignore non-coding (i.e. non protein-altering) variants in report, default: False
  --vep_n_forks VEP_N_FORKS
                        Number of forks (VEP option '--fork'), default: 4
  --vep_buffer_size VEP_BUFFER_SIZE
                        Variant buffer size (variants read into memory simultaneously, VEP option
                        '--buffer_size') - set lower to reduce memory usage, default: 500
  --vep_pick_order VEP_PICK_ORDER
                        Comma-separated string of ordered transcript/variant properties for selection of
                        primary variant consequence (option '--pick_order' in VEP), default:
                        mane_select,mane_plus_clinical,canonical,biotype,ccds,rank,tsl,appris,length
  --vep_no_intergenic   Skip intergenic variants during variant annotation (VEP option '--no_intergenic'),
                        default: False
  --vep_regulatory      Add VEP regulatory annotations (VEP option '--regulatory'), default: False
  --vep_gencode_basic   Consider basic GENCODE transcript set only (VEP option '--gencode_basic')

Tumor sample options:
  --sex {FEMALE,MALE,UNKNOWN}
                        Sex of cancer case/sample (default: UNKNOWN)
  --tumor_site TSITE    Optional integer code to specify primary tumor type/site of query sample,
                        choose any of the following identifiers:
                        0 = Any
                        1 = Adrenal Gland
                        2 = Ampulla of Vater
                        3 = Biliary Tract
                        4 = Bladder/Urinary Tract
                        5 = Bone
                        6 = Breast
                        7 = Cervix
                        8 = CNS/Brain
                        9 = Colon/Rectum
                        10 = Esophagus/Stomach
                        11 = Eye
                        12 = Head and Neck
                        13 = Kidney
                        14 = Liver
                        15 = Lung
                        16 = Lymphoid
                        17 = Myeloid
                        18 = Ovary/Fallopian Tube
                        19 = Pancreas
                        20 = Peripheral Nervous System
                        21 = Peritoneum
                        22 = Pleura
                        23 = Prostate
                        24 = Skin
                        25 = Soft Tissue
                        26 = Testis
                        27 = Thymus
                        28 = Thyroid
                        29 = Uterus
                        30 = Vulva/Vagina
                        (default: 0 - any tumor type)
  --tumor_purity TUMOR_PURITY
                        Estimated tumor purity (between 0 and 1) (default: None)
  --tumor_ploidy TUMOR_PLOIDY
                        Estimated tumor ploidy (default: None)

Allelic support options:
  --tumor_dp_tag TUMOR_DP_TAG
                        Specify VCF INFO tag for sequencing depth (tumor, must be Type=Integer,
                        default: _NA_)
  --tumor_af_tag TUMOR_AF_TAG
                        Specify VCF INFO tag for variant allelic fraction (tumor, must be Type=Float,
                        default: _NA_)
  --control_dp_tag CONTROL_DP_TAG
                        Specify VCF INFO tag for sequencing depth (control, must be Type=Integer,
                        default: _NA_)
  --control_af_tag CONTROL_AF_TAG
                        Specify VCF INFO tag for variant allelic fraction (control, must be Type=Float,
                        default: _NA_)
  --call_conf_tag CALL_CONF_TAG
                        Specify VCF INFO tag for somatic variant call confidence (must be categorical,
                        e.g. Type=String, default: _NA_)
  --tumor_dp_min TUMOR_DP_MIN
                        Minimum sequencing depth (tumor) for inclusion in report (default: None)
  --tumor_ad_min TUMOR_AD_MIN
                        Minimum allelic depth (tumor, reads supporting alternate allele) for inclusion
                        in report (default: None)
  --tumor_af_min TUMOR_AF_MIN
                        Minimum allelic fraction (tumor) for inclusion in report (default: None)
  --control_dp_min CONTROL_DP_MIN
                        Minimum sequencing depth (control) for inclusion in report (default: None)
  --control_ad_max CONTROL_AD_MAX
                        Maximum allelic depth (control, reads supporting alternate allele) for inclusion
                        in report (default: None)
  --control_af_max CONTROL_AF_MAX
                        Maximum allelic fraction (control) for inclusion in report (default: None)

Tumor-only filtering options:
  --pon_vcf PON_VCF     VCF file with germline calls from Panel of Normals (PON) - blacklisted variants
                        (default: None)
  --gnomad_popmax_af_tolerated GNOMAD_POPMAX_AF_TOLERATED
                        Exclude variants in tumor (SNVs/InDels, tumor-only mode) with gnomAD popmax
                        MAF greater than this value (default: 0.001)
  --exclude_pon         Exclude variants occurring in PoN (if --pon_vcf provided), default: False
  --exclude_likely_hom_germline
                        Exclude likely homozygous germline variants (allelic fraction of 1.0 for
                        alternate allele in tumor), default: False
  --exclude_likely_het_germline
                        Exclude likely heterozygous germline variants (0.4-0.6 allelic fraction, AND
                        presence in dbSNP + gnomAD, AND not in COSMIC/TCGA), default: False
  --exclude_clinvar_germline
                        Exclude variants found in ClinVar (germline origin only), default: False
  --exclude_dbsnp_nonsomatic
                        Exclude variants found in dbSNP (except those in ClinVar (somatic), TCGA, or
                        COSMIC), default: False
  --exclude_nonexonic   Exclude non-exonic variants, default: False

Tumor mutational burden (TMB) and MSI options:
  --estimate_tmb        Estimate tumor mutational burden from total somatic mutations and target region
                        size, default: False
  --tmb_display {coding_and_silent,coding_non_silent,missense_only}
                        Type of TMB measure to show in report, default: coding_and_silent
  --tmb_dp_min TMB_DP_MIN
                        Minimum sequencing depth (tumor) for TMB calculation (default: None)
  --tmb_af_min TMB_AF_MIN
                        Minimum allelic fraction (tumor) for TMB calculation (default: None)
  --tmb_ad_min TMB_AD_MIN
                        Minimum allelic depth (tumor) for TMB calculation (default: None)
  --estimate_msi        Predict microsatellite instability status from somatic mutation patterns,
                        default: False

Mutational signature options:
  --estimate_signatures
                        Estimate relative contributions of reference mutational signatures (re-fitting),
                        default: False
  --min_mutations_signatures MIN_MUTATIONS_SIGNATURES
                        Minimum number of SNVs required for signature re-fitting (SBS)
                        (default: 200, minimum n = 100)
  --all_reference_signatures
                        Use all reference mutational signatures (SBS) during re-fitting rather than
                        only those attributed to the tumor type (default: False)
  --include_artefact_signatures
                        Include sequencing artefacts in the reference signature collection
                        (default: False)
  --prevalence_reference_signatures PREVALENCE_REFERENCE_SIGNATURES
                        Minimum tumor-type prevalence (in percent) of reference signatures to be
                        included in refitting (default: 0.1)

Somatic CNA analysis options:
  --cna_threshold_mode {absolute,relative,combined}
                        Thresholding mode for CNA tier assignment: 'absolute' (total copy number),
                        'relative' (fold-change over tumor ploidy), or 'combined' (both criteria must
                        be met) (default: absolute)
  --cna_amp_threshold_absolute CNA_AMP_THRESHOLD_ABSOLUTE
                        Absolute total copy number threshold for amplifications (default: 5)
  --cna_amp_threshold_relative CNA_AMP_THRESHOLD_RELATIVE
                        Relative fold-change over tumor ploidy for amplifications (default: 2.5)
  --cna_gain_threshold_absolute CNA_GAIN_THRESHOLD_ABSOLUTE
                        Absolute total copy number threshold for gains (default: 3)
  --cna_gain_threshold_relative CNA_GAIN_THRESHOLD_RELATIVE
                        Relative fold-change over tumor ploidy for gains (default: 1.5)
  --cna_del_threshold_absolute CNA_DEL_THRESHOLD_ABSOLUTE
                        Absolute total copy number at or below which (but above zero) a segment is
                        considered a heterozygous deletion (default: 1)
  --cna_del_threshold_relative CNA_DEL_THRESHOLD_RELATIVE
                        Relative fold-change below tumor ploidy for heterozygous deletions
                        (default: 0.5)
  --cna_transcript_overlap_pct CNA_TRANSCRIPT_OVERLAP_PCT
                        Minimum mean percent overlap between copy number segment and gene transcripts
                        for reporting of gains/losses in oncogenes/tumor suppressors (default: 50)

RNA expression and fusion options:
  --fusion_min_split_reads FUSION_MIN_SPLIT_READS
                        Minimum number of split reads supporting a fusion event
                        (default: 3, minimum: 2)
  --expression_sim      Compare expression profile of tumor sample to known expression profiles
                        (default: False)
  --expression_sim_db EXPRESSION_SIM_DB
                        Comma-separated string of databases for RNA expression similarity analysis,
                        default: tcga,depmap,treehouse

Germline variant options:
  --input_cpsr INPUT_CPSR
                        CPSR-classified germline calls
                        (file '<cpsr_sample_id>.cpsr.<genome_assembly>.classification.tsv.gz')
  --input_cpsr_yaml INPUT_CPSR_YAML
                        CPSR YAML configuration file
                        (file '<cpsr_sample_id>.cpsr.<genome_assembly>.conf.yaml')
  --cpsr_ignore_vus     Do not show variants of uncertain significance (VUS) in the germline section
                        of the HTML report (default: False)

Biomarker and tiering options:
  --oncokb_api_token ONCOKB_API_TOKEN
                        OncoKB API token for querying the OncoKB precision oncology knowledge base
                        (default: None)
  --oncokb_oncotree_code ONCOKB_ONCOTREE_CODE
                        OncoTree code specifying the tumor type for OncoKB queries (default: None)
  --oncokb_exclusive    Limit biomarker reporting to OncoKB only - skip CIViC and CGI sources
                        (default: False)
  --oncokb_maf_query_all
                        Query OncoKB for all variant classes, including non-coding variants
                        (IGR, Intron, UTR, flanking regions). By default, non-coding variants
                        are filtered out before OncoKB annotation to reduce processing time.
                        Intended for TARGETED/WES assays only - enabling for WGS may result in
                        very long MafAnnotator.py runtimes (default: False)

Other options:
  --force_overwrite     Force overwrite of existing result files (default: False)
  --version             Show program's version number and exit
  --no_reporting        Run VEP/vcfanno annotation only; skip tier assignment, MSI, TMB, signatures,
                        and report generation (default: False)
  --no_html             Do not generate HTML report (default: False)
  --debug               Print full commands to log
  --pcgrr_conda PCGRR_CONDA
                        pcgrr conda environment name (default: pcgrr)
```

## Example run

The *examples* folder contains input VCF files from two tumor samples
sequenced within TCGA. A molecular interpretation report for a
colorectal tumor sample can be generated by running the following
command in your terminal window (this assumes you have installed PCGR
through Conda, and that you have downloaded the necessary reference data
files, VEP cache, see
[Installation](https://sigven.github.io/pcgr/dev/articles/installation.md)):

``` bash
$ (base) conda activate pcgr
$ (pcgr)
pcgr \
    --refdata_dir /Users/you/dir2/data \
    --vep_dir /Users/you/dir3/vep/.vep \
    --output_dir /Users/you/dir4/pcgr_outputs \
    --sample_id T001-COAD \
    --genome_assembly grch37 \
    --input_vcf /Users/you/pcgr/examples/T001-COAD.grch37.vcf.gz \
    --tumor_dp_tag TDP \
    --tumor_af_tag TVAF \
    --tumor_site 9 \
    --input_cna /Users/you/pcgr/examples/T001-COAD.grch37.cna.tsv \
    --input_rna_fusion /Users/you/pcgr/examples/T001-COAD.grch37.fusions.tsv \
    --sex MALE \
    --tumor_purity 0.9 \
    --tumor_ploidy 2.0 \
    --assay WES \
    --vcf2maf \
    --estimate_signatures \
    --estimate_msi \
    --estimate_tmb \
    --force_overwrite
```

This command will run the Conda-based PCGR workflow and produce the
following files in the */Users/you/dir4/pcgr_outputs* folder:

| N | File | Description |
|----|----|----|
| 1 | **\<sample_id\>.pcgr.grch37.html** | An interactive HTML report for clinical interpretation (quarto-based) |
| 2 | **\<sample_id\>.pcgr.grch37.xlsx** | An excel workbook with multiple sheets of annotations (Assay & sample info/SNVs & InDels/CNAs/biomarkers/TMB/MSI), suitable for aggregation analysis across multiple samples |
| 3 | **\<sample_id\>.pcgr.grch37.vcf.gz (.tbi)** | Bgzipped VCF file with rich set of variant annotations to support interpretation |
| 4 | **\<sample_id\>.pcgr.grch37.pass.vcf.gz (.tbi)** | Bgzipped VCF file with rich set of variant annotations to support interpretation (PASS variants only) |
| 5 | **\<sample_id\>.pcgr.grch37.pass.tsv.gz** | Compressed vcf2tsv-converted file with rich set of variant annotations to support interpretation |
| 6 | **\<sample_id\>.pcgr.grch37.conf.yaml** | PCGR configuration data file (YAML), as generated by pre-reporting annotation workflow |
| 7 | **\<sample_id\>.pcgr.grch37.cna_segment.tsv.gz** | Tab-separated values file with raw copy number segments, per affected transcript |
| 8 | **\<sample_id\>.pcgr.grch37.cna_gene_ann.tsv.gz** | Tab-separated values file with annotated gene copy number alteration events, including actionability assessment |
| 9 | **\<sample_id\>.pcgr.grch37.fusion_ann.tsv.gz** | Tab-separated values file with annotated RNA fusion events, including actionability assessment |
| 10 | **\<sample_id\>.pcgr.grch37.tmb.tsv** | Tab-separated values file with information on tumor mutational burden (TMB) estimates |
| 11 | **\<sample_id\>.pcgr.grch37.maf** | Annotated SNVs/InDels, converted to the Mutation Annotation Format (MAF) |
| 12 | **\<sample_id\>.pcgr.grch37.msigs.tsv.gz** | Tab-separated values file with information on the contribution of mutational signatures in the tumor sample |
| 13 | **\<sample_id\>.pcgr.grch37.snv_indel_ann.tsv.gz** | Tab-separated values file with key SNV/InDel variant annotations, including oncogenicity and clinical actionability assessment |
