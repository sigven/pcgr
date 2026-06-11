If bulk RNA-seq data is provided, PCGR reports a gene expression section covering three analyses:
outlier detection, sample similarity, and immune cell fraction estimation.

Expression input is expected as transcript-level TPM values (e.g. from Salmon or kallisto).
PCGR aggregates these to gene-level TPM (taking the maximum transcript TPM per gene) and
converts values to log2(TPM + 0.001) for comparability with reference data.

The sample's gene-level expression profile is compared against the TCGA project cohort that most
closely matches the user-specified primary tumor site. For each protein-coding gene, a
percentile rank is computed relative to the reference cohort, and genes with extreme
percentile values (high or low) are flagged as expression outliers. Outlier detection
is currently limited to tumor types with a corresponding TCGA cohort.

Spearman correlations are computed between the input sample's expression profile and profiles
from one or more reference collections:

- **[TCGA](https://portal.gdc.cancer.gov/)** — across all available tumor type cohorts
- **[DepMap](https://depmap.org/portal/)** — cancer cell line expression profiles
- **[Treehouse Childhood Cancer Data Initiative](https://treehousegenomics.soe.ucsc.edu/)** — pediatric tumor expression profiles (where available)

The top-correlated reference samples are reported along with their metadata (tumor type,
tissue of origin), providing an independent check on the sample's transcriptomic identity.

#### Immune cell fraction estimation

Immune cell deconvolution is performed on the TPM profile to estimate the relative fractions
of immune cell types infiltrating the tumor microenvironment.
