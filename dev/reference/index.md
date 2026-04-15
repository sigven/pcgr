# Package index

## All functions

- [`af_distribution()`](https://sigven.github.io/pcgr/dev/reference/af_distribution.md)
  : Function that plots a histogram of the the variant allelic support
  (tumor)
- [`append_annotation_links()`](https://sigven.github.io/pcgr/dev/reference/append_annotation_links.md)
  : Function that appends multiple HTML annotation links to variant
  identifiers e.g. COSMIC, CLINVAR, REFSEQ etc
- [`append_cancer_association_ranks()`](https://sigven.github.io/pcgr/dev/reference/append_cancer_association_ranks.md)
  : Function that appends cancer gene evidence links
- [`append_cancer_gene_evidence()`](https://sigven.github.io/pcgr/dev/reference/append_cancer_gene_evidence.md)
  : Function that appends cancer gene evidence links
- [`append_dbmts_var_link()`](https://sigven.github.io/pcgr/dev/reference/append_dbmts_var_link.md)
  : Function that adds miRNA target annotations (dbMTS) to genetic
  variant identifiers
- [`append_dbnsfp_var_link()`](https://sigven.github.io/pcgr/dev/reference/append_dbnsfp_var_link.md)
  : Function that assigns HTML links to dbNSFP prediction entries
- [`append_drug_var_link()`](https://sigven.github.io/pcgr/dev/reference/append_drug_var_link.md)
  : Function that adds link to targeted drugs (on and off-label) for a
  list of variants and associated targeted
- [`append_gwas_citation_phenotype()`](https://sigven.github.io/pcgr/dev/reference/append_gwas_citation_phenotype.md)
  : Function that adds GWAS citation/phenotype to GWAS hit found through
  PCGR annotation
- [`append_oncogenicity_docs()`](https://sigven.github.io/pcgr/dev/reference/append_oncogenicity_docs.md)
  : Function that adds oncogenicity documentation from codes
- [`append_targeted_drug_annotations()`](https://sigven.github.io/pcgr/dev/reference/append_targeted_drug_annotations.md)
  : Function that adds link to targeted drugs (on and off-label) for a
  list of variants and associated targeted
- [`append_tcga_var_link()`](https://sigven.github.io/pcgr/dev/reference/append_tcga_var_link.md)
  : Function that adds TCGA annotations (cohort, frequency etc.) to
  variant identifiers
- [`append_tfbs_annotation()`](https://sigven.github.io/pcgr/dev/reference/append_tfbs_annotation.md)
  : Function that adds TFBS annotations (dbMTS) to genetic variant
  identifiers
- [`assign_amp_asco_tiers()`](https://sigven.github.io/pcgr/dev/reference/assign_amp_asco_tiers.md)
  : Function that assigns tier classifications to somatic CNA segments
  and SNVs/InDels, based on the presence of biomarker evidence found in
  the variant set
- [`assign_germline_popfreq_status()`](https://sigven.github.io/pcgr/dev/reference/assign_germline_popfreq_status.md)
  : Function that sets STATUS_POPFREQ_1KG_ABOVE_TOLERATED/
  STATUS_POPFREQ_GNOMAD_ABOVE_TOLERATED to TRUE for variants if any
  population frequency exceeds max_tolerated_af
- [`assign_mutation_type()`](https://sigven.github.io/pcgr/dev/reference/assign_mutation_type.md)
  : Function that assigns one of six mutation types to a list of
  mutations
- [`assign_somatic_classification()`](https://sigven.github.io/pcgr/dev/reference/assign_somatic_classification.md)
  : Function that assigns a SOMATIC_CLASSIFICATION to variants based on
  evidence found in variant set, potentially limited by user-defined
  options
- [`assign_somatic_germline_evidence()`](https://sigven.github.io/pcgr/dev/reference/assign_somatic_germline_evidence.md)
  : Function that appends several tags denoting evidence for
  somatic/germline status of variants
- [`biomarker_evidence`](https://sigven.github.io/pcgr/dev/reference/biomarker_evidence.md)
  : Fixed data types/categories used for biomarker evidence, e.g.
  'types','levels' etc.
- [`cancer_phenotypes_regex`](https://sigven.github.io/pcgr/dev/reference/cancer_phenotypes_regex.md)
  : Regular expression of terms indicative of cancer-related phenotypes
  and syndromes
- [`check_common_colnames()`](https://sigven.github.io/pcgr/dev/reference/check_common_colnames.md)
  : Function that checks whether a set of column names are present in
  two different data frames
- [`check_file_exists()`](https://sigven.github.io/pcgr/dev/reference/check_file_exists.md)
  : Function that checks the existence of a file
- [`clinvar_germline_status()`](https://sigven.github.io/pcgr/dev/reference/clinvar_germline_status.md)
  : Function that assigns a logical to STATUS_CLINVAR_GERMLINE based on
  whether a ClinVar entry of germline origin is found for a given
  variant (for entries in a data frame)
- [`color_palette`](https://sigven.github.io/pcgr/dev/reference/color_palette.md)
  : Color encodings for report elements of PCGR/CPSR
- [`cosmic_sbs_signatures`](https://sigven.github.io/pcgr/dev/reference/cosmic_sbs_signatures.md)
  : List of COSMIC reference mutational signatures (SBS, v3.4)
- [`cosmic_somatic_status()`](https://sigven.github.io/pcgr/dev/reference/cosmic_somatic_status.md)
  : Function that assigns a logical (STATUS_COSMIC) reflecting whether a
  variant co-incides with an entry in COSMIC (germline)
- [`data_coltype_defs`](https://sigven.github.io/pcgr/dev/reference/data_coltype_defs.md)
  : List of coltype definitions for input files to pcgrr (e.g.
  VCF-converted TSV, CNA TVS etc.)
- [`dbsnp_germline_status()`](https://sigven.github.io/pcgr/dev/reference/dbsnp_germline_status.md)
  : Function that assigns a logical (STATUS_DBSNP) reflecting whether a
  variant co-incides with an entry in dbSNP (germline)
- [`deduplicate_eitems()`](https://sigven.github.io/pcgr/dev/reference/deduplicate_eitems.md)
  : Function that removes redundancy in variant evidence items (i.e. if
  a variant is assicated with evidence at the codon level, evidence at
  the exon/gene level is ignored)
- [`detect_vcf_sample_name()`](https://sigven.github.io/pcgr/dev/reference/detect_vcf_sample_name.md)
  : A function that detects whether the sample name in variant data
  frame is unique (as present in column name VCF_SAMPLE_ID), throws an
  error if multiple sample names are present for the CPSR workflow
- [`df_string_replace()`](https://sigven.github.io/pcgr/dev/reference/df_string_replace.md)
  : Function that performs stringr::str_replace on strings of multiple
  string columns of a dataframe
- [`dt_display`](https://sigven.github.io/pcgr/dev/reference/dt_display.md)
  : DT Display
- [`effect_prediction_algos`](https://sigven.github.io/pcgr/dev/reference/effect_prediction_algos.md)
  : List of URLs for a range of variant effect prediction algorithms
- [`exclude_non_chrom_variants()`](https://sigven.github.io/pcgr/dev/reference/exclude_non_chrom_variants.md)
  : Function that excludes genomic aberrations from non-nuclear
  chromosomes
- [`expand_biomarker_items()`](https://sigven.github.io/pcgr/dev/reference/expand_biomarker_items.md)
  : Function that expands biomarker evidence items with variant
  annotations
- [`export_quarto_evars()`](https://sigven.github.io/pcgr/dev/reference/export_quarto_evars.md)
  : Export Quarto Environment Variables
- [`filter_eitems_by_site()`](https://sigven.github.io/pcgr/dev/reference/filter_eitems_by_site.md)
  : Function that filters clinical evidence items by tumor type/primary
  site
- [`filter_maf_file()`](https://sigven.github.io/pcgr/dev/reference/filter_maf_file.md)
  : Function that takes a MAF file generated with vcf2maf and filters
  out variants that are presumably germline (tumor-only run)
- [`filter_read_support()`](https://sigven.github.io/pcgr/dev/reference/filter_read_support.md)
  : Function that filters variant set on depth and allelic fraction
  according to settings provided by user (tumor and control)
- [`generate_annotation_link()`](https://sigven.github.io/pcgr/dev/reference/generate_annotation_link.md)
  : A function that generates a HTML link for selected identifiers
  (DBSNP, COSMIC, CLINVAR, ENTREZ)
- [`generate_report()`](https://sigven.github.io/pcgr/dev/reference/generate_report.md)
  : Function that generates all contents of the cancer genome report
  (PCGR)
- [`generate_report_data_expression()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_expression.md)
  : Function that generates expression data for PCGR report
- [`generate_report_data_kataegis()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_kataegis.md)
  : Function that generates data frame with potential kataegis events
- [`generate_report_data_msi()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_msi.md)
  : Function that generates MSI prediction data for PCGR report
- [`generate_report_data_rainfall()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_rainfall.md)
  : Function that generates data for rainfall plot (mutation density
  along genome, considering SNVs only)
- [`generate_report_data_signatures()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_signatures.md)
  : Function that generates mutational signatures data for PCGR report
- [`generate_report_data_tmb()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_tmb.md)
  : Function that reads TSV file with TMB estimates from sample
- [`generate_report_data_trials()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_trials.md)
  : Function that retrieves relevant (interventional based on molecular
  target) clinical trials for a given tumor type
- [`generate_tier_tsv()`](https://sigven.github.io/pcgr/dev/reference/generate_tier_tsv.md)
  : Function that generates dense and tiered annotated variant datasets
- [`germline_filter_levels`](https://sigven.github.io/pcgr/dev/reference/germline_filter_levels.md)
  : Data frame with germline filtering criteria
- [`get_clin_assocs_cna()`](https://sigven.github.io/pcgr/dev/reference/get_clin_assocs_cna.md)
  : Function that retrieves clinical evidence items (CIVIC, CBMDB) for
  CNA aberrations
- [`get_dt_tables()`](https://sigven.github.io/pcgr/dev/reference/get_dt_tables.md)
  : Function that gathers data tables on actionable variants for display
  in report (tier 1 + tier 2)
- [`get_excel_sheets()`](https://sigven.github.io/pcgr/dev/reference/get_excel_sheets.md)
  : Function that produces the contents of sheets for an Excel report of
  PCGR output
- [`get_genome_obj()`](https://sigven.github.io/pcgr/dev/reference/get_genome_obj.md)
  : Get BSgenome Object
- [`get_oncogenic_cna_events()`](https://sigven.github.io/pcgr/dev/reference/get_oncogenic_cna_events.md)
  : Get oncogenic copy number events
- [`get_prevalent_site_signatures()`](https://sigven.github.io/pcgr/dev/reference/get_prevalent_site_signatures.md)
  : Function that retrieves prevalent signatures for a given tumor
  type/primary site Data is collected from COSMIC v3.4.
- [`get_tumor_only_filtering_criteria()`](https://sigven.github.io/pcgr/dev/reference/get_tumor_only_filtering_criteria.md)
  : Function that generates a string with filtering criteria for
  callsets coming from tumor-only sequencing
- [`get_valid_chromosomes()`](https://sigven.github.io/pcgr/dev/reference/get_valid_chromosomes.md)
  : Checks for valid chromosome names in data frame of variants
- [`get_variant_statistics()`](https://sigven.github.io/pcgr/dev/reference/get_variant_statistics.md)
  : Function that computes various variant statistics from a data frame
  with variant records
- [`het_af_germline_status()`](https://sigven.github.io/pcgr/dev/reference/het_af_germline_status.md)
  : Function that assigns a logical
  (STATUS_LIKELY_GERMLINE_HETEROZYGOUS) reflecting whether a variant is
  likely heterozygous (germline) - based on allelic fraction
  (VAF_TUMOR), presence in gnomAD and dbSNP, and no presence in TCGA and
  COSMIC
- [`hex_to_rgba()`](https://sigven.github.io/pcgr/dev/reference/hex_to_rgba.md)
  : Convert Hex Color to RGBA
- [`hom_af_status()`](https://sigven.github.io/pcgr/dev/reference/hom_af_status.md)
  : Function that assigns a logical (STATUS_LIKELY_GERMLINE_HOMOZYGOUS)
  reflecting whether a variant is likely homozygous (germline) - based
  on allelic fraction (VAF_TUMOR)
- [`immune_celltypes`](https://sigven.github.io/pcgr/dev/reference/immune_celltypes.md)
  : Data frame with immune cell types
- [`init_cna_vstats()`](https://sigven.github.io/pcgr/dev/reference/init_cna_vstats.md)
  : Function that initiates report element with CNA information
- [`init_expression_content()`](https://sigven.github.io/pcgr/dev/reference/init_expression_content.md)
  : Function that initiates report element with expression information
- [`init_germline_content()`](https://sigven.github.io/pcgr/dev/reference/init_germline_content.md)
  : Function that initiates report element with germline variant
  information (CPSR)
- [`init_kataegis_content()`](https://sigven.github.io/pcgr/dev/reference/init_kataegis_content.md)
  : Function that initiates report element with kataegis information
- [`init_m_signature_content()`](https://sigven.github.io/pcgr/dev/reference/init_m_signature_content.md)
  : Function that initiates report element with mutational signatures
  information
- [`init_msi_content()`](https://sigven.github.io/pcgr/dev/reference/init_msi_content.md)
  : Function that initiates report element with MSI classification
- [`init_rainfall_content()`](https://sigven.github.io/pcgr/dev/reference/init_rainfall_content.md)
  : Function that initiates report element with rainfall information
- [`init_report()`](https://sigven.github.io/pcgr/dev/reference/init_report.md)
  : Function that initiates PCGR/CPSR report object
- [`init_snv_indel_vstats()`](https://sigven.github.io/pcgr/dev/reference/init_snv_indel_vstats.md)
  : Function that initiates report element with SNV/InDel statistics
  information
- [`init_tmb_content()`](https://sigven.github.io/pcgr/dev/reference/init_tmb_content.md)
  : Function that initiates report element with TMB information
- [`init_tumor_only_content()`](https://sigven.github.io/pcgr/dev/reference/init_tumor_only_content.md)
  : Function that initiates report element with tumor-only information
- [`init_var_content()`](https://sigven.github.io/pcgr/dev/reference/init_var_content.md)
  : Function that initiates report element with variant data
- [`kataegis_detect()`](https://sigven.github.io/pcgr/dev/reference/kataegis_detect.md)
  : Function that detects kataegis events from a data frame with genomic
  cooordinates of mutations
- [`kataegis_input()`](https://sigven.github.io/pcgr/dev/reference/kataegis_input.md)
  : Function that detects kataegis events from a data frame with genomic
  cooordinates of mutations
- [`load_all_eitems()`](https://sigven.github.io/pcgr/dev/reference/load_all_eitems.md)
  : Function that loads all evidence items from CIViC and CGI, and
  combines them in a unified data.frame
- [`load_cpsr_classified_variants()`](https://sigven.github.io/pcgr/dev/reference/load_cpsr_classified_variants.md)
  : Function that reads CPSR-classified variants from a TSV file
- [`load_dna_variants()`](https://sigven.github.io/pcgr/dev/reference/load_dna_variants.md)
  : Function that reads and validates CNA or SNV/InDel TSV files file
  from PCGR/CPSR pre-report (Python) pipeline
- [`load_eitems()`](https://sigven.github.io/pcgr/dev/reference/load_eitems.md)
  : Function that loads specific set of clinical variant evidence items
  (CIViC + CGI) based on given parameters (mutation type, variant
  origin, tumor type etc)
- [`load_expression_csq()`](https://sigven.github.io/pcgr/dev/reference/load_expression_csq.md)
  : Load expression consequence settings
- [`load_expression_outliers()`](https://sigven.github.io/pcgr/dev/reference/load_expression_outliers.md)
  : Load expression outlier results
- [`load_expression_similarity()`](https://sigven.github.io/pcgr/dev/reference/load_expression_similarity.md)
  : Load expression similarity results
- [`load_reference_data()`](https://sigven.github.io/pcgr/dev/reference/load_reference_data.md)
  : Function that parses and loads reference data from files in the
  assembly-specific PCGR bundle directory
- [`load_somatic_cna()`](https://sigven.github.io/pcgr/dev/reference/load_somatic_cna.md)
  : Function that reads and validates fully annotated CNA data (segments
  and genes) from PCGR pre-reporting pipeline
- [`load_somatic_snv_indel()`](https://sigven.github.io/pcgr/dev/reference/load_somatic_snv_indel.md)
  : Function that reads and validates an annotated somatic SNV/InDel
  file from PCGR pre-reporting pipeline
- [`load_yaml()`](https://sigven.github.io/pcgr/dev/reference/load_yaml.md)
  : Function that loads YAML data with settings and file paths to
  annotated molecular profiles
- [`log4r_debug()`](https://sigven.github.io/pcgr/dev/reference/log4r_debug.md)
  : Write messages to logs at a given priority level
- [`log4r_fatal()`](https://sigven.github.io/pcgr/dev/reference/log4r_fatal.md)
  : Write messages to logs at a given priority level
- [`log4r_info()`](https://sigven.github.io/pcgr/dev/reference/log4r_info.md)
  : Write messages to logs at a given priority level
- [`log4r_warn()`](https://sigven.github.io/pcgr/dev/reference/log4r_warn.md)
  : Write messages to logs at a given priority level
- [`log_var_eitem_stats()`](https://sigven.github.io/pcgr/dev/reference/log_var_eitem_stats.md)
  : Function that logs the number of evidence items found, for different
  levels of resolution
- [`max_af_gnomad()`](https://sigven.github.io/pcgr/dev/reference/max_af_gnomad.md)
  : Function that assigns a maximum value to a variable (MAX_AF_GNOMAD)
  reflecting the maximum allele frequency for a given variant across
  gnomAD populations
- [`mkdir()`](https://sigven.github.io/pcgr/dev/reference/mkdir.md) :
  Create directory
- [`msi_indel_fraction_plot()`](https://sigven.github.io/pcgr/dev/reference/msi_indel_fraction_plot.md)
  : Function that plots the indel fraction for a given sample and
  contrasts this with the distribution for MSI-H/MSS samples from TCGA
- [`msi_indel_load_plot()`](https://sigven.github.io/pcgr/dev/reference/msi_indel_load_plot.md)
  : Function that plots the indel load for a given sample and contrasts
  this with the distribution for MSI-H/MSS samples from TCGA
- [`oncogenicity_criteria`](https://sigven.github.io/pcgr/dev/reference/oncogenicity_criteria.md)
  : Oncogenicity criteria (ClinGen/CGC/VICC)
- [`order_variants()`](https://sigven.github.io/pcgr/dev/reference/order_variants.md)
  : Function that orders genomic aberrations according to order of
  chromosomes and chromosomal position
- [`pkg_exists()`](https://sigven.github.io/pcgr/dev/reference/pkg_exists.md)
  : Does R Package Exist
- [`plot_cna_segments()`](https://sigven.github.io/pcgr/dev/reference/plot_cna_segments.md)
  : Plot allele-specific copy number segments
- [`plot_filtering_stats_exonic()`](https://sigven.github.io/pcgr/dev/reference/plot_filtering_stats_exonic.md)
  : Function that generates a pie chart for exonic/non-exonic variant
  statistics (for callsets coming from tumor-only sequencing)
- [`plot_filtering_stats_germline()`](https://sigven.github.io/pcgr/dev/reference/plot_filtering_stats_germline.md)
  : Function that makes input data for an UpSet plot
  (filtering/intersection results) for the somatic-germline
  classification procedure
- [`plot_signature_contributions()`](https://sigven.github.io/pcgr/dev/reference/plot_signature_contributions.md)
  : Function that makes plots of mutational signature contributions in a
  given sample (both ggplot and plotly)
- [`plot_tmb_primary_site_tcga()`](https://sigven.github.io/pcgr/dev/reference/plot_tmb_primary_site_tcga.md)
  : Function that makes a plot with TMB boxplots for reference cohorts,
  highlighting the TMB estimate for a given sample and the
  cohort/primary site of interest
- [`plot_value_boxes()`](https://sigven.github.io/pcgr/dev/reference/plot_value_boxes.md)
  : Function that plots four value boxes with the most important
  findings in the cancer genome
- [`plotly_pie_chart()`](https://sigven.github.io/pcgr/dev/reference/plotly_pie_chart.md)
  : Plotly Pie Chart - variant statistics
- [`pon_status()`](https://sigven.github.io/pcgr/dev/reference/pon_status.md)
  : Function that assigns a logical (STATUS_PON) reflecting whether a
  variant is co-inciding with a variant present in a panel-of-normals
  database (PANEL_OF_NORMALS column is TRUE)
- [`predict_msi_status()`](https://sigven.github.io/pcgr/dev/reference/predict_msi_status.md)
  : Function that predicts MSI status based on fraction of indels among
  calls
- [`qc_var_eitems()`](https://sigven.github.io/pcgr/dev/reference/qc_var_eitems.md)
  : Function that matches variants to evidence items
- [`remove_cols_from_df()`](https://sigven.github.io/pcgr/dev/reference/remove_cols_from_df.md)
  : Function that removes column(s) from data frame
- [`sort_chromosomal_segments()`](https://sigven.github.io/pcgr/dev/reference/sort_chromosomal_segments.md)
  : Function that sorts chromosomal segments according to chromosome and
  chromosomal start/end position
- [`strip_html()`](https://sigven.github.io/pcgr/dev/reference/strip_html.md)
  : Strip HTML tags
- [`structure_var_eitems()`](https://sigven.github.io/pcgr/dev/reference/structure_var_eitems.md)
  : Function that structures variant evidence items according to
  strength of evidence
- [`tcga_cohorts`](https://sigven.github.io/pcgr/dev/reference/tcga_cohorts.md)
  : Data frame with all TCGA cohorts
- [`tcga_somatic_status()`](https://sigven.github.io/pcgr/dev/reference/tcga_somatic_status.md)
  : Function that assigns a logical (STATUS_TCGA_SOMATIC) reflecting
  whether a variant co-incides with an entry in TCGA (somatic)
- [`tier_af_distribution()`](https://sigven.github.io/pcgr/dev/reference/tier_af_distribution.md)
  : Function that plots a histogram of the the variant allelic support
  (tumor) - grouped by tiers
- [`tsv_cols`](https://sigven.github.io/pcgr/dev/reference/tsv_cols.md)
  : TSV columns
- [`update_report()`](https://sigven.github.io/pcgr/dev/reference/update_report.md)
  : Function that updates a PCGR/CPSR report object structure
- [`vaf_plot()`](https://sigven.github.io/pcgr/dev/reference/vaf_plot.md)
  : Function that generates a VAF distribution plot for a given PCGR
  report object
- [`variant_db_url`](https://sigven.github.io/pcgr/dev/reference/variant_db_url.md)
  : List of URLS and variant identifiers for variant/gene/protein domain
  databases
- [`variant_stats_report()`](https://sigven.github.io/pcgr/dev/reference/variant_stats_report.md)
  : Function that generate stats for a given variant set, considering
  number of variants/genes affected across tiers, types of variants ()
- [`write_processed_vcf()`](https://sigven.github.io/pcgr/dev/reference/write_processed_vcf.md)
  : Function that writes a VCF intended for mutational signature
  analysis
- [`write_report_excel()`](https://sigven.github.io/pcgr/dev/reference/write_report_excel.md)
  : Function that writes key datasets of PCGR object to an Excel
  workbook
- [`write_report_quarto_html()`](https://sigven.github.io/pcgr/dev/reference/write_report_quarto_html.md)
  : Function that writes contents of PCGR object to an HTML report
  (quarto-based)
- [`write_report_tsv()`](https://sigven.github.io/pcgr/dev/reference/write_report_tsv.md)
  : Function that writes contents of PCGR object to various output
  formats (Rmarkdown/flexdashboard HTML reports, JSON, tab-separated
  etc)
