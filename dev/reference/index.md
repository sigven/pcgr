# Package index

## All functions

- [`actionability_doc_note()`](https://sigven.github.io/pcgr/dev/reference/actionability_doc_note.md)
  : Get documentation string for AMP/ASCO/CAP clinical actionability
  tiers
- [`af_distribution()`](https://sigven.github.io/pcgr/dev/reference/af_distribution.md)
  : Function that plots a histogram of the the variant allelic support
  (tumor)
- [`append_alteration_name()`](https://sigven.github.io/pcgr/dev/reference/append_alteration_name.md)
  : Function that appends informative alteration names to sample
  variants based on gene symbol, consequence and HGVS annotations
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
- [`append_protein_domains()`](https://sigven.github.io/pcgr/dev/reference/append_protein_domains.md)
  : Function that adds protein domain annotations (PFAM) for a list of
  variants and associated targeted
- [`append_styled_cna_vclass()`](https://sigven.github.io/pcgr/dev/reference/append_styled_cna_vclass.md)
  : Function that styles CNA variant classes (amplification, gain,
  hetdel, homdel/hemdel) with more informative names for display in PCGR
  reports and other outputs
- [`append_targeted_drug_annotations()`](https://sigven.github.io/pcgr/dev/reference/append_targeted_drug_annotations.md)
  : Function that adds link to targeted drugs (on and off-label) for a
  list of variants and associated targeted
- [`append_tcga_var_link()`](https://sigven.github.io/pcgr/dev/reference/append_tcga_var_link.md)
  : Function that adds TCGA annotations (cohort, frequency etc.) to
  variant identifiers
- [`append_tfbs_annotation()`](https://sigven.github.io/pcgr/dev/reference/append_tfbs_annotation.md)
  : Function that adds TFBS annotations (dbMTS) to genetic variant
  identifiers
- [`assign_amp_asco_cap_tiers()`](https://sigven.github.io/pcgr/dev/reference/assign_amp_asco_cap_tiers.md)
  : AMP/ASCO/CAP tier classification for somatic variants in cancer
- [`assign_bm_tier_support_ttagnostic()`](https://sigven.github.io/pcgr/dev/reference/assign_bm_tier_support_ttagnostic.md)
  : For all biomarker evidence items with specified confidence level,
  assign whether these are tier-defining or providing additional
  support - tumor-type agnostic query
- [`assign_bm_tier_support_ttspecific()`](https://sigven.github.io/pcgr/dev/reference/assign_bm_tier_support_ttspecific.md)
  : For all biomarker evidence items with specified confidence level,
  assign whether these are tier-defining or providing additional
  support - tumor-type specific query
- [`assign_germline_popfreq_status()`](https://sigven.github.io/pcgr/dev/reference/assign_germline_popfreq_status.md)
  : Function that sets gnomAD_AF_ABOVE_TOLERATED to TRUE for variants if
  any gnomAD population frequency exceeds max_tolerated_af
- [`assign_germline_popfreq_status_old()`](https://sigven.github.io/pcgr/dev/reference/assign_germline_popfreq_status_old.md)
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
- [`assign_variant_tiers_cna()`](https://sigven.github.io/pcgr/dev/reference/assign_variant_tiers_cna.md)
  : Assign tiers of clinical significance (AMP/ASCO/CAP framework) to
  somatic CNAs
- [`assign_variant_tiers_fusion()`](https://sigven.github.io/pcgr/dev/reference/assign_variant_tiers_fusion.md)
  : Assign tiers of clinical significance (AMP/ASCO/CAP framework) to
  RNA fusions
- [`assign_variant_tiers_snv_indel()`](https://sigven.github.io/pcgr/dev/reference/assign_variant_tiers_snv_indel.md)
  : Assign tiers of clinical significance (AMP/ASCO/CAP framework) to
  somatic SNVs/InDels
- [`assign_variant_top_tiers_ttagnostic()`](https://sigven.github.io/pcgr/dev/reference/assign_variant_top_tiers_ttagnostic.md)
  : Assign top tiers of clinical significance (AMP/ASCO/CAP framework)
  to variants based on data from biomarker evidence items (tumor-type
  agnostic query)
- [`assign_variant_top_tiers_ttspecific()`](https://sigven.github.io/pcgr/dev/reference/assign_variant_top_tiers_ttspecific.md)
  : Assign top tiers of clinical significance (AMP/ASCO/CAP framework)
  to variants based on data from biomarker evidence items (tumor
  type-specific query)
- [`biomarker_evidence`](https://sigven.github.io/pcgr/dev/reference/biomarker_evidence.md)
  : Fixed data types and levels for biomarker evidence items
- [`bm_categories`](https://sigven.github.io/pcgr/dev/reference/bm_categories.md)
  : Named list of biomarker categories used for variant classification
- [`bm_evidence`](https://sigven.github.io/pcgr/dev/reference/bm_evidence.md)
  : Fixed data types/categories used for biomarker evidence, e.g.
  'types','levels' etc.
- [`bp_junction_transcript_overlap()`](https://sigven.github.io/pcgr/dev/reference/bp_junction_transcript_overlap.md)
  : Find transcripts covering given splice junction breakpoints
- [`build_oncogenicity_col_defs()`](https://sigven.github.io/pcgr/dev/reference/build_oncogenicity_col_defs.md)
  : Build column definitions for oncogenicity table with category-aware
  styling
- [`build_rt_row_details()`](https://sigven.github.io/pcgr/dev/reference/build_rt_row_details.md)
  : JS function for reactable row details with exclusions and card
  styling (Generated by Claude Opus 4.6 with some manual tweaks)
- [`build_twohit_display_data()`](https://sigven.github.io/pcgr/dev/reference/build_twohit_display_data.md)
  : Build display data for potential two-hit events (nested reactable)
- [`callout_biomarker_scope()`](https://sigven.github.io/pcgr/dev/reference/callout_biomarker_scope.md)
  : Emit the shared "biomarker types and report scope" callout note
- [`cancer_phenotypes_regex`](https://sigven.github.io/pcgr/dev/reference/cancer_phenotypes_regex.md)
  : Regular expression of terms indicative of cancer-related phenotypes
  and syndromes
- [`check_common_colnames()`](https://sigven.github.io/pcgr/dev/reference/check_common_colnames.md)
  : Function that checks whether a set of column names are present in
  two different data frames
- [`check_file_exists()`](https://sigven.github.io/pcgr/dev/reference/check_file_exists.md)
  : Function that checks the existence of a file
- [`clean_gnomad_annotations()`](https://sigven.github.io/pcgr/dev/reference/clean_gnomad_annotations.md)
  : Clean gnomAD VCF annotations
- [`clean_oncokb_evidence()`](https://sigven.github.io/pcgr/dev/reference/clean_oncokb_evidence.md)
  : Clean extracted evidence data from OncoKB
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
- [`detect_vcf_sample_name()`](https://sigven.github.io/pcgr/dev/reference/detect_vcf_sample_name.md)
  : A function that detects whether the sample name in variant data
  frame is unique (as present in column name VCF_SAMPLE_ID), throws an
  error if multiple sample names are present for the CPSR workflow
- [`df_string_replace()`](https://sigven.github.io/pcgr/dev/reference/df_string_replace.md)
  : Function that performs stringr::str_replace on strings of multiple
  string columns of a dataframe
- [`effect_prediction_algos`](https://sigven.github.io/pcgr/dev/reference/effect_prediction_algos.md)
  : List of URLs for a range of variant effect prediction algorithms
- [`exclude_non_chrom_variants()`](https://sigven.github.io/pcgr/dev/reference/exclude_non_chrom_variants.md)
  : Function that excludes genomic aberrations from non-nuclear
  chromosomes
- [`exonic_filter_levels`](https://sigven.github.io/pcgr/dev/reference/exonic_filter_levels.md)
  : Exonic filter levels
- [`export_quarto_evars()`](https://sigven.github.io/pcgr/dev/reference/export_quarto_evars.md)
  : Export Quarto Environment Variables
- [`expression_doc_note()`](https://sigven.github.io/pcgr/dev/reference/expression_doc_note.md)
  : Get documentation string for RNA expression analysis (outlier,
  similarity)
- [`extract_complete_annotation()`](https://sigven.github.io/pcgr/dev/reference/extract_complete_annotation.md)
  : Extract mutation effect information
- [`extract_diagnostic_evidence()`](https://sigven.github.io/pcgr/dev/reference/extract_diagnostic_evidence.md)
  : Extract diagnostic implications from OncoKB annotation
- [`extract_prognostic_evidence()`](https://sigven.github.io/pcgr/dev/reference/extract_prognostic_evidence.md)
  : Extract prognostic implications from OncoKB annotation
- [`extract_therapeutic_evidence()`](https://sigven.github.io/pcgr/dev/reference/extract_therapeutic_evidence.md)
  : Extract therapeutic evidence items from OncoKB annotation
- [`fetch_oncokb_cna_annotation()`](https://sigven.github.io/pcgr/dev/reference/fetch_oncokb_cna_annotation.md)
  : Fetch OncoKB annotation for copy number alteration
- [`fetch_oncokb_fusion_annotation()`](https://sigven.github.io/pcgr/dev/reference/fetch_oncokb_fusion_annotation.md)
  : Fetch OncoKB annotation for gene fusion
- [`fetch_oncokb_genomic_annotation()`](https://sigven.github.io/pcgr/dev/reference/fetch_oncokb_genomic_annotation.md)
  : Fetch OncoKB annotation for SNV/InDel via genomic change
- [`fetch_oncokb_hgvsp_annotation()`](https://sigven.github.io/pcgr/dev/reference/fetch_oncokb_hgvsp_annotation.md)
  : Fetch OncoKB annotation for SNV/InDel via protein change
- [`filter_maf_file()`](https://sigven.github.io/pcgr/dev/reference/filter_maf_file.md)
  : Function that takes a MAF file generated with vcf2maf and filters
  out variants that are presumably germline (tumor-only run)
- [`filter_read_support()`](https://sigven.github.io/pcgr/dev/reference/filter_read_support.md)
  : Function that filters variant set on depth and allelic fraction
  according to settings provided by user (tumor and control)
- [`fusion_doc_note()`](https://sigven.github.io/pcgr/dev/reference/fusion_doc_note.md)
  : Get documentation string for RNA fusion analysis
- [`generate_annotation_link()`](https://sigven.github.io/pcgr/dev/reference/generate_annotation_link.md)
  : A function that generates a HTML link for selected identifiers
  (DBSNP, COSMIC, CLINVAR, ENTREZ)
- [`generate_report()`](https://sigven.github.io/pcgr/dev/reference/generate_report.md)
  : Function that generates all contents of the cancer genome report
  (PCGR)
- [`generate_report_data_expression()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_expression.md)
  : Function that generates expression data for PCGR report
- [`generate_report_data_fusion()`](https://sigven.github.io/pcgr/dev/reference/generate_report_data_fusion.md)
  : Function that generates fusion data for PCGR report
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
- [`generate_tier_tsv()`](https://sigven.github.io/pcgr/dev/reference/generate_tier_tsv.md)
  : Function that generates dense and tiered annotated variant datasets
- [`germline_filter_levels`](https://sigven.github.io/pcgr/dev/reference/germline_filter_levels.md)
  : Data frame with germline filtering criteria
- [`get_data_versions_sheet()`](https://sigven.github.io/pcgr/dev/reference/get_data_versions_sheet.md)
  : Build the DATA_VERSIONS sheet for the Excel workbook
- [`get_druggable_fusion_partner()`](https://sigven.github.io/pcgr/dev/reference/get_druggable_fusion_partner.md)
  : Identify druggable fusion partners
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
- [`get_settings_sheet()`](https://sigven.github.io/pcgr/dev/reference/get_settings_sheet.md)
  : Build the SETTINGS sheet for the Excel workbook
- [`get_tumor_only_filtering_criteria()`](https://sigven.github.io/pcgr/dev/reference/get_tumor_only_filtering_criteria.md)
  : Function that generates a string with filtering criteria for
  callsets coming from tumor-only sequencing
- [`get_valid_chromosomes()`](https://sigven.github.io/pcgr/dev/reference/get_valid_chromosomes.md)
  : Checks for valid chromosome names in data frame of variants
- [`grpmax_faf_nc_gnomad()`](https://sigven.github.io/pcgr/dev/reference/grpmax_faf_nc_gnomad.md)
  : Function that assigns a maximum value to a variable
  (gnomAD_NC_FAF_GRPMAX) reflecting the filter allele frequency (GrpMAX)
  for a given variant in the non-cancer gnomAD subset (v3.1)
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
- [`init_biomarker_content()`](https://sigven.github.io/pcgr/dev/reference/init_biomarker_content.md)
  : Function that initiates report element with biomarker evidence
  information
- [`init_expression_content()`](https://sigven.github.io/pcgr/dev/reference/init_expression_content.md)
  : Function that initiates report element with expression information
- [`init_fusion_content()`](https://sigven.github.io/pcgr/dev/reference/init_fusion_content.md)
  : Function that initiates report element with fusion information
- [`init_germline_content()`](https://sigven.github.io/pcgr/dev/reference/init_germline_content.md)
  : Function that initiates report element with germline variant
  information (CPSR)
- [`init_kataegis_content()`](https://sigven.github.io/pcgr/dev/reference/init_kataegis_content.md)
  : Function that initiates report element with kataegis information
- [`init_msi_content()`](https://sigven.github.io/pcgr/dev/reference/init_msi_content.md)
  : Function that initiates report element with MSI classification
- [`init_mutsignature_content()`](https://sigven.github.io/pcgr/dev/reference/init_mutsignature_content.md)
  : Function that initiates report element with mutational signatures
  information
- [`init_rainfall_content()`](https://sigven.github.io/pcgr/dev/reference/init_rainfall_content.md)
  : Function that initiates report element with rainfall information
- [`init_report()`](https://sigven.github.io/pcgr/dev/reference/init_report.md)
  : Function that initiates PCGR/CPSR report object
- [`init_tmb_content()`](https://sigven.github.io/pcgr/dev/reference/init_tmb_content.md)
  : Function that initiates report element with TMB information
- [`init_tumor_only_content()`](https://sigven.github.io/pcgr/dev/reference/init_tumor_only_content.md)
  : Function that initiates report element with tumor-only information
- [`init_var_content()`](https://sigven.github.io/pcgr/dev/reference/init_var_content.md)
  : Function that initiates report element with variant data
- [`init_vstats_actionable()`](https://sigven.github.io/pcgr/dev/reference/init_vstats_actionable.md)
  : Function that initiates report element with actionable variant
  statistics information
- [`init_vstats_cna()`](https://sigven.github.io/pcgr/dev/reference/init_vstats_cna.md)
  : Function that initiates CNA statistics
- [`init_vstats_fusion()`](https://sigven.github.io/pcgr/dev/reference/init_vstats_fusion.md)
  : Function that initiates RNA fusion statistics
- [`init_vstats_snv_indel()`](https://sigven.github.io/pcgr/dev/reference/init_vstats_snv_indel.md)
  : Function that initiates report element with SNV/InDel statistics
  information
- [`kataegis_detect()`](https://sigven.github.io/pcgr/dev/reference/kataegis_detect.md)
  : Function that detects kataegis events from a data frame with genomic
  cooordinates of mutations
- [`kataegis_input()`](https://sigven.github.io/pcgr/dev/reference/kataegis_input.md)
  : Function that detects kataegis events from a data frame with genomic
  cooordinates of mutations
- [`list_of_list_to_df()`](https://sigven.github.io/pcgr/dev/reference/list_of_list_to_df.md)
  : Helper function to convert list of lists to data.frame
- [`load_cpsr_classified_variants()`](https://sigven.github.io/pcgr/dev/reference/load_cpsr_classified_variants.md)
  : Function that reads CPSR-classified variants from a TSV file
- [`load_dna_variants()`](https://sigven.github.io/pcgr/dev/reference/load_dna_variants.md)
  : Function that reads and validates CNA or SNV/InDel TSV files file
  from PCGR/CPSR pre-report (Python) pipeline
- [`load_expression_csq()`](https://sigven.github.io/pcgr/dev/reference/load_expression_csq.md)
  : Load expression consequence settings
- [`load_expression_outliers()`](https://sigven.github.io/pcgr/dev/reference/load_expression_outliers.md)
  : Load expression outlier results
- [`load_expression_similarity()`](https://sigven.github.io/pcgr/dev/reference/load_expression_similarity.md)
  : Load expression similarity results
- [`load_reference_data()`](https://sigven.github.io/pcgr/dev/reference/load_reference_data.md)
  : Function that parses and loads reference data from files in the
  assembly-specific PCGR bundle directory
- [`load_rna_fusions()`](https://sigven.github.io/pcgr/dev/reference/load_rna_fusions.md)
  : Load RNA fusion results
- [`load_somatic_cna()`](https://sigven.github.io/pcgr/dev/reference/load_somatic_cna.md)
  : Function that reads and validates fully annotated CNA data (segments
  and genes) from PCGR pre-reporting pipeline
- [`load_somatic_snv_indel()`](https://sigven.github.io/pcgr/dev/reference/load_somatic_snv_indel.md)
  : Function that reads and validates an annotated somatic SNV/InDel
  file from PCGR pre-reporting pipeline
- [`load_yaml()`](https://sigven.github.io/pcgr/dev/reference/load_yaml.md)
  : Function that loads YAML data with settings and file paths to
  annotated molecular profiles
- [`lof_doc_note()`](https://sigven.github.io/pcgr/dev/reference/lof_doc_note.md)
  : Get documentation string for loss-of-function annotation
- [`log4r_debug()`](https://sigven.github.io/pcgr/dev/reference/log4r_debug.md)
  : Write messages to logs at a given priority level
- [`log4r_fatal()`](https://sigven.github.io/pcgr/dev/reference/log4r_fatal.md)
  : Write messages to logs at a given priority level
- [`log4r_info()`](https://sigven.github.io/pcgr/dev/reference/log4r_info.md)
  : Write messages to logs at a given priority level
- [`log4r_warn()`](https://sigven.github.io/pcgr/dev/reference/log4r_warn.md)
  : Write messages to logs at a given priority level
- [`map_biomarker_data()`](https://sigven.github.io/pcgr/dev/reference/map_biomarker_data.md)
  : Function that maps biomarker identifiers from VCF variant annotation
  to full biomarker data tables for display in report
- [`max_af_gnomad()`](https://sigven.github.io/pcgr/dev/reference/max_af_gnomad.md)
  : Function that assigns a maximum value to a variable (MAX_AF_GNOMAD)
  reflecting the maximum allele frequency for a given variant across
  gnomAD populations
- [`mkdir()`](https://sigven.github.io/pcgr/dev/reference/mkdir.md) :
  Create directory
- [`msi_doc_note()`](https://sigven.github.io/pcgr/dev/reference/msi_doc_note.md)
  : Get documentation string for MSI status prediction
- [`msi_indel_fraction_plot()`](https://sigven.github.io/pcgr/dev/reference/msi_indel_fraction_plot.md)
  : Function that plots the indel fraction for a given sample and
  contrasts this with the distribution for MSI-H/MSS samples from TCGA
- [`msi_indel_load_plot()`](https://sigven.github.io/pcgr/dev/reference/msi_indel_load_plot.md)
  : Function that plots the indel load for a given sample and contrasts
  this with the distribution for MSI-H/MSS samples from TCGA
- [`mutational_signatures_doc_note()`](https://sigven.github.io/pcgr/dev/reference/mutational_signatures_doc_note.md)
  : Get documentation string for mutational signatures analysis
- [`oncogenicity_criteria`](https://sigven.github.io/pcgr/dev/reference/oncogenicity_criteria.md)
  : Oncogenicity criteria (ClinGen/CGC/VICC)
- [`oncogenicity_doc_note()`](https://sigven.github.io/pcgr/dev/reference/oncogenicity_doc_note.md)
  : Get documentation string for oncogenicity annotation
- [`oncokb_annotations`](https://sigven.github.io/pcgr/dev/reference/oncokb_annotations.md)
  : Character vector with OncoKB annotations coming from the
  MafAnnotator / FusionAnnotator / CnaAnnnotator tools in the PCGR
  Python workflow. These annotations are used for variant classification
  and reporting in PCGR.
- [`oncokb_base_api_url`](https://sigven.github.io/pcgr/dev/reference/oncokb_base_api_url.md)
  : OncoKB base API URL
- [`order_variants()`](https://sigven.github.io/pcgr/dev/reference/order_variants.md)
  : Function that orders genomic aberrations according to order of
  chromosomes and chromosomal position
- [`pkg_exists()`](https://sigven.github.io/pcgr/dev/reference/pkg_exists.md)
  : Does R Package Exist
- [`plot_cna_segments_absolute()`](https://sigven.github.io/pcgr/dev/reference/plot_cna_segments_absolute.md)
  : Plot allele-specific copy number segments (absolute copies)
- [`plot_cna_segments_relative()`](https://sigven.github.io/pcgr/dev/reference/plot_cna_segments_relative.md)
  : Plot copy number segments (relative log2 fold change)
- [`plot_filtering_stats_exonic()`](https://sigven.github.io/pcgr/dev/reference/plot_filtering_stats_exonic.md)
  : Function that generates a pie chart for exonic/non-exonic variant
  statistics (for callsets coming from tumor-only sequencing)
- [`plot_filtering_stats_germline()`](https://sigven.github.io/pcgr/dev/reference/plot_filtering_stats_germline.md)
  : Function that generates a pie chart for germline filtering
  statistics for callsets coming from tumor-only sequencing
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
- [`popmax_af_gnomad()`](https://sigven.github.io/pcgr/dev/reference/popmax_af_gnomad.md)
  : Function that assigns a maximum value to a variable
  (gnomAD_AF_POPMAX) reflecting the maximum allele frequency for a given
  variant across gnomAD populations
- [`predict_msi_status()`](https://sigven.github.io/pcgr/dev/reference/predict_msi_status.md)
  : Function that predicts MSI status based on fraction of indels among
  calls
- [`prep_actble_display_tbl()`](https://sigven.github.io/pcgr/dev/reference/prep_actble_display_tbl.md)
  : Function that gathers data tables on actionable variants for display
  in report (tier 1 + tier 2)
- [`prep_diagn_display_tbl()`](https://sigven.github.io/pcgr/dev/reference/prep_diagn_display_tbl.md)
  : Function that gathers data tables on diagnostic variants
- [`prep_progn_display_tbl()`](https://sigven.github.io/pcgr/dev/reference/prep_progn_display_tbl.md)
  : Function that gathers data tables on prognostic variants for display
  in report
- [`process_oncokb_cna()`](https://sigven.github.io/pcgr/dev/reference/process_oncokb_cna.md)
  : Process OncoKB CNA output file and fetch complete annotations
- [`process_oncokb_fusion()`](https://sigven.github.io/pcgr/dev/reference/process_oncokb_fusion.md)
  : Process OncoKB fusion output file and fetch complete annotations
- [`process_oncokb_maf()`](https://sigven.github.io/pcgr/dev/reference/process_oncokb_maf.md)
  : Process OncoKB MAF output files (both HGVSp and HGVSg) and fetch
  complete annotations
- [`remove_cols_from_df()`](https://sigven.github.io/pcgr/dev/reference/remove_cols_from_df.md)
  : Function that removes column(s) from data frame
- [`render_actble_bm_table()`](https://sigven.github.io/pcgr/dev/reference/render_actble_bm_table.md)
  : Build biomarker reactable with category-aware styling Combines tier
  1 and tier 2 records in one table. Header uses tier 1 color;
  THERAPY_MATCH cell background reflects the row's tier (1 or 2).
- [`render_alteration_cell()`](https://sigven.github.io/pcgr/dev/reference/render_alteration_cell.md)
  : Render sample alteration with confidence icon Uses filled circle +
  bold for high confidence, hollow circle + muted for medium confidence
- [`render_alteration_clinvar_style()`](https://sigven.github.io/pcgr/dev/reference/render_alteration_clinvar_style.md)
  : Style alteration cell background based on ClinVar classification
  (Generated by Claude Opus 4.6 with some manual tweaks)
- [`render_bar_cell()`](https://sigven.github.io/pcgr/dev/reference/render_bar_cell.md)
  : Render a value between 0 and 1 as a horizontal bar (Generated by
  Claude Opus 4.6 with some manual tweaks)
- [`render_diagn_bm_table()`](https://sigven.github.io/pcgr/dev/reference/render_diagn_bm_table.md)
  : Build biomarker reactable with category-aware styling Combines tier
  1 and tier 2 records in one table. Header uses tier 1 color;
  THERAPY_MATCH cell background reflects the row's tier (1 or 2).
- [`render_diagnosis()`](https://sigven.github.io/pcgr/dev/reference/render_diagnosis.md)
  : Render diagnosis cell with background color based on actionability
  tier (Generated by Claude Opus 4.6 with some manual tweaks)
- [`render_evidence_desc_cell()`](https://sigven.github.io/pcgr/dev/reference/render_evidence_desc_cell.md)
  : Render evidence description with truncation and tooltip (Generated
  by Claude Opus 4.6 with some manual tweaks)
- [`render_evidence_level_cell()`](https://sigven.github.io/pcgr/dev/reference/render_evidence_level_cell.md)
  : Render evidence level badge (Generated by Claude Opus 4.6 with some
  manual tweaks)
- [`render_oncogenicity_cell()`](https://sigven.github.io/pcgr/dev/reference/render_oncogenicity_cell.md)
  : Style oncogenicity cell background based on value (Generated by
  Claude Opus 4.6 with some manual tweaks)
- [`render_progn_bm_table()`](https://sigven.github.io/pcgr/dev/reference/render_progn_bm_table.md)
  : Build biomarker reactable with category-aware styling Combines tier
  1 and tier 2 records in one table. Header uses tier 1 color;
  THERAPY_MATCH cell background reflects the row's tier (1 or 2).
- [`render_prognostic_outcome()`](https://sigven.github.io/pcgr/dev/reference/render_prognostic_outcome.md)
  : Render prognostic outcome cell with background color based on
  actionability tier (Generated by Claude Opus 4.6 with some manual
  tweaks)
- [`render_source_logos()`](https://sigven.github.io/pcgr/dev/reference/render_source_logos.md)
  : Render one or more source logos from a pipe-separated string
  (Generated by Claude Opus 4.6 with some manual tweaks)
- [`render_symbol_assoc_style()`](https://sigven.github.io/pcgr/dev/reference/render_symbol_assoc_style.md)
  : Style alteration cell background based on cancer association rank
  (Generated by Claude Opus 4.6 with some manual tweaks)
- [`render_therapy_style()`](https://sigven.github.io/pcgr/dev/reference/render_therapy_style.md)
  : Render therapy match cell with background color based on
  actionability tier (Generated by Claude Opus 4.6 with some manual
  tweaks)
- [`render_tier_cell()`](https://sigven.github.io/pcgr/dev/reference/render_tier_cell.md)
  : Render actionability tier as a styled badge with tier-specific
  colors
- [`rna_fusion_recurrence_mitdb()`](https://sigven.github.io/pcgr/dev/reference/rna_fusion_recurrence_mitdb.md)
  : Function that checks for recurrence of RNA fusions in Mitelman
  database and annotates with clinical and cytogenetic features
- [`rt_theme()`](https://sigven.github.io/pcgr/dev/reference/rt_theme.md)
  : Predefined reactable theme for PCGR/CPSR tables (Generated by Claude
  Opus 4.6 with some manual tweaks)
- [`sort_chromosomal_segments()`](https://sigven.github.io/pcgr/dev/reference/sort_chromosomal_segments.md)
  : Function that sorts chromosomal segments according to chromosome and
  chromosomal start/end position
- [`stats_report_cna()`](https://sigven.github.io/pcgr/dev/reference/stats_report_cna.md)
  : Function that computes various variant statistics for CNAs from a
  callset object
- [`stats_report_fusion()`](https://sigven.github.io/pcgr/dev/reference/stats_report_fusion.md)
  : Function that computes various variant statistics for fusions from a
  callset object
- [`stats_report_germline()`](https://sigven.github.io/pcgr/dev/reference/stats_report_germline.md)
  : Function that generate stats for a germline variant callset,
  including number of variants and number of variants with evidence
  items for each BM evidence type
- [`stats_report_snv_indel()`](https://sigven.github.io/pcgr/dev/reference/stats_report_snv_indel.md)
  : Function that computes various variant statistics for SNVs/InDels
  from a callset object
- [`stats_type_snv_indel()`](https://sigven.github.io/pcgr/dev/reference/stats_type_snv_indel.md)
  : Function that computes various variant statistics from a data frame
  with variant records
- [`strip_html()`](https://sigven.github.io/pcgr/dev/reference/strip_html.md)
  : Strip HTML tags
- [`table_display_cols`](https://sigven.github.io/pcgr/dev/reference/table_display_cols.md)
  : HTML Table display columns
- [`tcga_cohorts`](https://sigven.github.io/pcgr/dev/reference/tcga_cohorts.md)
  : Data frame with all TCGA cohorts
- [`tcga_somatic_status()`](https://sigven.github.io/pcgr/dev/reference/tcga_somatic_status.md)
  : Function that assigns a logical (STATUS_TCGA_SOMATIC) reflecting
  whether a variant co-incides with an entry in TCGA (somatic)
- [`tier_af_distribution()`](https://sigven.github.io/pcgr/dev/reference/tier_af_distribution.md)
  : Function that plots a histogram of the the variant allelic support
  (tumor) - grouped by tiers
- [`tmb_doc_note()`](https://sigven.github.io/pcgr/dev/reference/tmb_doc_note.md)
  : Get documentation string for tumor mutational burden (TMB)
  annotation
- [`tsv_cols`](https://sigven.github.io/pcgr/dev/reference/tsv_cols.md)
  : TSV columns
- [`tumor_sites`](https://sigven.github.io/pcgr/dev/reference/tumor_sites.md)
  : Character vector of tumor site names used in PCGR/CPSR reports
- [`twohit_doc_note()`](https://sigven.github.io/pcgr/dev/reference/twohit_doc_note.md)
  : Get documentation string for two-hit detection & VAF consistency
  assessment
- [`update_report()`](https://sigven.github.io/pcgr/dev/reference/update_report.md)
  : Function that updates a PCGR/CPSR report object structure
- [`vaf_plot()`](https://sigven.github.io/pcgr/dev/reference/vaf_plot.md)
  : Function that generates a VAF distribution plot for a given PCGR
  report object
- [`variant_db_url`](https://sigven.github.io/pcgr/dev/reference/variant_db_url.md)
  : List of URLS and variant identifiers for variant/gene/protein domain
  databases
- [`write_processed_vcf()`](https://sigven.github.io/pcgr/dev/reference/write_processed_vcf.md)
  : Function that writes a VCF intended for mutational signature
  analysis
- [`write_report_excel()`](https://sigven.github.io/pcgr/dev/reference/write_report_excel.md)
  : Function that writes key datasets of PCGR object to an Excel
  workbook
- [`write_report_html()`](https://sigven.github.io/pcgr/dev/reference/write_report_html.md)
  : Function that writes contents of PCGR object to an HTML report
  (quarto-based)
- [`write_report_tsv()`](https://sigven.github.io/pcgr/dev/reference/write_report_tsv.md)
  : Function that writes contents of PCGR object to various output
  formats (Rmarkdown/flexdashboard HTML reports, JSON, tab-separated
  etc)
