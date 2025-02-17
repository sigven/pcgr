#!/usr/bin/env python

from pcgr import pcgr_vars
from pcgr.utils import check_file_exists, error_message

import pandas as pd
import os
import gzip
import csv
import re

def create_config(arg_dict, workflow = "PCGR"):
    
    conf_options = {}
    if workflow == "PCGR" or workflow == "CPSR":
        conf_options = {
            'sample_id': arg_dict['sample_id'],
            'genome_assembly': arg_dict['genome_assembly'],
            'debug': arg_dict['debug'],
            'pcgrr_conda': arg_dict['pcgrr_conda'],
            'output_prefix': "None"                
    	}
        
        conf_options['vep'] = {
            'vep_buffer_size': int(arg_dict['vep_buffer_size']),
            'vep_pick_order': str(arg_dict['vep_pick_order']),
            'vep_n_forks': int(arg_dict['vep_n_forks']),
            'vep_no_intergenic': int(arg_dict['vep_no_intergenic']),
            'vep_regulatory': 1,
            'vep_gencode_basic': int(arg_dict['vep_gencode_basic'])            
        }
            
        conf_options['other'] = {
            'vcfanno_n_proc': int(arg_dict['vcfanno_n_proc']),                                          
            'no_reporting': int(arg_dict['no_reporting']),
            'no_html': int(arg_dict['no_html']),
            'retained_vcf_info_tags': str(arg_dict['retained_info_tags']),
            'show_noncoding': not int(arg_dict['ignore_noncoding']),
            'force_overwrite': int(arg_dict['force_overwrite'])
        }
        conf_options['sample_properties'] = {}
        conf_options['molecular_data'] = {}
        conf_options['molecular_data']['fname_mut_vcf'] = "None"
        conf_options['molecular_data']['fname_mut_tsv'] = "None"
        
    if workflow == 'PCGR':
        conf_options['assay_properties'] = {}
        conf_options['sample_properties']['tumor_purity'] = 'NA'
        conf_options['sample_properties']['tumor_ploidy'] = 'NA'
        conf_options['sample_properties']['site'] = str(pcgr_vars.tsites[arg_dict['tsite']])
        conf_options['sample_properties']['site2'] = str(pcgr_vars.tsites[arg_dict['tsite']]).replace(" ","_").replace("/","@")
        conf_options['sample_properties']['dp_control_detected'] = 0
        conf_options['sample_properties']['vaf_control_detected'] = 0
        conf_options['sample_properties']['dp_tumor_detected'] = 0
        conf_options['sample_properties']['vaf_tumor_detected'] = 0
        conf_options['assay_properties']['type'] = str(arg_dict['assay'])
        conf_options['assay_properties']['vcf_tumor_only'] = 0
        conf_options['assay_properties']['mode'] = "Tumor-Control"
        conf_options['assay_properties']['effective_target_size_mb'] = float(arg_dict['effective_target_size_mb'])
        if int(arg_dict['tumor_only']) == 1:
            conf_options['assay_properties']['vcf_tumor_only'] = 1
            conf_options['assay_properties']['mode'] = "Tumor-Only"
        
        if not arg_dict['tumor_purity'] is None:
            conf_options['sample_properties']['tumor_purity'] = float(arg_dict['tumor_purity'])
        if not arg_dict['tumor_ploidy'] is None:
            conf_options['sample_properties']['tumor_ploidy'] = float(arg_dict['tumor_ploidy'])
            
        #conf_options['clinicaltrials'] = {
        #    'run': int(arg_dict['include_trials'])
        #}
        conf_options['other']['vcf2maf'] = int(arg_dict['vcf2maf'])
        
        conf_options['somatic_cna'] = {            
            'cna_overlap_pct': float(arg_dict['cna_overlap_pct']),
            'n_copy_gain': int(arg_dict['n_copy_gain'])
        }
        
        conf_options['germline'] = {
            'show': 0,
            'ignore_vus': int(arg_dict['cpsr_ignore_vus'])
        }
        if not arg_dict['input_cpsr'] is None:
            conf_options['germline']['show'] = 1
            
        conf_options['expression'] = {}
        conf_options['expression']['run'] = int(not arg_dict['input_rna_exp'] is None)
        conf_options['expression']['similarity_analysis'] = int(arg_dict['expression_sim'])
        conf_options['expression']['similarity_db'] = {}
        for db in arg_dict['expression_sim_db'].split(','):
            conf_options['expression']['similarity_db'][db] = 1
            if db == 'tcga':
                conf_options['expression']['similarity_db']['tcga'] = {}
                for cohort in pcgr_vars.DISEASE_COHORTS:
                    conf_options['expression']['similarity_db']['tcga'][cohort] = 1
                    
        conf_options['somatic_snv'] = {}
        conf_options['somatic_snv']['allelic_support'] = {
            'tumor_dp_min': int(arg_dict['tumor_dp_min']),
            'control_dp_min': int(arg_dict['control_dp_min']),
            'tumor_af_min': float(arg_dict['tumor_af_min']),
            'control_af_max': float(arg_dict['control_af_max']),
            'control_dp_tag': str(arg_dict['control_dp_tag']),
            'control_af_tag': str(arg_dict['control_af_tag']),
            'tumor_dp_tag': str(arg_dict['tumor_dp_tag']),
            'tumor_af_tag': str(arg_dict['tumor_af_tag']),
            'call_conf_tag': str(arg_dict['call_conf_tag'])
        }

        conf_options['somatic_snv']['tumor_only'] = {
            'exclude_pon': int(arg_dict['exclude_pon']),
            'exclude_likely_hom_germline': int(arg_dict['exclude_likely_hom_germline']),
            'exclude_likely_het_germline': int(arg_dict['exclude_likely_het_germline']),
            'exclude_clinvar_germline': int(arg_dict['exclude_clinvar_germline']),
            'exclude_dbsnp_nonsomatic': int(arg_dict['exclude_dbsnp_nonsomatic']),
            'exclude_nonexonic': int(arg_dict['exclude_nonexonic'])
        }
        conf_options['somatic_snv']['msi'] = {
            'run': int(arg_dict['estimate_msi'])
        }
        conf_options['somatic_snv']['tmb'] = {
            'run': int(arg_dict['estimate_tmb']), 
            'tmb_display': arg_dict['tmb_display'],           
            'tmb_dp_min': arg_dict['tmb_dp_min'],
            'tmb_af_min': arg_dict['tmb_af_min']
        }

        for pop in ['nfe', 'fin', 'amr', 'eas', 'sas', 'asj', 'oth', 'afr', 'global']:
            tag = 'maf_gnomad_' + str(pop)
            if arg_dict[tag]:
                conf_options['somatic_snv']['tumor_only'][tag] = float(arg_dict[tag])
        
        conf_options['somatic_snv']['mutational_signatures'] = {
            'run': int(arg_dict['estimate_signatures']),
            'mutation_limit': int(arg_dict['min_mutations_signatures']),
            'all_reference_signatures': int(arg_dict['all_reference_signatures']),
            'include_artefact_signatures': int(arg_dict['include_artefact_signatures']),
            'prevalence_reference_signatures': float(arg_dict['prevalence_reference_signatures'])
        }
        
        conf_options['molecular_data']['fname_cna_gene_tsv'] = "None"
        conf_options['molecular_data']['fname_cna_segment_tsv'] = "None"
        conf_options['molecular_data']['fname_expression_tsv'] = "None"
        conf_options['molecular_data']['fname_expression_outliers_tsv'] = "None"
        conf_options['molecular_data']['fname_maf_tsv'] = "None"
        conf_options['molecular_data']['fname_germline_tsv'] = "None"
        conf_options['molecular_data']['fname_germline_yaml'] = "None"
        #conf_options['molecular_data']['fname_expression_csq_tsv'] = "None"
        conf_options['molecular_data']['fname_expression_similarity_tsv'] = "None"
        conf_options['molecular_data']['fname_tmb_tsv'] = "None"
    
    
    if workflow == "CPSR":        
        conf_options['sample_properties']['phenotype'] = 'None'
        conf_options['sample_properties']['site'] = 'Hereditary (blood)'
        conf_options['sample_properties']['gt_detected'] = 0
        conf_options['sample_properties']['dp_detected'] = 0
        
        #conf_options['visual_reporting']['table_display'] = str(arg_dict['report_table_display'])
        conf_options['gene_panel'] = {
            'panel_id': str(arg_dict['virtual_panel_id']),
            'description': 'Exploratory virtual gene panel (panel 0)',
            'description_trait': 'None',
            'url': 'None',
            'custom_list_tsv': str(arg_dict['custom_list']),
            'custom_list_name': str(arg_dict['custom_list_name']),
            'custom_list_bed': 'None',
            'diagnostic_grade_only': int(arg_dict['diagnostic_grade_only'])
        }
            
        conf_options['variant_classification'] = {
            'gwas_findings': int(arg_dict['gwas_findings']),
            'secondary_findings': int(arg_dict['secondary_findings']),
            'pgx_findings': int(arg_dict['pgx_findings']),
            'pop_gnomad': str(arg_dict['pop_gnomad']),
            'maf_upper_threshold': float(arg_dict['maf_upper_threshold']),
            'classify_all': int(arg_dict['classify_all']),
            'clinvar_report_noncancer': int(arg_dict['clinvar_report_noncancer'])
        }
            

    return conf_options


def populate_config_data(conf_options: dict, refdata_assembly_dir: str, workflow = "PCGR", logger=None):
    
    conf_data = {}
    
    ## move/set some parameters to highest level in YAML/conf_data object
    conf_data['sample_id'] = conf_options['sample_id']
    conf_data['output_dir'] = conf_options['output_dir']
    conf_data['output_prefix'] = conf_options['output_prefix']
    conf_data['workflow'] = workflow
    conf_data['genome_assembly'] = conf_options['genome_assembly']
    conf_data['software'] = {}
    conf_data['software']['pcgr_version'] = pcgr_vars.PCGR_VERSION
    conf_data['software']['cpsr_version'] = pcgr_vars.PCGR_VERSION
    conf_data['molecular_data'] = conf_options['molecular_data']
    
    ## remove the core parameters set above from 'conf_options', 
    ## (now placed at the highest level in YAML/conf_data object)
    for e in ['sample_id','genome_assembly','output_dir',
              'molecular_data','output_prefix']:
        del conf_options[e]

    conf_data['conf'] = conf_options
    
    conf_data['reference_data'] = {}
    conf_data['reference_data']['version'] = pcgr_vars.DB_VERSION
    conf_data['reference_data']['path'] = \
        refdata_assembly_dir
    
    metadata_pd = pd.DataFrame()
    conf_data['reference_data']['source_metadata'] = {}
    conf_data['conf']['sample_properties']['phenotype'] = {}
    
    ## add metadata information for each data source
    
    cpsr_sources_regex = r'^(gepa|cpg_other|maxwell2016|acmg_sf|dbmts|woods_dnarepair|gerp|tcga_pancan_2018|gwas_catalog|cpic)'
    pcgr_sources_regex = r'^(cytoband|mitelmandb|tcga|nci|intogen|opentargets|depmap|treehouse|dgidb|pubchem|cosmic_mutsigs)'
    sources_skip_regex = r'^(illumina|foundation_one)'
    
    for dtype in ['gene','phenotype','biomarker','drug','gwas','hotspot','other']:
        metadata_fname = os.path.join(
            refdata_assembly_dir,
            ".METADATA", "tsv", dtype + "_metadata.tsv")
        if check_file_exists(metadata_fname, logger):
            metadata_df = pd.read_csv(metadata_fname, sep="\t", na_values=".")
            metadata_df["source_type"] = dtype
            metadata_df['wflow'] = 'pcgr_cpsr'
            metadata_df.loc[metadata_df['source_abbreviation'].str.match(cpsr_sources_regex), 'wflow'] = 'cpsr' 
            metadata_df.loc[metadata_df['source_abbreviation'].str.match(pcgr_sources_regex), 'wflow'] = 'pcgr'
            metadata_df.loc[metadata_df['source_abbreviation'].str.match(sources_skip_regex), 'wflow'] = 'skip'
            metadata_pd = metadata_pd._append(metadata_df, ignore_index=True)
    
    conf_data['reference_data']['source_metadata'] = metadata_pd.to_dict(orient='records')

    oncotree_fname = os.path.join(refdata_assembly_dir,
                                  "phenotype","tsv","phenotype_onco.tsv.gz")
    
    if check_file_exists(oncotree_fname, logger):
        oncotree_df = pd.read_csv(oncotree_fname, sep="\t", na_values=".")
        if workflow == "CPSR":
            conf_data['conf']['sample_properties']['phenotype'] = \
                oncotree_df[oncotree_df['primary_site'].isna()].to_dict(orient='records')
        else:
            if conf_data['conf']['sample_properties']['site'] != "Any":
                oncotree_site = oncotree_df.loc[oncotree_df['primary_site'] == conf_data['conf']['sample_properties']['site'],:].copy()
                conf_data['conf']['sample_properties']['phenotype'] = oncotree_site.to_dict(orient='records')
            else:
                conf_data['conf']['sample_properties']['phenotype'] = 'None'
    
    if workflow == "CPSR":
        if conf_data['conf']['gene_panel']['panel_id'] is not None:
                
            if conf_data['conf']['gene_panel']['panel_id'] == "-1":
                conf_data['conf']['gene_panel']['description'] = 'User-defined panel (custom geneset from panel 0)'
                conf_data['conf']['gene_panel']['description_trait'] = 'User-defined panel (custom geneset from panel 0)'
            else:
                if ',' in conf_data['conf']['gene_panel']['panel_id']:
                    conf_data['conf']['gene_panel']['description'] = \
                        'Genomics England PanelApp - multiple panels (' + conf_data['conf']['gene_panel']['panel_id'] + ')'
                    conf_data['conf']['gene_panel']['description2'] = \
                        'Genomics England PanelApp - multiple panels (' + conf_data['conf']['gene_panel']['panel_id'] + ')'
                else:
                    if conf_data['conf']['gene_panel']['panel_id'] != "0":
                        conf_data['conf']['gene_panel']['description'] = \
                            'Genomics England PanelApp - panel ' + conf_data['conf']['gene_panel']['panel_id']
        
            conf_data['conf']['gene_panel']['panel_genes'] = set_virtual_target_genes(
                conf_data['conf']['gene_panel']['panel_id'], 
                refdata_assembly_dir,
                conf_data['conf']['gene_panel']['diagnostic_grade_only'], 
                conf_data['conf']['gene_panel']['custom_list_tsv'],
                bool(conf_data['conf']['variant_classification']['secondary_findings']),
                bool(conf_data['conf']['variant_classification']['pgx_findings']),
                logger)
            
            if conf_data['conf']['gene_panel']['panel_id'] != "-1":
                if not ',' in conf_data['conf']['gene_panel']['panel_id']:                     
                    conf_data['conf']['gene_panel']['url'] = str(conf_data['conf']['gene_panel']['panel_genes'][0]['panel_url'])
                    conf_data['conf']['gene_panel']['description_trait'] = str(conf_data['conf']['gene_panel']['panel_genes'][0]['panel_name'])
                    
                
        

    return(conf_data)

def set_virtual_target_genes(panel_id: str, refdata_assembly_dir: str, diagnostic_grade_only: bool, 
                             custom_list_tsv: str, secondary_findings: bool, pgx_findings: bool, logger=None):
    
    all_panels_fname = os.path.join(
        refdata_assembly_dir,
        "gene","tsv","gene_virtual_panel", 
        "gene_virtual_panel.tsv.gz")
    
    all_cpg_fname = os.path.join(
        refdata_assembly_dir,
        "gene","tsv","gene_cpg", 
        "gene_cpg.tsv.gz")
        
    all_virtual_panels = pd.DataFrame()
    all_sf_pgx_targets = pd.DataFrame()
    panel_targets = pd.DataFrame()
    if check_file_exists(all_panels_fname, logger):    
        all_virtual_panels = pd.read_csv(all_panels_fname, sep="\t", na_values=".")
        all_virtual_panels['id'] = all_virtual_panels['id'].astype(str)
        all_virtual_panels = \
            all_virtual_panels.drop(['gepa_phenotype', 'genome_build','gepa_penetrance'], axis=1)
        all_virtual_panels = \
            all_virtual_panels.rename(
                columns={'gepa_moi': 'moi',
                         'gepa_mod': 'mod',
                         'gepa_panel_version': 'panel_version', 
                         'gepa_confidence_level': 'confidence_level',
                         'gepa_panel_url': 'panel_url',
                         'gepa_panel_name': 'panel_name',
                         'gepa_panel_id': 'panel_id'})
        all_virtual_panels['primary_target'] = True
        for f in ['panel_url', 'panel_name','ensembl_gene_id','moi','mod','symbol']:
            all_virtual_panels[f] = all_virtual_panels[f].astype(str)
    
    if secondary_findings is True or pgx_findings is True:
        if check_file_exists(all_cpg_fname, logger):
            all_sf_pgx_targets = pd.read_csv(all_cpg_fname, sep="\t", na_values=".")
            regex = r'(ACMG_SF|CPIC_PGX_ONCOLOGY)'
            if secondary_findings is False:
                regex = r'(CPIC_PGX_ONCOLOGY)'
            else:
                if pgx_findings is False:
                    regex = r'(ACMG_SF)'
                
            all_sf_pgx_targets = \
                all_sf_pgx_targets[all_sf_pgx_targets.cpg_source.str.match(regex)]    
            all_sf_pgx_targets = \
                all_sf_pgx_targets.drop(
                    ['cpg_phenotypes', 'cpg_cancer_cui','cpg_syndrome_cui'], axis=1)
            all_sf_pgx_targets = \
                all_sf_pgx_targets.rename(
                    columns={'cpg_moi': 'moi',
                             'cpg_mod': 'mod'}) 
            all_sf_pgx_targets['confidence_level'] = -1
            all_sf_pgx_targets['panel_id'] = float(3.2)
            all_sf_pgx_targets['panel_url'] = "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"
            all_sf_pgx_targets['panel_version'] = float(3.2)
            all_sf_pgx_targets['panel_name'] = all_sf_pgx_targets['cpg_source']
            all_sf_pgx_targets['primary_target'] = False
            all_sf_pgx_targets['id'] = str("-2")
            for f in ['panel_url', 'panel_name','moi','mod','symbol','ensembl_gene_id']:
                all_sf_pgx_targets[f] = all_sf_pgx_targets[f].astype(str)


    if panel_id == "-1":
        if not custom_list_tsv == 'None':
            custom_ensembl_gene_ids = {}
            if check_file_exists(custom_list_tsv, logger):
                custom_genes = csv.DictReader(open(custom_list_tsv,'r'), delimiter='\n', fieldnames=['ensembl_gene_id'])                
                for row in custom_genes:
                    custom_ensembl_gene_ids[row['ensembl_gene_id']] = 1            
            
            panel_targets = all_virtual_panels[all_virtual_panels['id'] == "0"].copy()
            panel_targets = panel_targets[panel_targets['ensembl_gene_id'].isin(custom_ensembl_gene_ids)]
            if len(panel_targets) == 0:
                err_msg = "Custom list of genes from CPSR superpanel (panel 0) should be provided as Ensembl " + \
                    "gene identifiers, '" + str(custom_ensembl_gene_ids) + "' were not find in panel 0"
                error_message(err_msg, logger)
            else:
                panel_targets.loc[:,'confidence_level'] = -1
                panel_targets.loc[:,'panel_id'] = -5
                panel_targets.loc[:,'panel_url'] = 'None'
                panel_targets.loc[:,'panel_version'] = 1.0
                panel_targets.loc[:,'panel_name'] = "CustomPanel"
                panel_targets.loc[:,'primary_target'] = True
                for f in ['panel_url', 'panel_name','moi','mod','symbol','ensembl_gene_id']:
                    panel_targets[f] = panel_targets[f].astype(str)
            
    else:
        panel_ids = panel_id.split(',')            
        panel_targets = all_virtual_panels[all_virtual_panels['id'].isin(panel_ids)].copy()
        if diagnostic_grade_only:
            panel_targets = panel_targets[panel_targets['confidence_level'] >= 3]
            if len(panel_targets) == 0:
                err_msg = \
                    f"No set of target genes in PanelApp panel(s) '{panel_ids}' that are diagnostic-grade (green level) - consider skipping 'diagnostic_grade_only'"
                error_message(err_msg, logger)
            
        if len(panel_ids) > 1:
           panel_targets.loc[:,'confidence_level'] = 5  

    all_targets = panel_targets
    
    if len(all_sf_pgx_targets) > 0:
        df1 = pd.DataFrame(all_sf_pgx_targets['ensembl_gene_id'])
        df2 = pd.DataFrame(panel_targets['ensembl_gene_id'])
        outer = df1.merge(df2, on = "ensembl_gene_id", how='outer', indicator=True) ## find genes in SF-PGX that are not in panel
        df1_only = outer[(outer._merge=='left_only')].drop('_merge', axis=1)
        all_sf_pgx_targets = all_sf_pgx_targets[all_sf_pgx_targets['ensembl_gene_id'].isin(df1_only['ensembl_gene_id'])]        
        all_targets = pd.concat([panel_targets, all_sf_pgx_targets], axis=0)
        
    
    return all_targets.to_dict(orient='records')
    
