#!/usr/bin/env python

from pcgr import pcgr_vars
from pcgr.utils import check_file_exists, error_message, warn_message

import pandas as pd
import os
import csv
import requests

def fetch_oncokb_info(api_token: str, logger=None) -> dict:
    """Fetch OncoKB data version, release date, and API version from the OncoKB info endpoint."""
    result = {
        'data_version': None,
        'data_release_date': None,
        'api_version': None
    }
    try:
        response = requests.get(
            "https://www.oncokb.org/api/v1/info",
            headers={"Authorization": f"Bearer {api_token}"},
            timeout=10
        )
        if response.status_code == 200:
            info = response.json()
            result['data_version'] = info.get('dataVersion', {}).get('version')
            result['data_release_date'] = info.get('dataVersion', {}).get('date')
            result['api_version'] = info.get('apiVersion', {}).get('version')
        else:
            warn_message(
                f"OncoKB info endpoint returned HTTP {response.status_code} - "
                "data_version/api_version will not be recorded", logger)
    except Exception as e:
        warn_message(f"Could not retrieve OncoKB version info: {e}", logger)
    return result


def create_config(arg_dict, workflow = "PCGR", logger=None):
    
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
        conf_options['sample_properties']['tumor_purity'] = 'NA' if \
            arg_dict['tumor_purity'] is None else float(arg_dict['tumor_purity'])
        conf_options['sample_properties']['tumor_ploidy'] = 'NA' if \
            arg_dict['tumor_ploidy'] is None else float(arg_dict['tumor_ploidy'])
        conf_options['sample_properties']['tumor_ploidy_source'] = 'NA'
        conf_options['sample_properties']['sex'] = 'UNKNOWN'
        conf_options['sample_properties']['site'] = str(pcgr_vars.tsites[arg_dict['tsite']])
        conf_options['sample_properties']['site2'] = str(pcgr_vars.tsites[arg_dict['tsite']]).replace(" ","_").replace("/","@")
        for prop in ['dp_control_detected', 'vaf_control_detected', 'dp_tumor_detected', 'vaf_tumor_detected']:
            conf_options['sample_properties'][prop] = 0
        conf_options['assay_properties']['type'] = str(arg_dict['assay'])
        conf_options['assay_properties']['vcf_tumor_only'] = 0
        conf_options['assay_properties']['mode'] = "Tumor-Control"
        conf_options['assay_properties']['effective_target_size_mb'] = float(arg_dict['effective_target_size_mb'])
        if int(arg_dict['tumor_only']) == 1:
            conf_options['assay_properties']['vcf_tumor_only'] = 1
            conf_options['assay_properties']['mode'] = "Tumor-Only"
                
        if arg_dict['sex'] is not None:
            sex = str(arg_dict['sex'])
            ## Breast, Ovary/Fallopian Tube, Uterus, Vulva/Vagina
            if arg_dict['tsite'] == 6 or arg_dict['tsite'] == 18 or arg_dict['tsite'] == 29 or arg_dict['tsite'] == 30:
                if sex == 'MALE':
                    if arg_dict['tsite'] == 6:
                        warn_message(f"Tumor site '{conf_options['sample_properties']['site']}' is typically observed in females - ensure '--sex' option is correctly set", logger)
                    else:
                        error_message(f"Tumor site '{conf_options['sample_properties']['site']}' is not observed in males - please check the '--sex' option", logger)
            ## Prostate, Testis
            if arg_dict['tsite'] == 23 or arg_dict['tsite'] == 26:
                if sex == 'FEMALE':
                    error_message(f"Tumor site '{conf_options['sample_properties']['site']}' is not observed in females - please check the '--sex' option", logger)

            conf_options['sample_properties']['sex'] = str(arg_dict['sex'])

            
        #conf_options['clinicaltrials'] = {
        #    'run': int(arg_dict['include_trials'])
        #}
        conf_options['other']['vcf2maf'] = int(arg_dict['vcf2maf'])

        conf_options['oncokb'] = {
            'api_token': str(arg_dict['oncokb_api_token']) if arg_dict['oncokb_api_token'] is not None else None,
            'oncotree_code': str(arg_dict['oncokb_oncotree_code']) if arg_dict['oncokb_oncotree_code'] is not None else None,
            'exclusive': int(arg_dict['oncokb_exclusive']),
            'data_version': None,
            'data_release_date': None,
            'api_version': None
        }
        if conf_options['oncokb']['api_token'] is not None:
            oncokb_info = fetch_oncokb_info(conf_options['oncokb']['api_token'], logger)
            conf_options['oncokb'].update(oncokb_info)
            conf_options['oncokb']['run'] = 1        

        conf_options['somatic_cna'] = {            
            'cna_transcript_overlap_pct': float(arg_dict['cna_transcript_overlap_pct']),
            'threshold_mode': str(arg_dict['cna_threshold_mode']),
            'amp_threshold_absolute': int(arg_dict['cna_amp_threshold_absolute']),
            'amp_threshold_relative': float(arg_dict['cna_amp_threshold_relative']),
            'amp_threshold_effective': float(arg_dict['cna_amp_threshold_absolute']) if \
                arg_dict['cna_threshold_mode'] == 'absolute' else \
                    float(arg_dict['cna_amp_threshold_relative']) * 2,
            'gain_threshold_absolute': int(arg_dict['cna_gain_threshold_absolute']),
            'gain_threshold_relative': float(arg_dict['cna_gain_threshold_relative']),
            'del_threshold_absolute': int(arg_dict['cna_del_threshold_absolute']),
            'del_threshold_relative': float(arg_dict['cna_del_threshold_relative']),
        }
        
        conf_options['germline'] = {
            'show': 0,
            'ignore_vus': int(arg_dict['cpsr_ignore_vus'])
        }
        if arg_dict['input_cpsr'] is not None:
            conf_options['germline']['show'] = 1
            
        conf_options['expression'] = {}
        conf_options['expression']['run'] = int(arg_dict['input_rna_exp'] is not None)
        conf_options['expression']['similarity_analysis'] = int(arg_dict['expression_sim'])
        conf_options['expression']['similarity_db'] = {}
        for db in arg_dict['expression_sim_db'].split(','):
            conf_options['expression']['similarity_db'][db] = 1
            if db == 'tcga':
                conf_options['expression']['similarity_db']['tcga'] = {}
                for cohort in pcgr_vars.DISEASE_COHORTS:
                    conf_options['expression']['similarity_db']['tcga'][cohort] = 1
                    
        conf_options['rna_fusion'] = {
            'min_split_reads': int(arg_dict['fusion_min_split_reads'])
        }

        conf_options['somatic_snv'] = {}
        conf_options['somatic_snv']['allelic_support'] = {
            'tumor_dp_min': arg_dict['tumor_dp_min'],
            'control_dp_min': arg_dict['control_dp_min'],
            'tumor_af_min': arg_dict['tumor_af_min'],
            'control_af_max': arg_dict['control_af_max'],
            'tumor_ad_min': arg_dict['tumor_ad_min'],
            'control_ad_max': arg_dict['control_ad_max'],
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
            'exclude_nonexonic': int(arg_dict['exclude_nonexonic']),
            'gnomad_popmax_af_tolerated': float(arg_dict['gnomad_popmax_af_tolerated'])
        }
        conf_options['somatic_snv']['msi'] = {
            'run': int(arg_dict['estimate_msi'])
        }

        tmb_dp_min = 0 if arg_dict['tmb_dp_min'] is None else int(arg_dict['tmb_dp_min'])
        tmb_af_min = 0.0 if arg_dict['tmb_af_min'] is None else float(arg_dict['tmb_af_min'])
        tmb_ad_min = 0 if arg_dict['tmb_ad_min'] is None else int(arg_dict['tmb_ad_min'])
        conf_options['somatic_snv']['tmb'] = {
            'run': int(arg_dict['estimate_tmb']), 
            'tmb_display': arg_dict['tmb_display'],           
            'tmb_dp_min': tmb_dp_min,
            'tmb_af_min': tmb_af_min,
            'tmb_ad_min': tmb_ad_min
        }

        #conf_options['somatic_snv']['tumor_only']['popmax_af_gnomad'] = float(arg_dict['popmax_af_gnomad'])
        # for pop in ['nfe', 'fin', 'amr', 'eas', 'sas', 'asj', 'oth', 'afr', 'global']:
        #     tag = 'maf_gnomad_' + str(pop)
        #     if arg_dict[tag]:
        #         conf_options['somatic_snv']['tumor_only'][tag] = float(arg_dict[tag])
        
        conf_options['somatic_snv']['mutational_signatures'] = {
            'run': int(arg_dict['estimate_signatures']),
            'mutation_limit': int(arg_dict['min_mutations_signatures']),
            'all_reference_signatures': int(arg_dict['all_reference_signatures']),
            'include_artefact_signatures': int(arg_dict['include_artefact_signatures']),
            'prevalence_reference_signatures': float(arg_dict['prevalence_reference_signatures'])
        }
        
        for fname in ['fname_cna_gene_tsv', 'fname_cna_segment_tsv', 
                      'fname_expression_tsv', 'fname_expression_outliers_tsv', 
                      'fname_rna_fusion_tsv', 'fname_maf_tsv', 
                      'fname_germline_tsv', 'fname_germline_yaml', 
                      'fname_expression_similarity_tsv', 'fname_tmb_tsv']:
            conf_options['molecular_data'][fname] = "None"
    
    
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
            #'pop_gnomad': str(arg_dict['pop_gnomad']),
            'max_af_gnomad': float(arg_dict['max_af_gnomad']),
            #'maf_upper_threshold': float(arg_dict['maf_upper_threshold']),
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

                if ',' not in conf_data['conf']['gene_panel']['panel_id']:
                    conf_data['conf']['gene_panel']['url'] = str(conf_data['conf']['gene_panel']['panel_genes'][0]['panel_url'])
                    conf_data['conf']['gene_panel']['description_trait'] = str(conf_data['conf']['gene_panel']['panel_genes'][0]['panel_name'])
                else:
                    names = set([str(x['panel_name']) for x in conf_data['conf']['gene_panel']['panel_genes']])
                    names2 = []
                    for n in names:
                        if 'ACMG' not in n and 'CPIC' not in n:
                            names2.append(n)
                    conf_data['conf']['gene_panel']['description'] = "Genomics England PanelApp - multiple panels (" + ', '.join(names2) + ")"
                    
                    
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


def verify_oncotree_code(oncotree_code, tumor_site, refdata_assembly_dir, logger=None):
    """Verify that a user-provided OncoTree code is valid and matches the
    selected PCGR tumor site. Returns the validated oncotree_code on success,
    or None if validation fails."""

    oncotree_fname = os.path.join(
        refdata_assembly_dir, "phenotype", "tsv", "phenotype_oncotree.tsv.gz")

    if not check_file_exists(oncotree_fname, logger):
        warn_message(
            f"OncoTree annotation file not found: {oncotree_fname} "
            "- skipping OncoTree code validation", logger)
        return oncotree_code

    oncotree_df = pd.read_csv(oncotree_fname, sep="\t", na_values=".")

    if 'ot_code' not in oncotree_df.columns or 'primary_site' not in oncotree_df.columns:
        warn_message(
            "OncoTree annotation file missing required columns "
            "('ot_code', 'primary_site') - skipping validation", logger)
        return oncotree_code

    valid_codes = oncotree_df['ot_code'].dropna().unique()
    if oncotree_code not in valid_codes:
        warn_message(
            f"User-provided OncoTree code '{oncotree_code}' is not a recognized "
            f"code in {oncotree_fname} - ignoring oncotree code for OncoKB annotation", logger)
        return None

    code_rows = oncotree_df.loc[oncotree_df['ot_code'] == oncotree_code]
    code_primary_sites = code_rows['primary_site'].dropna().unique()

    if len(code_primary_sites) == 0:
        warn_message(
            f"OncoTree code '{oncotree_code}' has no associated primary_site "
            "- proceeding with the provided code", logger)
        return oncotree_code

    if tumor_site not in code_primary_sites:
        warn_message(
            f"OncoTree code '{oncotree_code}' maps to primary site(s) "
            f"'{', '.join(code_primary_sites)}', which does not match the "
            f"selected PCGR tumor site '{tumor_site}' - ignoring oncotree code "
            "for OncoKB annotation", logger)
        return None

    return oncotree_code

