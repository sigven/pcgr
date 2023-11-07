#!/usr/bin/env python

from pcgr import pcgr_vars
from pcgr.utils import check_file_exists, error_message

import pandas as pd
import os
import gzip
import csv

def create_config(arg_dict, workflow = "PCGR"):
    
    config_options = {}
    if workflow == "PCGR" or workflow == "CPSR":
        config_options = {
            'sample_id': arg_dict['sample_id'],
            'genome_assembly': arg_dict['genome_assembly'],
            'debug': arg_dict['debug'],
            'pcgrr_conda': arg_dict['pcgrr_conda'],                
    	}
            
        config_options['vep'] = {
            'vep_buffer_size': int(arg_dict['vep_buffer_size']),
            'vep_pick_order': str(arg_dict['vep_pick_order']),
            'vep_n_forks': int(arg_dict['vep_n_forks']),
            'vep_no_intergenic': int(arg_dict['vep_no_intergenic']),
            'vep_regulatory': 1,
            'vep_gencode_all': int(arg_dict['vep_gencode_all']),            
        }
            
        config_options['visual_reporting'] = {
            'visual_theme': str(arg_dict['report_theme']),
            'nonfloating_toc': int(arg_dict['report_nonfloating_toc']),
        }
            
        config_options['other'] = {
            'vcfanno_n_proc': int(arg_dict['vcfanno_n_proc']),                                          
            'no_reporting': int(arg_dict['no_reporting']),
            'retained_info_tags': str(arg_dict['retained_info_tags']),
            'show_noncoding': int(arg_dict['show_noncoding']),
            'force_overwrite': int(arg_dict['force_overwrite'])
        }
    
    if workflow == 'PCGR':
        config_options['tumor_purity'] = 'NA'
        config_options['tumor_ploidy'] = 'NA'
        config_options['tumor_type'] = {
            'type': str(pcgr_vars.tsites[arg_dict['tsite']])
        }
        config_options['assay'] = str(arg_dict['assay'])
        config_options['clinicaltrials'] = {
            'run': int(arg_dict['include_trials'])
        }
        config_options['msi'] = {
            'run': int(arg_dict['estimate_msi_status'])
        }
        
        if not arg_dict['tumor_purity'] is None:
            config_options['tumor_purity'] = float(arg_dict['tumor_purity'])
        if not arg_dict['tumor_ploidy'] is None:
            config_options['tumor_ploidy'] = float(arg_dict['tumor_ploidy'])

        config_options['other']['vcf2maf'] = int(arg_dict['vcf2maf'])
        
        config_options['msigs'] = {
            'run': int(arg_dict['estimate_signatures']),
            'mutation_limit': int(arg_dict['min_mutations_signatures']),
            'all_reference_signatures': int(arg_dict['all_reference_signatures']),
            'include_artefact_signatures': int(arg_dict['include_artefact_signatures']),
            'prevalence_reference_signatures': int(arg_dict['prevalence_reference_signatures']),
        }

        config_options['tmb'] = {
            'run': int(arg_dict['estimate_tmb']),
            'target_size_mb': arg_dict['target_size_mb'],
            'algorithm': arg_dict['tmb_algorithm'],
        }

        config_options['cna'] = {
            'logR_homdel': float(arg_dict['logr_homdel']),
            'logR_gain': float(arg_dict['logr_gain']),
            'cna_overlap_pct': float(arg_dict['cna_overlap_pct']),
            'n_copy_ampl': int(arg_dict['n_copy_ampl'])
        }

        config_options['allelic_support'] = {
            'tumor_dp_min': int(arg_dict['tumor_dp_min']),
            'control_dp_min': int(arg_dict['control_dp_min']),
            'tumor_af_min': float(arg_dict['tumor_af_min']),
            'control_af_max': float(arg_dict['control_af_max']),
            'control_dp_tag': str(arg_dict['control_dp_tag']),
            'control_af_tag': str(arg_dict['control_af_tag']),
            'tumor_dp_tag': str(arg_dict['tumor_dp_tag']),
            'tumor_af_tag': str(arg_dict['tumor_af_tag']),
            'call_conf_tag': str(arg_dict['call_conf_tag']),
        }

        config_options['tumor_only'] = {
            'tumor_only': int(arg_dict['tumor_only']),
            'cell_line': int(arg_dict['cell_line']),
            'exclude_pon': int(arg_dict['exclude_pon']),
            'exclude_likely_hom_germline': int(arg_dict['exclude_likely_hom_germline']),
            'exclude_likely_het_germline': int(arg_dict['exclude_likely_het_germline']),
            'exclude_dbsnp_nonsomatic': int(arg_dict['exclude_dbsnp_nonsomatic']),
            'exclude_nonexonic': int(arg_dict['exclude_nonexonic']),
        }

        for pop in ['nfe', 'fin', 'amr', 'eas', 'sas', 'asj', 'oth', 'afr', 'global']:
            tag = 'maf_gnomad_' + str(pop)
            if arg_dict[tag]:
                config_options['tumor_only'][tag] = float(arg_dict[tag])
    
    if workflow == "CPSR":
        
        config_options['visual_reporting']['table_display'] = str(arg_dict['report_table_display'])
        config_options['gene_panel'] = {
            'panel_id': str(arg_dict['virtual_panel_id']),
            'description': 'Exploratory virtual gene panel (panel 0)',
            'custom_list_tsv': str(arg_dict['custom_list']),
            'custom_list_name': str(arg_dict['custom_list_name']),
            'custom_list_bed': 'None',
            'diagnostic_grade_only': int(arg_dict['diagnostic_grade_only'])
        }
            
        config_options['classification'] = {
            'gwas_findings': int(arg_dict['gwas_findings']),
            'secondary_findings': int(arg_dict['secondary_findings']),
            'pop_gnomad': str(arg_dict['pop_gnomad']),
            'maf_upper_threshold': float(arg_dict['maf_upper_threshold']),
            'classify_all': int(arg_dict['classify_all']),
            'clinvar_ignore_noncancer': int(arg_dict['clinvar_ignore_noncancer'])
        }
            

    return config_options


def populate_config_data(config_options: dict, db_dir: str, workflow = "PCGR", logger=None):
    
    conf_data = {}
    conf_data['sample_id'] = config_options['sample_id']
    conf_data['genome_assembly'] = config_options['genome_assembly']
    conf_data['software'] = {}
    conf_data['software']['pcgr_version'] = pcgr_vars.PCGR_VERSION
    conf_data['software']['cpsr_version'] = pcgr_vars.PCGR_VERSION
    
    ## add paths to annotated files (VCF/TSV, CNA)
    conf_data['pre_report_output'] = {}
    conf_data['pre_report_output']['path'] = config_options['output_dir']
    conf_data['pre_report_output']['fname_mut_vcf'] = config_options['annotated_vcf']
    conf_data['pre_report_output']['fname_mut_tsv'] = config_options['annotated_tsv']
    if workflow == "PCGR" and config_options['annotated_cna'] is not None:
        conf_data['pre_report_output']['fname_cna_tsv'] = config_options['annotated_cna']
        del config_options['annotated_cna']
    
    genome_assembly = config_options['genome_assembly']
    del config_options['sample_id']
    del config_options['genome_assembly']
    del config_options['output_dir']
    del config_options['annotated_vcf']
    del config_options['annotated_tsv']
    conf_data['config'] = config_options
    
    conf_data['reference_data'] = {}
    conf_data['reference_data']['version'] = pcgr_vars.DB_VERSION
    conf_data['reference_data']['path'] = \
        os.path.join(db_dir, "data", genome_assembly)
    
    metadata_pd = pd.DataFrame()
    conf_data['reference_data']['source_metadata'] = {}
    conf_data['phenotype'] = {}
    conf_data['phenotype']['oncotree'] = {}
    
    ## add metadata information for each data source
    for dtype in ['gene','phenotype','biomarker','drug','gwas','hotspot','other']:
        metadata_fname = os.path.join(
            db_dir, "data", conf_data['genome_assembly'],
            ".METADATA", "tsv", dtype + "_metadata.tsv")
        if check_file_exists(metadata_fname, logger):
            metadata_df = pd.read_csv(metadata_fname, sep="\t", na_values=".")
            metadata_df["source_type"] = dtype
            metadata_pd = metadata_pd.append(metadata_df, ignore_index=True)
    
    conf_data['reference_data']['source_metadata'] = metadata_pd.to_dict(orient='records')

    oncotree_fname = os.path.join(db_dir, "data", genome_assembly,
                                  "phenotype","tsv","phenotype_onco.tsv.gz")
    
    if check_file_exists(oncotree_fname, logger):
        oncotree_df = pd.read_csv(oncotree_fname, sep="\t", na_values=".")
        if workflow == "CPSR":
            conf_data['phenotype']['oncotree'] = \
                oncotree_df[oncotree_df['primary_site'].isna()].to_dict(orient='records')
        else:
            conf_data['phenotype']['oncotree'] = pd.DataFrame()
    
    if workflow == "CPSR":
        if conf_data['config']['gene_panel']['panel_id'] is not None and \
            conf_data['config']['gene_panel']['diagnostic_grade_only'] is not None:
                
            if conf_data['config']['gene_panel']['panel_id'] == "-1":
                conf_data['config']['gene_panel']['description'] = 'User-defined panel (custom geneset from panel 0)'
            else:
                if ',' in conf_data['config']['gene_panel']['panel_id']:
                    conf_data['config']['gene_panel']['description'] = \
                        'Genomics England PanelApp - multiple panels (' + conf_data['config']['gene_panel']['panel_id'] + ')'
                else:
                    if conf_data['config']['gene_panel']['panel_id'] != "0":
                        conf_data['config']['gene_panel']['description'] = \
                            'Genomics England PanelApp - panel ' + conf_data['config']['gene_panel']['panel_id']
        
            conf_data['config']['gene_panel']['panel_genes'] = set_virtual_target_genes(
                conf_data['config']['gene_panel']['panel_id'], 
                db_dir, 
                conf_data['genome_assembly'],
                conf_data['config']['gene_panel']['diagnostic_grade_only'], 
                conf_data['config']['gene_panel']['custom_list_bed'],
                logger)
        

    return(conf_data)

def set_virtual_target_genes(panel_id: str, db_dir: str, genome_assembly: str, diagnostic_grade_only: bool, custom_list_bed: str, logger=None):
    
    all_panels_fname = os.path.join(
        db_dir, "data", genome_assembly,
        "gene","tsv","gene_virtual_panel", 
        "gene_virtual_panel.tsv.gz")
    
    all_virtual_panels = pd.DataFrame()
    panel_targets = pd.DataFrame()
    if check_file_exists(all_panels_fname, logger):    
        all_virtual_panels = pd.read_csv(all_panels_fname, sep="\t", na_values=".")
        all_virtual_panels['id'] = all_virtual_panels['id'].astype(str)
        all_virtual_panels = \
            all_virtual_panels.drop(['gepa_phenotype', 'genome_build','gepa_penetrance','gepa_moi'], axis=1)
        all_virtual_panels = \
            all_virtual_panels.rename(
                columns={'gepa_panel_version': 'panel_version', 
                         'gepa_confidence_level': 'confidence_level',
                         'gepa_panel_url': 'panel_url',
                         'gepa_panel_name': 'panel_name',
                         'gepa_panel_id': 'panel_id'})


    if panel_id == "-1":
        if not custom_list_bed == 'None':
            custom_ensembl_gene_ids = {}
            if check_file_exists(custom_list_bed, logger):
                with open(custom_list_bed) as f:
                    reader = csv.reader(f, delimiter='\t')
                    for row in reader:
                        ensembl_gene_id = row[3].split('|')[1]
                        custom_ensembl_gene_ids[ensembl_gene_id] = 1
                f.close()
            
            panel_targets = all_virtual_panels[all_virtual_panels['id'] == "0"].copy()
            panel_targets = panel_targets[panel_targets['ensembl_gene_id'].isin(custom_ensembl_gene_ids)]
            if len(panel_targets) == 0:
                err_msg = "Custom list of genes from CPSR superpanel (panel 0) should be provided as Ensembl " + \
                    "gene identifiers, '" + str(custom_ensembl_gene_ids) + "' were not find in panel 0"
                error_message(err_msg, logger)
            else:
                panel_targets.loc[:,'confidence_level'] = -1
                panel_targets.loc[:,'panel_id'] = None
                panel_targets.loc[:,'panel_url'] = None
                panel_targets.loc[:,'panel_version'] = None
            
    else:
        panel_ids = panel_id.split(',')            
        panel_targets = all_virtual_panels[all_virtual_panels['id'].isin(panel_ids)].copy()
        if diagnostic_grade_only:
            panel_targets = panel_targets[panel_targets['confidence_level'] >= 3]
        if len(panel_targets) == 0:
            m = 1 ## issue error
        else:
            if len(panel_ids) > 0:
                panel_targets.loc[:,'confidence_level'] = 5  
                    
    return panel_targets.to_dict(orient='records')
    
