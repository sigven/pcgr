#!/usr/bin/env python

from pcgr import pcgr_vars
from pcgr.utils import check_file_exists, error_message

import pandas as pd
import os
import gzip
import csv

def create_config(arg_dict, workflow = "PCGR"):
    
    conf_options = {}
    if workflow == "PCGR" or workflow == "CPSR":
        conf_options = {
            'sample_id': arg_dict['sample_id'],
            'genome_assembly': arg_dict['genome_assembly'],
            'debug': arg_dict['debug'],
            'pcgrr_conda': arg_dict['pcgrr_conda']                
    	}
        
        conf_options['vep'] = {
            'vep_buffer_size': int(arg_dict['vep_buffer_size']),
            'vep_pick_order': str(arg_dict['vep_pick_order']),
            'vep_n_forks': int(arg_dict['vep_n_forks']),
            'vep_no_intergenic': int(arg_dict['vep_no_intergenic']),
            'vep_regulatory': 1,
            'vep_gencode_basic': int(arg_dict['vep_gencode_basic'])            
        }
            
        conf_options['visual_reporting'] = {
            'visual_theme': str(arg_dict['report_theme']),
            'nonfloating_toc': int(arg_dict['report_nonfloating_toc'])
        }
            
        conf_options['other'] = {
            'vcfanno_n_proc': int(arg_dict['vcfanno_n_proc']),                                          
            'no_reporting': int(arg_dict['no_reporting']),
            'retained_vcf_info_tags': str(arg_dict['retained_info_tags']),
            'show_noncoding': not int(arg_dict['ignore_noncoding']),
            'force_overwrite': int(arg_dict['force_overwrite'])
        }
        conf_options['sample_properties'] = {}

    
    if workflow == 'PCGR':
        conf_options['assay_properties'] = {}
        conf_options['sample_properties']['purity'] = 'NA'
        conf_options['sample_properties']['ploidy'] = 'NA'
        conf_options['sample_properties']['site'] = str(pcgr_vars.tsites[arg_dict['tsite']])
        conf_options['assay_properties']['type'] = str(arg_dict['assay'])
        conf_options['assay_properties']['vcf_tumor_only'] = 0
        conf_options['assay_properties']['mode'] = "Tumor-Control"
        conf_options['assay_properties']['cell_line'] = int(arg_dict['cell_line'])
        conf_options['assay_properties']['effective_target_size_mb'] = float(arg_dict['effective_target_size_mb'])
        if int(arg_dict['tumor_only']) == 1:
            conf_options['assay_properties']['vcf_tumor_only'] = 1
            conf_options['assay_properties']['mode'] = "Tumor-Only"
            if int(arg_dict['cell_line']) == 1:
                conf_options['assay_properties']['mode'] = "Cell line (Tumor-Only)"
        
        if not arg_dict['tumor_purity'] is None:
            conf_options['sample_properties']['tumor_purity'] = float(arg_dict['tumor_purity'])
        if not arg_dict['tumor_ploidy'] is None:
            conf_options['sample_properties']['tumor_ploidy'] = float(arg_dict['tumor_ploidy'])
            
        conf_options['clinicaltrials'] = {
            'run': int(arg_dict['include_trials'])
        }
        conf_options['other']['vcf2maf'] = int(arg_dict['vcf2maf'])
        conf_options['somatic_cna'] = {            
            'cna_overlap_pct': float(arg_dict['cna_overlap_pct']),
            'n_copy_gain': int(arg_dict['n_copy_gain'])
        }
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
            'prevalence_reference_signatures': int(arg_dict['prevalence_reference_signatures'])
        }

    
    if workflow == "CPSR":        
        conf_options['sample_properties']['phenotype'] = 'None'
        conf_options['sample_properties']['site'] = 'Hereditary (blood)'
        conf_options['sample_properties']['genotypes_available'] = 0
        conf_options['visual_reporting']['table_display'] = str(arg_dict['report_table_display'])
        conf_options['gene_panel'] = {
            'panel_id': str(arg_dict['virtual_panel_id']),
            'description': 'Exploratory virtual gene panel (panel 0)',
            'custom_list_tsv': str(arg_dict['custom_list']),
            'custom_list_name': str(arg_dict['custom_list_name']),
            'custom_list_bed': 'None',
            'diagnostic_grade_only': int(arg_dict['diagnostic_grade_only'])
        }
            
        conf_options['variant_classification'] = {
            'gwas_findings': int(arg_dict['gwas_findings']),
            'secondary_findings': int(arg_dict['secondary_findings']),
            'pop_gnomad': str(arg_dict['pop_gnomad']),
            'maf_upper_threshold': float(arg_dict['maf_upper_threshold']),
            'classify_all': int(arg_dict['classify_all']),
            'clinvar_report_noncancer': int(arg_dict['clinvar_report_noncancer'])
        }
            

    return conf_options


def populate_config_data(conf_options: dict, db_dir: str, workflow = "PCGR", logger=None):
    
    conf_data = {}
    conf_data['sample_id'] = conf_options['sample_id']
    conf_data['output_dir'] = conf_options['output_dir']
    conf_data['workflow'] = workflow
    conf_data['genome_assembly'] = conf_options['genome_assembly']
    conf_data['software'] = {}
    conf_data['software']['pcgr_version'] = pcgr_vars.PCGR_VERSION
    conf_data['software']['cpsr_version'] = pcgr_vars.PCGR_VERSION
    
    ## add paths to annotated files (VCF/TSV, CNA, TMB)
    conf_data['molecular_data'] = {}
    conf_data['molecular_data']['fname_mut_vcf'] = conf_options['annotated_vcf']
    conf_data['molecular_data']['fname_mut_tsv'] = conf_options['annotated_tsv']
    
    conf_data['molecular_data']['fname_cna_tsv'] = "None"
    if workflow == "PCGR" and conf_options['annotated_cna'] != "None":
        conf_data['molecular_data']['fname_cna_tsv'] = conf_options['annotated_cna']
        del conf_options['annotated_cna']
    
    conf_data['molecular_data']['fname_tmb'] = "None"
    if workflow == "PCGR" and conf_options['fname_tmb'] != "None":
        conf_data['molecular_data']['fname_tmb'] = conf_options['fname_tmb']
        del conf_options['fname_tmb']
    
    genome_assembly = conf_options['genome_assembly']
    del conf_options['sample_id']
    del conf_options['genome_assembly']
    del conf_options['output_dir']
    del conf_options['annotated_vcf']
    del conf_options['annotated_tsv']
    conf_data['conf'] = conf_options
    
    conf_data['reference_data'] = {}
    conf_data['reference_data']['version'] = pcgr_vars.DB_VERSION
    conf_data['reference_data']['path'] = \
        os.path.join(db_dir, "data", genome_assembly)
    
    metadata_pd = pd.DataFrame()
    conf_data['reference_data']['source_metadata'] = {}
    conf_data['conf']['sample_properties']['phenotype'] = {}
    
    ## add metadata information for each data source
    for dtype in ['gene','phenotype','biomarker','drug','gwas','hotspot','other']:
        metadata_fname = os.path.join(
            db_dir, "data", conf_data['genome_assembly'],
            ".METADATA", "tsv", dtype + "_metadata.tsv")
        if check_file_exists(metadata_fname, logger):
            metadata_df = pd.read_csv(metadata_fname, sep="\t", na_values=".")
            metadata_df["source_type"] = dtype
            metadata_pd = metadata_pd._append(metadata_df, ignore_index=True)
    
    conf_data['reference_data']['source_metadata'] = metadata_pd.to_dict(orient='records')

    oncotree_fname = os.path.join(db_dir, "data", genome_assembly,
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
            else:
                if ',' in conf_data['conf']['gene_panel']['panel_id']:
                    conf_data['conf']['gene_panel']['description'] = \
                        'Genomics England PanelApp - multiple panels (' + conf_data['conf']['gene_panel']['panel_id'] + ')'
                else:
                    if conf_data['conf']['gene_panel']['panel_id'] != "0":
                        conf_data['conf']['gene_panel']['description'] = \
                            'Genomics England PanelApp - panel ' + conf_data['conf']['gene_panel']['panel_id']
        
            conf_data['conf']['gene_panel']['panel_genes'] = set_virtual_target_genes(
                conf_data['conf']['gene_panel']['panel_id'], 
                db_dir, 
                conf_data['genome_assembly'],
                conf_data['conf']['gene_panel']['diagnostic_grade_only'], 
                conf_data['conf']['gene_panel']['custom_list_tsv'],
                bool(conf_data['conf']['variant_classification']['secondary_findings']),
                logger)
        

    return(conf_data)

def set_virtual_target_genes(panel_id: str, db_dir: str, genome_assembly: str, diagnostic_grade_only: bool, 
                             custom_list_tsv: str, secondary_findings: bool, logger=None):
    
    all_panels_fname = os.path.join(
        db_dir, "data", genome_assembly,
        "gene","tsv","gene_virtual_panel", 
        "gene_virtual_panel.tsv.gz")
    
    all_cpg_fname = os.path.join(
        db_dir, "data", genome_assembly,
        "gene","tsv","gene_cpg", 
        "gene_cpg.tsv.gz")
        
    all_virtual_panels = pd.DataFrame()
    all_secondary_finding_targets = pd.DataFrame()
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
    
    if secondary_findings is True:
        if check_file_exists(all_cpg_fname, logger):
            all_secondary_finding_targets = pd.read_csv(all_cpg_fname, sep="\t", na_values=".")    
            all_secondary_finding_targets = \
                all_secondary_finding_targets[all_secondary_finding_targets.cpg_source.str.match('ACMG_SF')]    
            all_secondary_finding_targets = \
                all_secondary_finding_targets.drop(
                    ['cpg_phenotypes', 'cpg_cancer_cui','cpg_source','cpg_syndrome_cui'], axis=1)
            all_secondary_finding_targets = \
                all_secondary_finding_targets.rename(
                    columns={'cpg_moi': 'moi',
                             'cpg_mod': 'mod'}) 
            all_secondary_finding_targets['confidence_level'] = -1
            all_secondary_finding_targets['panel_id'] = float(3.2)
            all_secondary_finding_targets['panel_url'] = "https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/"
            all_secondary_finding_targets['panel_version'] = float(3.2)
            all_secondary_finding_targets['panel_name'] = "ACMG_SF"
            all_secondary_finding_targets['primary_target'] = False
            all_secondary_finding_targets['id'] = str("-2")
            for f in ['panel_url', 'panel_name','moi','mod','symbol','ensembl_gene_id']:
                all_secondary_finding_targets[f] = all_secondary_finding_targets[f].astype(str)


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
    
    if len(all_secondary_finding_targets) > 0:        
        all_targets = pd.concat([panel_targets, all_secondary_finding_targets], axis=0)
    
    return all_targets.to_dict(orient='records')
    
