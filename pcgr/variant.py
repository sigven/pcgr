#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import warnings
from pandas.api.types import is_float_dtype, is_integer_dtype


from pcgr import pcgr_vars
from pcgr.annoutils import read_infotag_file, get_vcfanno_tracks, read_vcfanno_tag_file
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def append_annotations(vcf2tsv_gz_fname: str, refdata_assembly_dir: str, logger):
    
    """
    Function that appends various variant annotations to a vcf2tsv file
    It assumes that the VCF2TSV file has the following columns:
    CHROM, POS, REF, ALT, CLINVAR_MSID, PFAM_DOMAIN, ENTREZGENE, EFFECT_PREDICTIONS
    
    The function appends the following columns to the input VCF2TSV file:
    CLINVAR_TRAITS_ALL, PFAM_DOMAIN_NAME, GENENAME
    
    CLINVAR_TRAITS_ALL is a column that concatenates the ClinVar trait names
    associated with a given variant. PFAM_DOMAIN_NAME is a column that
    concatenates the PFAM domain names associated with a given variant.
    GENENAME is a column that provides the gene name associated with the
    ENTREZGENE ID associated with a given variant.
    
    The function requires the following files from the reference data directory:
    (1) variant/clinvar/clinvar.tsv.gz
    (2) misc/protein_domain/protein_domain.tsv.gz
    (3) gene/transcript_xref/gene_transcript_xref.tsv.gz
    
    Parameters
    ----------
    vcf2tsv_gz_fname : str
        Name of input VCF2TSV file
    refdata_assembly_dir : str
        Name of reference data directory
    logger : logging.Logger
        Logger object for logging messages
    
    Returns
    -------
    pd.DataFrame
        Input VCF2TSV DataFrame with appended annotations
    """
    clinvar_tsv_fname = os.path.join(refdata_assembly_dir, 'variant','tsv','clinvar', 'clinvar.tsv.gz')
    protein_domain_tsv_fname = os.path.join(refdata_assembly_dir, 'misc','tsv','protein_domain', 'protein_domain.tsv.gz') 
    gene_xref_tsv_fname = os.path.join(refdata_assembly_dir, 'gene','tsv','gene_transcript_xref', 'gene_transcript_xref.tsv.gz')
    vcf2tsv_df = None
    clinvar_data_df = None
    
    num_recs_with_clinvar_hits = 0
    
    if os.path.exists(vcf2tsv_gz_fname):
        vcf2tsv_df = pd.read_csv(
            vcf2tsv_gz_fname, skiprows=[0], sep="\t", na_values=".",
            low_memory = False)
        
        if {'CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN','ENTREZGENE','EFFECT_PREDICTIONS'}.issubset(vcf2tsv_df.columns):
            for elem in ['CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN','EFFECT_PREDICTIONS']:
                vcf2tsv_df = vcf2tsv_df.astype({elem:'string'})            
            vcf2tsv_df['CLINVAR_MSID'] = vcf2tsv_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df['PFAM_DOMAIN'] = vcf2tsv_df['PFAM_DOMAIN'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df['EFFECT_PREDICTIONS'] = \
                vcf2tsv_df['EFFECT_PREDICTIONS'].str.replace(
                    "\\.&|\\.$", "NA&", regex = True).str.replace(
                        "&$", "", regex = True).str.replace("&", ", ", regex = True)
            vcf2tsv_df["VAR_ID"] = vcf2tsv_df["CHROM"].str.cat(
                vcf2tsv_df["POS"], sep = "_").str.cat(
                    vcf2tsv_df["REF"], sep = "_").str.cat(
                        vcf2tsv_df["ALT"], sep = "_")
           
            ## check number of variants with ClinVar ID's
            num_recs_with_clinvar_hits = vcf2tsv_df["CLINVAR_MSID"].notna().sum()
            ## check number of variants with PFAM ID's
            num_recs_with_pfam_hits = vcf2tsv_df["PFAM_DOMAIN"].notna().sum()
            ## check number of variants with Ensembl gene ID's
            num_recs_with_entrez_hits = vcf2tsv_df["ENTREZGENE"].notna().sum()
    
            ## merge variant set with ClinVar trait and variant origin annotations
            if num_recs_with_clinvar_hits > 0:
                
                if {'CLINVAR_TRAITS_ALL'}.issubset(vcf2tsv_df.columns):
                    vcf2tsv_df.drop('CLINVAR_TRAITS_ALL', inplace=True, axis=1)
        
                if os.path.exists(clinvar_tsv_fname):
                    clinvar_data_df = pd.read_csv(
                        clinvar_tsv_fname, sep="\t", 
                        usecols=["variation_id","origin_simple","VAR_ID","trait","classification","conflicted"],
                        low_memory = False)
                    clinvar_data_df['classification'] = clinvar_data_df['classification'].str.replace("_", " ", regex = True)
                    clinvar_data_df.loc[(clinvar_data_df['classification'] == "CPP Unknown") & (clinvar_data_df['conflicted'] == 1),"classification"] = \
                        "Conflicting classifications of pathogenicity"
                    clinvar_data_df.loc[(clinvar_data_df['classification'] == "CPP Unknown") & (clinvar_data_df['conflicted'] == 0),"classification"] = \
                        "No classification"
                    clinvar_data_df['CLINVAR_TRAITS_ALL'] = clinvar_data_df['classification'].str.cat(
                        clinvar_data_df['origin_simple'].str.capitalize().str.cat(
                        clinvar_data_df['trait'], sep = " - "), sep = " - ")
                    clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['variation_id']
                    clinvar_data_df['CLINVAR_TRAITS'] = clinvar_data_df['trait']
                    clinvar_data_df = clinvar_data_df.astype({'CLINVAR_MSID':'string'})
                    clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                    clinvar_data_df = clinvar_data_df[['VAR_ID','CLINVAR_MSID','CLINVAR_TRAITS','CLINVAR_TRAITS_ALL']]
                    vcf2tsv_df = vcf2tsv_df.merge(
                        clinvar_data_df, left_on=["VAR_ID", "CLINVAR_MSID"], right_on=["VAR_ID", "CLINVAR_MSID"], how="left")
                else:
                    logger.error(f"Could not find {clinvar_tsv_fname} needed for ClinVar variant annotation - exiting")
            else:
                vcf2tsv_df['CLINVAR_TRAITS_ALL'] = '.'
                vcf2tsv_df['CLINVAR_MSID'] = '.'
                vcf2tsv_df['CLINVAR_TRAITS'] = '.'
                
            
            ## merge variant set with PFAM domain annotations
            if num_recs_with_pfam_hits > 0:
                vcf2tsv_df.drop(['PFAM_DOMAIN_NAME', 'PFAM_ENTRY_LOCATIONS'], axis=1, errors='ignore', inplace=True)
                if os.path.exists(protein_domain_tsv_fname):
                    pfam_data_df = pd.read_csv(
                        protein_domain_tsv_fname, sep="\t",
                        usecols=["pfam_id", "pfam_name", "pfam_entry_locations", "entrezgene"]
                    ).drop_duplicates()
                    pfam_data_df.rename(
                        columns={
                            'pfam_id': 'PFAM_DOMAIN',
                            'pfam_name': 'PFAM_DOMAIN_NAME',
                            'pfam_entry_locations': 'PFAM_ENTRY_LOCATIONS',
                            'entrezgene': 'ENTREZGENE'
                        }, inplace=True)
                    vcf2tsv_df = vcf2tsv_df.merge(
                        pfam_data_df,
                        on=["PFAM_DOMAIN", "ENTREZGENE"],
                        how="left"
                    )
                else:
                    logger.error(f"Could not find {protein_domain_tsv_fname} needed for PFAM domain annotation - exiting")
            else:
                vcf2tsv_df['PFAM_DOMAIN_NAME'] = '.'
                vcf2tsv_df['PFAM_ENTRY_LOCATIONS'] = '.'
            
            if num_recs_with_entrez_hits > 0:
                
                if {'GENENAME'}.issubset(vcf2tsv_df.columns):
                    vcf2tsv_df.drop('GENENAME', inplace=True, axis=1)
                
                if os.path.exists(gene_xref_tsv_fname):
                    gene_xref_df = pd.read_csv(
                        gene_xref_tsv_fname, sep="\t", na_values=".", 
                        low_memory=False, usecols=["entrezgene","name"])
                    gene_xref_df = gene_xref_df[gene_xref_df['entrezgene'].notnull()].drop_duplicates()
                    gene_xref_df.rename(columns = {'entrezgene':'ENTREZGENE', 'name':'GENENAME'}, inplace = True)                
                    vcf2tsv_df = vcf2tsv_df.merge(gene_xref_df, left_on=["ENTREZGENE"], right_on=["ENTREZGENE"], how="left")
                else:
                    logger.error(f"Could not find {gene_xref_tsv_fname} needed for gene name annotation - exiting")
            else:
                vcf2tsv_df['GENENAME'] = '.'
                
    return(vcf2tsv_df)

def set_allelic_support(variant_set: pd.DataFrame, allelic_support_tags: dict, logger: None) -> pd.DataFrame:
    """
    Set allelic support for variants
    """
    tag_to_info_val = {'control_dp_tag': 'DP_CONTROL', 
                        'control_af_tag': 'VAF_CONTROL', 
                        'tumor_dp_tag': 'DP_TUMOR', 
                        'tumor_af_tag': 'VAF_TUMOR'}
    
    for t in tag_to_info_val.keys():
        if tag_to_info_val[t] in variant_set.columns:
            variant_set[tag_to_info_val[t]] = np.nan
    
    variant_set['AD_TUMOR'] = np.nan
    variant_set['AD_CONTROL'] = np.nan

    for t in tag_to_info_val.keys():
        if allelic_support_tags[t] != "_NA_":
            if {allelic_support_tags[t], tag_to_info_val[t]}.issubset(variant_set.columns):
                if '_af_' in t:
                    if variant_set[variant_set[allelic_support_tags[t]].isna()].empty:
                        variant_set[tag_to_info_val[t]] = variant_set[allelic_support_tags[t]].astype(float).round(4)
                    else:
                        logger.warning(f"Unable to set allelic support for '{t}' - missing values deteted")
                else:
                    if variant_set[variant_set[allelic_support_tags[t]].isna()].empty:

                        variant_set[tag_to_info_val[t]] = variant_set[allelic_support_tags[t]].astype(int)
                    else:
                        logger.warning(f"Unable to set read depth support for '{t}' - missing values deteted")
    
    if variant_set[variant_set['DP_TUMOR'].isna()].empty and variant_set[variant_set['VAF_TUMOR'].isna()].empty:
        # Calculate AD_TUMOR (floor of DP * VAF)
        variant_set["AD_TUMOR"]   = np.floor(variant_set["DP_TUMOR"] * variant_set["VAF_TUMOR"]).astype(int)
    
    if variant_set[variant_set['DP_CONTROL'].isna()].empty and variant_set[variant_set['VAF_CONTROL'].isna()].empty:        
        # Calculate AD_CONTROL (floor of DP * VAF)
        variant_set["AD_CONTROL"] = np.floor(variant_set["DP_CONTROL"] * variant_set["VAF_CONTROL"]).astype(int)       

    return variant_set

def calculate_tmb(variant_set: pd.DataFrame, tumor_dp_min: int, tumor_af_min: float, tumor_ad_min: int,
                  target_size_mb: float, sample_id: str, tmb_fname: str, logger) -> int:
    """
    Calculate TMB for somatic variant set
    """
    logger.info(f"Calculating tumor mutational load/TMB from somatic call set - effective coding target size (Mb): {target_size_mb}")
    
    if {'CONSEQUENCE','VAF_TUMOR','DP_TUMOR','AD_TUMOR','VARIANT_CLASS'}.issubset(variant_set.columns):
        tmb_data_set = \
            variant_set[['CONSEQUENCE','DP_TUMOR','AD_TUMOR','VAF_TUMOR','VARIANT_CLASS']]
        
        n_rows_raw = len(tmb_data_set)
        active_filters = {}
        if float(tumor_af_min) > 0:
            active_filters['tumor_af'] = f'AF < {tumor_af_min}'
            ## filter variant set to those with VAF_TUMOR > tumor_af_min
            if variant_set[variant_set['VAF_TUMOR'].isna()].empty:
                tmb_data_set = tmb_data_set.loc[tmb_data_set['VAF_TUMOR'].astype(float) >= float(tumor_af_min),:]
                #logger.info(f'Removing n = {n_rows_raw - len(tmb_data_set)} variants with allelic fraction (tumor sample) < {tumor_af_min}')
            else:
                logger.warning("Cannot filter on sequencing depth - 'VAF_TUMOR' contains missing values")
        if int(tumor_dp_min) > 0:
            active_filters['tumor_dp'] = f'DP < {tumor_dp_min}'
            ## filter variant set to those with DP_TUMOR > tumor_dp_min
            if variant_set[variant_set['DP_TUMOR'].isna()].empty:
                tmb_data_set = tmb_data_set.loc[tmb_data_set['DP_TUMOR'].astype(int) >= int(tumor_dp_min),:]
            else:
                logger.warning("Cannot filter on sequencing depth - 'DP_TUMOR' contains missing values")

        if int(tumor_ad_min) > 0:
            active_filters['tumor_ad'] = f'AD < {tumor_ad_min}'
            ## filter variant set to those with AD_TUMOR > tumor_ad_min
            if variant_set[variant_set['AD_TUMOR'].isna()].empty:
                tmb_data_set = tmb_data_set.loc[tmb_data_set['AD_TUMOR'].astype(int) >= int(tumor_ad_min),:]                
            else:
                logger.warning("Cannot filter on allelic depth - 'AD_TUMOR' contains missing values")
        
        tmb_categories = ['missense_only','coding_and_silent','coding_non_silent']
        pct_tmb_filtered = round(float((n_rows_raw - len(tmb_data_set)) / n_rows_raw * 100),2)
        if len(active_filters.keys()) > 0:
            logger.info(f'Excluded n = {n_rows_raw - len(tmb_data_set)} ({pct_tmb_filtered}% of all) tumor '
                        f'variants with {", ".join(active_filters.values())} for TMB calculation')
        
        regex_data = {}
        regex_data['missense_only'] = pcgr_vars.CSQ_MISSENSE_PATTERN
        regex_data['coding_and_silent'] = pcgr_vars.CSQ_CODING_SILENT_PATTERN
        regex_data['coding_non_silent'] = pcgr_vars.CSQ_CODING_PATTERN
        
        tmb_count_data = {}
        tmb_estimates = {}
        
        for cat in tmb_categories:
            tmb_count_data[cat] = 0
            tmb_estimates[cat] = 0
        
        if tmb_data_set.empty is False:
            for cat in tmb_categories:
                tmb_count_data[cat] = len(tmb_data_set.loc[tmb_data_set['CONSEQUENCE'].str.match(regex_data[cat]),:])
                tmb_estimates[cat] = round(float(tmb_count_data[cat] / target_size_mb), 4)
            
            tmb_recs = [[sample_id, n_rows_raw, 'TMB_missense_only', regex_data['missense_only'], target_size_mb, 
                        tumor_dp_min, tumor_af_min, tmb_count_data['missense_only'], tmb_estimates['missense_only'], 'mutations/Mb'],
                        [sample_id, n_rows_raw, 'TMB_coding_and_silent', regex_data['coding_and_silent'], target_size_mb, 
                        tumor_dp_min, tumor_af_min, tmb_count_data['coding_and_silent'], tmb_estimates['coding_and_silent'], 'mutations/Mb'],
                        [sample_id, n_rows_raw, 'TMB_coding_non_silent', regex_data['coding_non_silent'], target_size_mb, 
                        tumor_dp_min, tumor_af_min, tmb_count_data['coding_non_silent'], tmb_estimates['coding_non_silent'], 'mutations/Mb']]
                
            df = pd.DataFrame(
                tmb_recs, 
                columns=['SAMPLE_ID',
                        'N_SOMATIC_VARIANTS',
                        'TMB_MEASURE',
                        'TMB_CSQ_REGEX',
                        'TMB_TARGET_SIZE_MB',
                        'TMB_DP_MIN',
                        'TMB_AF_MIN',
                        'TMB_N_VARIANTS',
                        'TMB_ESTIMATE',
                        'TMB_UNIT'])
            df.to_csv(tmb_fname, sep="\t", header=True, index=False)
        else:
            logger.warning("No somatic variants left after depth/VAF filtering - TMB calculation not possible")

def clean_annotations(variant_set: pd.DataFrame, yaml_data: dict, logger) -> pd.DataFrame:
    
    refdata_assembly_dir = yaml_data["reference_data"]["path"]    
            
    if {'Consequence','EFFECT_PREDICTIONS','CLINVAR_CONFLICTED'}.issubset(variant_set.columns):
        variant_set.rename(columns = {'Consequence':'CONSEQUENCE'}, inplace = True)
        variant_set['clinvar_conflicted_bool'] = True
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] == 1, "clinvar_conflicted_bool"] = True
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] != 1, "clinvar_conflicted_bool"] = False
        variant_set.drop('CLINVAR_CONFLICTED', inplace=True, axis=1)        
        variant_set.rename(columns = {'clinvar_conflicted_bool':'CLINVAR_CONFLICTED'}, inplace = True)
        
    
    if not {'VCF_SAMPLE_ID'}.issubset(variant_set.columns):
        variant_set['VCF_SAMPLE_ID'] = str(yaml_data['sample_id'])
    variant_set['SAMPLE_ID'] = str(yaml_data['sample_id'])
    variant_set['GENOME_VERSION'] = yaml_data['genome_assembly']
    if {'CHROM','POS','REF','ALT',}.issubset(variant_set.columns):      
        variant_set['GENOMIC_CHANGE'] = variant_set['CHROM'].astype(str) + ":g." + variant_set['POS'].astype(str) + \
            variant_set['REF'].astype(str) + ">" + variant_set['ALT'].astype(str)
    if not {'PANEL_OF_NORMALS'}.issubset(variant_set.columns):
        variant_set['PANEL_OF_NORMALS'] = False


    ## Read VCF INFO tag definitions - find all tags that are defined as Integer and ensure that they are formatted as such
    vcf_infotags = {}
    vcf_infotags['other'] = read_infotag_file(os.path.join(refdata_assembly_dir, 'vcf_infotags_other.tsv'), scope = yaml_data['workflow'].lower())
    vcf_infotags['vep'] = read_infotag_file(os.path.join(refdata_assembly_dir, 'vcf_infotags_vep.tsv'), scope = "vep")
    vcf_infotags['other'].update(vcf_infotags['vep'])

    vcf_tags_integer = []
    vcfanno_tracks = get_vcfanno_tracks(refdata_assembly_dir)
    for track in vcfanno_tracks['tags_fname'].keys():
        infotags_vcfanno = read_vcfanno_tag_file(vcfanno_tracks['tags_fname'][track], logger)
        for k in infotags_vcfanno.keys():
            if infotags_vcfanno[k]['type'] == 'Integer' and infotags_vcfanno[k]['number'] == '1':
                vcf_tags_integer.append(infotags_vcfanno[k]['tag'])

    for k in vcf_infotags['other'].keys():
        if vcf_infotags['other'][k]['type'] == 'Integer':
            vcf_tags_integer.append(vcf_infotags['other'][k]['tag'])            

    for vcf_info_tag in sorted(vcf_tags_integer):
        if vcf_info_tag in variant_set.columns:
            if is_float_dtype(variant_set[vcf_info_tag]) or is_integer_dtype(variant_set[vcf_info_tag]):
                variant_set.loc[variant_set[vcf_info_tag].isna(),vcf_info_tag] = -123456789           
                variant_set[vcf_info_tag] = variant_set[vcf_info_tag].apply(lambda x: str(int(x)) )
                variant_set.loc[variant_set[vcf_info_tag] == "-123456789", vcf_info_tag] = np.nan
        
    
    return variant_set


def reduce_variants_for_report(
        output_pass_tsv_gz: str,
        output_pass_raw_tsv_gz: str,
        yaml_data: dict,
        debug: bool,
        logger) -> None:
    """
    When the annotated variant count exceeds MAX_VARIANTS_FOR_REPORT, progressively
    strip low-impact consequence classes so the reporting step stays tractable:

      1. Remove intergenic and intronic variants.
      2. If still above the limit, also remove upstream_gene / downstream_gene variants.

    The original gzipped TSV is renamed to *output_pass_raw_tsv_gz* before the
    filtered version is written to *output_pass_tsv_gz*.  Does nothing when the
    variant count is within the limit.
    """
    from pcgr.utils import check_subprocess, pd_to_csv

    if not os.path.exists(output_pass_tsv_gz):
        return

    var_data = pd.read_csv(output_pass_tsv_gz, sep='\t', low_memory=False, na_values='.')
    n_total = len(var_data)

    if n_total <= pcgr_vars.MAX_VARIANTS_FOR_REPORT:
        return

    logger.info(
        f'Number of raw variants in input VCF ({n_total}) exceeds '
        f'{pcgr_vars.MAX_VARIANTS_FOR_REPORT} - '
        f'intergenic/intronic variants will be excluded prior to reporting'
    )

    # Step 1: drop intergenic and intronic
    var_filtered = var_data[
        ~var_data['CONSEQUENCE'].str.contains('^intron', na=False) &
        ~var_data['CONSEQUENCE'].str.contains('^intergenic', na=False)
    ]
    logger.info(f'Number of intergenic/intronic variants excluded: {n_total - len(var_filtered)}')

    # Step 2: drop upstream/downstream if still above limit
    if len(var_filtered) > pcgr_vars.MAX_VARIANTS_FOR_REPORT:
        var_filtered_prev = var_filtered
        var_filtered = var_filtered_prev[
            ~var_filtered_prev['CONSEQUENCE'].str.contains('^upstream_gene', na=False) &
            ~var_filtered_prev['CONSEQUENCE'].str.contains('^downstream_gene', na=False)
        ]
        logger.info(
            f'Number of upstream_gene/downstream_gene variants excluded: '
            f'{len(var_filtered_prev) - len(var_filtered)}'
        )

    # Rename original file and write the filtered set
    check_subprocess(logger, f'mv {output_pass_tsv_gz} {output_pass_raw_tsv_gz}', debug)
    var_filtered = clean_annotations(var_filtered, yaml_data, logger=logger)
    pd_to_csv(df=var_filtered, fname=output_pass_tsv_gz, col_sep='\t')
    logger.info(f'Number of variants in final output TSV: {len(var_filtered)}')
