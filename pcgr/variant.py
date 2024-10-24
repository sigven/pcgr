#!/usr/bin/env python

import os
import re
import pandas as pd
import numpy as np
import warnings

from pcgr import pcgr_vars
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def append_annotations(vcf2tsv_gz_fname: str, refdata_assembly_dir: str, logger):
    
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
            vcf2tsv_df['EFFECT_PREDICTIONS'] = vcf2tsv_df['EFFECT_PREDICTIONS'].str.replace("\\.&|\\.$", "NA&", regex = True)
            vcf2tsv_df['EFFECT_PREDICTIONS'] = vcf2tsv_df['EFFECT_PREDICTIONS'].str.replace("&$", "", regex = True)
            vcf2tsv_df['EFFECT_PREDICTIONS'] = vcf2tsv_df['EFFECT_PREDICTIONS'].str.replace("&", ", ", regex = True)
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
    
            #print(str(num_recs_with_entrez_hits))
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
                    clinvar_data_df.loc[(clinvar_data_df['classification'] == "VUS") & (clinvar_data_df['conflicted'] == 1),"classification"] = \
                        "VUS/Conflicting evidence"
                    clinvar_data_df['CLINVAR_TRAITS_ALL'] = clinvar_data_df['classification'].str.cat(
                        clinvar_data_df['origin_simple'].str.capitalize().str.cat(
                        clinvar_data_df['trait'], sep = " - "), sep = " - ")
                    clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['variation_id']
                    clinvar_data_df = clinvar_data_df.astype({'CLINVAR_MSID':'string'})
                    clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                    clinvar_data_df = clinvar_data_df[['VAR_ID','CLINVAR_MSID','CLINVAR_TRAITS_ALL']]
                    
                    vcf2tsv_df = vcf2tsv_df.merge(
                        clinvar_data_df, left_on=["VAR_ID", "CLINVAR_MSID"], right_on=["VAR_ID", "CLINVAR_MSID"], how="left")
                else:
                    logger.error(f"Could not find {clinvar_tsv_fname} needed for ClinVar variant annotation - exiting")
            else:
                vcf2tsv_df['CLINVAR_TRAITS_ALL'] = '.'
                
            
            ## merge variant set with PFAM domain annotations
            if num_recs_with_pfam_hits > 0:
                
                if {'PFAM_DOMAIN_NAME'}.issubset(vcf2tsv_df.columns):
                    vcf2tsv_df.drop('PFAM_DOMAIN_NAME', inplace=True, axis=1)
                
                if os.path.exists(protein_domain_tsv_fname):
                    prot_domains_data_df = pd.read_csv(
                        protein_domain_tsv_fname, sep="\t", usecols=["pfam_id","pfam_name"]).drop_duplicates()
                    prot_domains_data_df.rename(columns = {'pfam_id':'PFAM_DOMAIN', 'pfam_name':'PFAM_DOMAIN_NAME'}, inplace = True)                                        
                    vcf2tsv_df = vcf2tsv_df.merge(prot_domains_data_df, left_on=["PFAM_DOMAIN"], right_on=["PFAM_DOMAIN"], how="left")
                else:
                    logger.error(f"Could not find {protein_domain_tsv_fname} needed for PFAM domain annotation - exiting")
            else:
                vcf2tsv_df['PFAM_DOMAIN_NAME'] = '.'
            
            if num_recs_with_entrez_hits > 0:
                
                if {'GENENAME'}.issubset(vcf2tsv_df.columns):
                    vcf2tsv_df.drop('GENENAME', inplace=True, axis=1)
                
                if os.path.exists(gene_xref_tsv_fname):
                    gene_xref_df = pd.read_csv(
                        gene_xref_tsv_fname, sep="\t", na_values=".", 
                        usecols=["entrezgene","name"])
                    gene_xref_df = gene_xref_df[gene_xref_df['entrezgene'].notnull()].drop_duplicates()
                    gene_xref_df.rename(columns = {'entrezgene':'ENTREZGENE', 'name':'GENENAME'}, inplace = True)                
                    vcf2tsv_df = vcf2tsv_df.merge(gene_xref_df, left_on=["ENTREZGENE"], right_on=["ENTREZGENE"], how="left")
                else:
                    logger.error(f"Could not find {gene_xref_tsv_fname} needed for gene name annotation - exiting")
            else:
                vcf2tsv_df['GENENAME'] = '.'
                
    return(vcf2tsv_df)

def set_allelic_support(variant_set: pd.DataFrame, allelic_support_tags: dict) -> pd.DataFrame:
    """
    Set allelic support for variants
    """
    for e in ['DP_CONTROL','DP_TUMOR','VAF_CONTROL','VAF_TUMOR']:
        if e in variant_set.columns:
            variant_set[e] = np.nan        

    if allelic_support_tags['control_dp_tag'] != "_NA_":
        if {allelic_support_tags['control_dp_tag'],'DP_CONTROL'}.issubset(variant_set.columns):
            variant_set['DP_CONTROL'] = variant_set[allelic_support_tags['control_dp_tag']].astype(int)
    
    if allelic_support_tags['tumor_dp_tag'] != "_NA_":
        if {allelic_support_tags['tumor_dp_tag'],'DP_TUMOR'}.issubset(variant_set.columns):
            variant_set['DP_TUMOR'] = variant_set[allelic_support_tags['tumor_dp_tag']].astype(int)
    
    if allelic_support_tags['control_af_tag'] != "_NA_":
        if {allelic_support_tags['control_af_tag'],'VAF_CONTROL'}.issubset(variant_set.columns):
            variant_set['VAF_CONTROL'] = variant_set[allelic_support_tags['control_af_tag']].astype(float).round(4)
    
    if allelic_support_tags['tumor_af_tag'] != "_NA_":
        if {allelic_support_tags['tumor_af_tag'],'VAF_TUMOR'}.issubset(variant_set.columns):
            variant_set['VAF_TUMOR'] = variant_set[allelic_support_tags['tumor_af_tag']].astype(float).round(4)
    
    return variant_set

def calculate_tmb(variant_set: pd.DataFrame, tumor_dp_min: int, tumor_af_min: float, 
                  target_size_mb: float, sample_id: str, tmb_fname: str, logger) -> int:
    """
    Calculate TMB for somatic variant set
    """
    logger.info(f"Calculating tumor mutational load/TMB from somatic call set - effective coding target size (Mb): {target_size_mb}")
    
    counts_missense_only = 0
    counts_coding_and_silent = 0
    counts_coding_non_silent = 0
    
    if {'CONSEQUENCE','VAF_TUMOR','DP_TUMOR','VARIANT_CLASS'}.issubset(variant_set.columns):
        tmb_data_set = \
            variant_set[['CONSEQUENCE','DP_TUMOR','VAF_TUMOR','VARIANT_CLASS']]
        
        n_rows_raw = len(tmb_data_set)
        if float(tumor_af_min) > 0:
            ## filter variant set to those with VAF_TUMOR > tumor_af_min
            if variant_set[variant_set['VAF_TUMOR'].isna() == True].empty is True:
                tmb_data_set = tmb_data_set.loc[tmb_data_set['VAF_TUMOR'].astype(float) >= float(tumor_af_min),:]
                n_removed_af = n_rows_raw - len(tmb_data_set)
                logger.info(f'Removing n = {n_removed_af} variants with allelic fraction (tumor sample) < {tumor_af_min}')
            else:
                logger.info(f"Cannot filter on sequencing depth - 'VAF_TUMOR' contains missing values")
        if int(tumor_dp_min) > 0:
            ## filter variant set to those with DP_TUMOR > tumor_dp_min
            if variant_set[variant_set['DP_TUMOR'].isna() == True].empty is True:
                tmb_data_set = tmb_data_set.loc[tmb_data_set['DP_TUMOR'].astype(int) >= int(tumor_dp_min),:]
                n_removed_dp = n_rows_raw - len(tmb_data_set)
                logger.info(f"Removing n = {n_removed_dp} variants with sequencing depth (tumor sample) < {tumor_dp_min}")
            else:
                logger.info(f"Cannot filter on sequencing depth - 'DP_TUMOR' contains missing values")
        
        regex_data = {}
        regex_data['missense_only'] = pcgr_vars.CSQ_MISSENSE_PATTERN
        regex_data['coding_and_silent'] = pcgr_vars.CSQ_CODING_SILENT_PATTERN
        regex_data['coding_non_silent'] = pcgr_vars.CSQ_CODING_PATTERN
        
        tmb_count_data = {}
        tmb_count_data['missense_only'] = 0
        tmb_count_data['coding_and_silent'] = 0
        tmb_count_data['coding_non_silent'] = 0
        
        tmb_estimates = {}
        tmb_estimates['missense_only'] = 0
        tmb_estimates['coding_and_silent'] = 0
        tmb_estimates['coding_non_silent'] = 0
        
        if tmb_data_set.empty is False:
            counts_missense_only = len(tmb_data_set.loc[tmb_data_set['CONSEQUENCE'].str.match(pcgr_vars.CSQ_MISSENSE_PATTERN),:])
            counts_coding_and_silent = len(tmb_data_set.loc[tmb_data_set['CONSEQUENCE'].str.match(pcgr_vars.CSQ_CODING_SILENT_PATTERN),:])
            counts_coding_non_silent = len(tmb_data_set.loc[tmb_data_set['CONSEQUENCE'].str.match(pcgr_vars.CSQ_CODING_PATTERN),:])
            
            tmb_count_data = {}
            tmb_count_data['missense_only'] = counts_missense_only
            tmb_count_data['coding_and_silent'] = counts_coding_and_silent
            tmb_count_data['coding_non_silent'] = counts_coding_non_silent
            
            tmb_estimates = {}
            tmb_estimates['missense_only'] = round(float(counts_missense_only / target_size_mb), 4)
            tmb_estimates['coding_and_silent'] = round(float(counts_coding_and_silent / target_size_mb), 4)
            tmb_estimates['coding_non_silent'] = round(float(counts_coding_non_silent / target_size_mb), 4)
        
            
        tmb_recs = [[sample_id, n_rows_raw, 'TMB_missense_only', regex_data['missense_only'], target_size_mb, 
                    tumor_dp_min, tumor_af_min, tmb_count_data['missense_only'], tmb_estimates['missense_only'], 'mutations/Mb'],
                    [sample_id, n_rows_raw, 'TMB_coding_and_silent', regex_data['coding_and_silent'], target_size_mb, 
                    tumor_dp_min, tumor_af_min, tmb_count_data['coding_and_silent'], tmb_estimates['coding_and_silent'], 'mutations/Mb'],
                    [sample_id, n_rows_raw, 'TMB_coding_non_silent', regex_data['coding_non_silent'], target_size_mb, 
                    tumor_dp_min, tumor_af_min, tmb_count_data['coding_non_silent'], tmb_estimates['coding_non_silent'], 'mutations/Mb']]
            
        df = pd.DataFrame(
            tmb_recs, 
            columns=['sample_id', 
                    'n_somatic_variants', 
                    'tmb_measure', 
                    'tmb_csq_regex',
                    'tmb_target_size_mb', 
                    'tmb_dp_min', 
                    'tmb_af_min',
                    'tmb_n_variants', 
                    'tmb_estimate',
                    'tmb_unit'])
        df.to_csv(tmb_fname, sep="\t", header=True, index=False)

def clean_annotations(variant_set: pd.DataFrame, yaml_data: dict, logger) -> pd.DataFrame:
    
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
    
    for pop in ['AFR','AMR','ASJ','FIN','EAS','NFE','SAS','OTH','']:
        for tag in ['AN','AC','NHOMALT']:
            vcf_info_tag = 'gnomADe_non_cancer_' + str(pop) + '_' + str(tag)
            if pop == '':
                vcf_info_tag = 'gnomADe_non_cancer_' + str(tag)
            if vcf_info_tag in variant_set.columns:
                variant_set.loc[variant_set[vcf_info_tag].isna(),vcf_info_tag] = -123456789           
                variant_set[vcf_info_tag] = variant_set[vcf_info_tag].apply(lambda x: str(int(x)) )
                variant_set.loc[variant_set[vcf_info_tag] == "-123456789", vcf_info_tag] = np.nan
                
    ## Make sure that specific tags are formatted as integers (not float) when exported to TSV
    ## Trick - convert to str and then back to int
    for vcf_info_tag in ['ONCOGENE_RANK',
                         'TSG_RANK',
                         'TCGA_PANCANCER_COUNT',
                         'CGC_TIER',
                         'ENTREZGENE',
                         'DP_TUMOR',
                         'DP_CONTROL',
                         'DISTANCE',
                         'EXON_AFFECTED',
                         'INTRON_POSITION',
                         'EXON_POSITION',
                         'AMINO_ACID_END',
                         'AMINO_ACID_START',
                         'ONCOGENICITY_SCORE',
                         'CLINVAR_NUM_SUBMITTERS',
                         'CLINVAR_ALLELE_ID',
                         'CLINVAR_ENTREZGENE',
                         'CLINVAR_REVIEW_STATUS_STARS']:
        if vcf_info_tag in variant_set.columns: 
            variant_set.loc[variant_set[vcf_info_tag].isna(),vcf_info_tag] = -123456789           
            variant_set[vcf_info_tag] = variant_set[vcf_info_tag].apply(lambda x: str(int(x)) )
            variant_set.loc[variant_set[vcf_info_tag] == "-123456789", vcf_info_tag] = np.nan
    
    return variant_set

def reverse_complement_dna(dna_string = "C"):
    pairs = {
        "A":"T",
        "C":"G",
        "G":"C",
        "T":"A",
    }
    reverse_complement = ""
    i = len(dna_string) - 1
    while i >= 0:
        base = str(dna_string[i]).upper()
        if base in pairs:
            complement = pairs[base]
        else:
            complement = base
        reverse_complement += complement
        i = i - 1
    return reverse_complement     
    