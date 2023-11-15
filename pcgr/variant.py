#!/usr/bin/env python

import os
import re
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

def append_annotations(vcf2tsv_gz_fname: str, pcgr_db_dir: str, logger):
    
    clinvar_tsv_fname = os.path.join(pcgr_db_dir, 'variant','tsv','clinvar', 'clinvar.tsv.gz')
    protein_domain_tsv_fname = os.path.join(pcgr_db_dir, 'misc','tsv','protein_domain', 'protein_domain.tsv.gz')    
    vcf2tsv_df = None
    clinvar_data_df = None
    
    num_recs_with_clinvar_hits = 0
    
    if os.path.exists(vcf2tsv_gz_fname):
        vcf2tsv_df = pd.read_csv(
            vcf2tsv_gz_fname, skiprows=[0], sep="\t", na_values=".",
            low_memory = False)
        if {'CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN'}.issubset(vcf2tsv_df.columns):
            for elem in ['CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN']:
                vcf2tsv_df = vcf2tsv_df.astype({elem:'string'})
            vcf2tsv_df['CLINVAR_MSID'] = vcf2tsv_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df['PFAM_DOMAIN'] = vcf2tsv_df['PFAM_DOMAIN'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df["VAR_ID"] = vcf2tsv_df["CHROM"].str.cat(
                vcf2tsv_df["POS"], sep = "_").str.cat(
                    vcf2tsv_df["REF"], sep = "_").str.cat(
                        vcf2tsv_df["ALT"], sep = "_")

            vcf2tsv_df.drop('CLINVAR_TRAITS_ALL', inplace=True, axis=1)
        
        ## check number of variants with ClinVar ID's
        num_recs_with_clinvar_hits = vcf2tsv_df["CLINVAR_MSID"].notna().sum()
        ## check number of variants with PFAM ID's
        num_recs_with_pfam_hits = vcf2tsv_df["PFAM_DOMAIN"].notna().sum()
    
        ## merge variant set with ClinVar trait and variant origin annotations
        if num_recs_with_clinvar_hits > 0:
            if os.path.exists(clinvar_tsv_fname):
                clinvar_data_df = pd.read_csv(
                    clinvar_tsv_fname, sep="\t", 
                    usecols=["variation_id","origin_simple","VAR_ID","trait"],
                    low_memory = False)
                clinvar_data_df['CLINVAR_TRAITS_ALL'] = clinvar_data_df['origin_simple'].str.capitalize().str.cat(
                    clinvar_data_df['trait'], sep = " - ")
                clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['variation_id']
                clinvar_data_df = clinvar_data_df.astype({'CLINVAR_MSID':'string'})
                clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                clinvar_data_df = clinvar_data_df[['VAR_ID','CLINVAR_MSID','CLINVAR_TRAITS_ALL']]
                
                vcf2tsv_df = vcf2tsv_df.merge(
                    clinvar_data_df, left_on=["VAR_ID", "CLINVAR_MSID"], right_on=["VAR_ID", "CLINVAR_MSID"], how="left")
                vcf2tsv_df = vcf2tsv_df.fillna('.')
            else:
                logger.error(f"Could not find {clinvar_tsv_fname} needed for ClinVar variant annotation - exiting")
        else:
            vcf2tsv_df['CLINVAR_TRAITS_ALL'] = '.'
            vcf2tsv_df = vcf2tsv_df.fillna('.')
            
        
        ## merge variant set with PFAM domain annotations
        if num_recs_with_pfam_hits > 0:
            
            vcf2tsv_df.drop('PFAM_DOMAIN_NAME', inplace=True, axis=1)
            
            if os.path.exists(protein_domain_tsv_fname):
                prot_domains_data_df = pd.read_csv(protein_domain_tsv_fname, sep="\t", usecols=["pfam_id","pfam_name"])
                prot_domains_data_df['PFAM_DOMAIN'] = prot_domains_data_df['pfam_id']
                prot_domains_data_df['PFAM_DOMAIN_NAME'] = prot_domains_data_df['pfam_name']
                prot_domains_data_df = prot_domains_data_df[['PFAM_DOMAIN','PFAM_DOMAIN_NAME']]
                
                vcf2tsv_df = vcf2tsv_df.merge(prot_domains_data_df, left_on=["PFAM_DOMAIN"], right_on=["PFAM_DOMAIN"], how="left")
                vcf2tsv_df = vcf2tsv_df.fillna('.')
            else:
                logger.error(f"Could not find {protein_domain_tsv_fname} needed for PFAM domain annotation - exiting")
        else:
            vcf2tsv_df['PFAM_DOMAIN_NAME'] = '.'
            vcf2tsv_df = vcf2tsv_df.fillna('.')
    
    return(vcf2tsv_df)

def set_allelic_support(variant_set: pd.DataFrame, allelic_support_tags: dict) -> pd.DataFrame:
    """
    Set allelic support for variants
    """
    
    if allelic_support_tags['control_dp_tag'] != "_NA_":
        if {allelic_support_tags['control_dp_tag'],'DP_CONTROL'}.issubset(variant_set.columns):
            variant_set['DP_CONTROL'] = variant_set[allelic_support_tags['control_dp_tag']].astype(int)
    
    if allelic_support_tags['tumor_dp_tag'] != "_NA_":
        if {allelic_support_tags['tumor_dp_tag'],'DP_TUMOR'}.issubset(variant_set.columns):
            variant_set['DP_TUMOR'] = variant_set[allelic_support_tags['tumor_dp_tag']].astype(int)
    
    if allelic_support_tags['control_af_tag'] != "_NA_":
        if {allelic_support_tags['control_af_tag'],'AF_CONTROL'}.issubset(variant_set.columns):
            variant_set['AF_CONTROL'] = variant_set[allelic_support_tags['control_af_tag']].astype(float).round(3)
    
    if allelic_support_tags['tumor_af_tag'] != "_NA_":
        if {allelic_support_tags['tumor_af_tag'],'AF_TUMOR'}.issubset(variant_set.columns):
            variant_set['AF_TUMOR'] = variant_set[allelic_support_tags['tumor_af_tag']].astype(float).round(3)
    
    return variant_set

def calculate_tmb(variant_set: pd.DataFrame, tmb_fname: str) -> int:
    """
    Calculate TMB for somatic variant set
    """
    
    if tmb_fname != "_NA_":
        if {tmb_fname,'TMB'}.issubset(variant_set.columns):
            variant_set['TMB'] = variant_set[tmb_fname].astype(float)
    
    
    
    