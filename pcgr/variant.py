#!/usr/bin/env python

import os
import re
import pandas as pd
import numpy as np
import warnings

from pcgr import pcgr_vars
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def set_genotype(variant_set: pd.DataFrame, logger) -> pd.DataFrame:
    """
    Set verbose genotype (homozygous, heterozygous) for each variant
    """
    variant_set['GENOTYPE'] = '.'    
    if {'GT'}.issubset(variant_set.columns):
        logger.info("Assignment of genotype (homozygous, heterozygous) for each variant based on 'GT' tag")
        heterozygous_states = []
        ref_allele_index = 0
        while ref_allele_index < 20:
            alt_allele_index = ref_allele_index + 1
            while alt_allele_index <= 20:
                phased_gt_1 = str(ref_allele_index) + "|" + str(alt_allele_index)
                phased_gt_2 = str(alt_allele_index) + "|" + str(ref_allele_index)
                unphased_gt_1 = str(ref_allele_index) + "/" + str(alt_allele_index)
                unphased_gt_2 = str(alt_allele_index) + "/" + str(ref_allele_index)
                gts = [phased_gt_1, phased_gt_2, unphased_gt_1, unphased_gt_2]
                heterozygous_states.extend(gts)
                
                alt_allele_index = alt_allele_index + 1
            ref_allele_index = ref_allele_index + 1
            
        homozygous_states = []    
        hom_allele_index = 1
        while hom_allele_index <= 20:
            phased_gt = str(hom_allele_index) + "|" + str(hom_allele_index)
            unphased_gt = str(hom_allele_index) + "/" + str(hom_allele_index)
            homozygous_states.extend([phased_gt, unphased_gt])
            hom_allele_index = hom_allele_index + 1
        
        variant_set.loc[variant_set['GT'].isin(homozygous_states), "GENOTYPE"] = "homozygous"
        variant_set.loc[variant_set['GT'].isin(heterozygous_states),"GENOTYPE"] = "heterozygous"
    else:
        variant_set['GENOTYPE'] = "undefined"
        
    return(variant_set)

def append_annotations(vcf2tsv_gz_fname: str, pcgr_db_dir: str, logger):
    
    clinvar_tsv_fname = os.path.join(pcgr_db_dir, 'variant','tsv','clinvar', 'clinvar.tsv.gz')
    protein_domain_tsv_fname = os.path.join(pcgr_db_dir, 'misc','tsv','protein_domain', 'protein_domain.tsv.gz') 
    gene_xref_tsv_fname = os.path.join(pcgr_db_dir, 'gene','tsv','gene_transcript_xref', 'gene_transcript_xref.tsv.gz')
    vcf2tsv_df = None
    clinvar_data_df = None
    
    num_recs_with_clinvar_hits = 0
    
    if os.path.exists(vcf2tsv_gz_fname):
        vcf2tsv_df = pd.read_csv(
            vcf2tsv_gz_fname, skiprows=[0], sep="\t", na_values=".",
            low_memory = False)
        if {'CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN','ENTREZGENE'}.issubset(vcf2tsv_df.columns):
            for elem in ['CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN','ENTREZGENE']:
                vcf2tsv_df = vcf2tsv_df.astype({elem:'string'})            
            vcf2tsv_df['CLINVAR_MSID'] = vcf2tsv_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df['PFAM_DOMAIN'] = vcf2tsv_df['PFAM_DOMAIN'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df["VAR_ID"] = vcf2tsv_df["CHROM"].str.cat(
                vcf2tsv_df["POS"], sep = "_").str.cat(
                    vcf2tsv_df["REF"], sep = "_").str.cat(
                        vcf2tsv_df["ALT"], sep = "_")

            if {'CLINVAR_TRAITS_ALL'}.issubset(vcf2tsv_df.columns):
                vcf2tsv_df.drop('CLINVAR_TRAITS_ALL', inplace=True, axis=1)
        
            ## check number of variants with ClinVar ID's
            num_recs_with_clinvar_hits = vcf2tsv_df["CLINVAR_MSID"].notna().sum()
            ## check number of variants with PFAM ID's
            num_recs_with_pfam_hits = vcf2tsv_df["PFAM_DOMAIN"].notna().sum()
            ## check number of variants with Ensembl gene ID's
            num_recs_with_entrez_hits = vcf2tsv_df["ENTREZGENE"].notna().sum()
    
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
                    prot_domains_data_df = pd.read_csv(
                        protein_domain_tsv_fname, sep="\t", usecols=["pfam_id","pfam_name"]).drop_duplicates()
                    prot_domains_data_df.rename(columns = {'pfam_id':'PFAM_DOMAIN', 'pfam_name':'PFAM_DOMAIN_NAME'}, inplace = True)                                        
                    vcf2tsv_df = vcf2tsv_df.merge(prot_domains_data_df, left_on=["PFAM_DOMAIN"], right_on=["PFAM_DOMAIN"], how="left")
                    vcf2tsv_df = vcf2tsv_df.fillna('.')
                else:
                    logger.error(f"Could not find {protein_domain_tsv_fname} needed for PFAM domain annotation - exiting")
            else:
                vcf2tsv_df['PFAM_DOMAIN_NAME'] = '.'
                vcf2tsv_df = vcf2tsv_df.fillna('.')
            
            if num_recs_with_entrez_hits > 0:
                
                if {'GENENAME'}.issubset(vcf2tsv_df.columns):
                    vcf2tsv_df.drop('GENENAME', inplace=True, axis=1)
                
                if os.path.exists(gene_xref_tsv_fname):
                    gene_xref_df = pd.read_csv(
                        gene_xref_tsv_fname, sep="\t", na_values=".", 
                        usecols=["entrezgene","name"])
                    gene_xref_df = gene_xref_df[gene_xref_df['entrezgene'].notnull()].drop_duplicates()
                    gene_xref_df["entrezgene"] = gene_xref_df["entrezgene"].astype("int64").astype("string")
                    #print(gene_xref_df.head)
                    gene_xref_df.rename(columns = {'entrezgene':'ENTREZGENE', 'name':'GENENAME'}, inplace = True)
                    vcf2tsv_df = vcf2tsv_df.merge(gene_xref_df, left_on=["ENTREZGENE"], right_on=["ENTREZGENE"], how="left")
                    vcf2tsv_df["ENTREZGENE"] = vcf2tsv_df['ENTREZGENE'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                    vcf2tsv_df = vcf2tsv_df.fillna('.')
                else:
                    logger.error(f"Could not find {gene_xref_tsv_fname} needed for gene name annotation - exiting")
            else:
                vcf2tsv_df['GENENAME'] = '.'
                vcf2tsv_df = vcf2tsv_df.fillna('.')
    
    return(vcf2tsv_df)

def set_allelic_support(variant_set: pd.DataFrame, allelic_support_tags: dict) -> pd.DataFrame:
    """
    Set allelic support for variants
    """
    for e in ['DP_CONTROL','DP_TUMOR','AF_CONTROL','AF_TUMOR']:
        if e in variant_set.columns:
            variant_set[e] = np.nan

    if allelic_support_tags['control_dp_tag'] != "_NA_":
        if {allelic_support_tags['control_dp_tag'],'DP_CONTROL'}.issubset(variant_set.columns):
            variant_set['DP_CONTROL'] = variant_set[allelic_support_tags['control_dp_tag']].astype(int)
    
    if allelic_support_tags['tumor_dp_tag'] != "_NA_":
        if {allelic_support_tags['tumor_dp_tag'],'DP_TUMOR'}.issubset(variant_set.columns):
            variant_set['DP_TUMOR'] = variant_set[allelic_support_tags['tumor_dp_tag']].astype(int)
    
    if allelic_support_tags['control_af_tag'] != "_NA_":
        if {allelic_support_tags['control_af_tag'],'AF_CONTROL'}.issubset(variant_set.columns):
            variant_set['AF_CONTROL'] = variant_set[allelic_support_tags['control_af_tag']].astype(float).round(4)
    
    if allelic_support_tags['tumor_af_tag'] != "_NA_":
        if {allelic_support_tags['tumor_af_tag'],'AF_TUMOR'}.issubset(variant_set.columns):
            variant_set['AF_TUMOR'] = variant_set[allelic_support_tags['tumor_af_tag']].astype(float).round(4)
    
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
    
    if {'CONSEQUENCE','AF_TUMOR','DP_TUMOR','VARIANT_CLASS'}.issubset(variant_set.columns):
        tmb_data_set = \
            variant_set[['CONSEQUENCE','DP_TUMOR','AF_TUMOR','VARIANT_CLASS']]
            
        n_rows_raw = len(tmb_data_set)
        if float(tumor_af_min) > 0:
            ## filter variant set to those with AF_TUMOR > af_min
            if variant_set[variant_set['AF_TUMOR'].isna() == True].empty is True:
                tmb_data_set = tmb_data_set.loc[tmb_data_set['AF_TUMOR'] >= float(tumor_af_min),:]
                n_removed_af = n_rows_raw - len(tmb_data_set)
                logger.info(f'Removing n = {n_removed_af} variants with allelic fraction (tumor sample) < {tumor_af_min}')
            else:
                logger.info(f"Cannot filter on sequencing depth - 'AF_TUMOR' contains missing values")
        if int(tumor_dp_min) > 0:
            ## filter variant set to those with AF_TUMOR > af_min
            if variant_set[variant_set['DP_TUMOR'].isna() == True].empty is True:
                tmb_data_set = tmb_data_set.loc[tmb_data_set['DP_TUMOR'] >= int(tumor_dp_min),:]
                n_removed_dp = n_rows_raw - len(tmb_data_set)
                logger.info(f"Removing n = {n_removed_dp} variants with sequencing depth (tumor sample) < {tumor_af_min}")
            else:
                logger.info(f"Cannot filter on sequencing depth - 'DP_TUMOR' contains missing values")
        
        if tmb_data_set.empty is False:
            counts_missense_only = len(tmb_data_set.loc[tmb_data_set['CONSEQUENCE'].str.match(pcgr_vars.CSQ_MISSENSE_PATTERN),:])
            counts_coding_and_silent = len(tmb_data_set.loc[tmb_data_set['CONSEQUENCE'].str.match(pcgr_vars.CSQ_CODING_SILENT_PATTERN),:])
            counts_coding_non_silent = len(tmb_data_set.loc[tmb_data_set['CONSEQUENCE'].str.match(pcgr_vars.CSQ_CODING_PATTERN),:])
            
            tmb_count_data = {}
            tmb_count_data['missense_only'] = counts_missense_only
            tmb_count_data['coding_and_silent'] = counts_coding_and_silent
            tmb_count_data['coding_non_silent'] = counts_coding_non_silent
            
            regex_data = {}
            regex_data['missense_only'] = pcgr_vars.CSQ_MISSENSE_PATTERN
            regex_data['coding_and_silent'] = pcgr_vars.CSQ_CODING_SILENT_PATTERN
            regex_data['coding_non_silent'] = pcgr_vars.CSQ_CODING_PATTERN
            
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

def clean_annotations(variant_set: pd.DataFrame, yaml_data: dict, germline: bool, logger) -> pd.DataFrame:
    
    if {'Consequence','EFFECT_PREDICTIONS','CLINVAR_CONFLICTED'}.issubset(variant_set.columns):
        variant_set.rename(columns = {'Consequence':'CONSEQUENCE'}, inplace = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("\\.&|\\.$", "NA&", regex = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("&$", "", regex = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("&", ", ", regex = True)
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] == 1, "CLINVAR_CONFLICTED"] = True
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] != 1, "CLINVAR_CONFLICTED"] = False
    
    if not {'VCF_SAMPLE_ID'}.issubset(variant_set.columns):
        variant_set['VCF_SAMPLE_ID'] = yaml_data['sample_id'].astype(str)
    variant_set['SAMPLE_ID'] = str(yaml_data['sample_id'])
    variant_set['GENOME_VERSION'] = yaml_data['genome_assembly']
    if {'CHROM','POS','REF','ALT',}.issubset(variant_set.columns):      
        variant_set['GENOMIC_CHANGE'] = variant_set['CHROM'].astype(str) + ":g." + variant_set['POS'].astype(str) + \
            variant_set['REF'].astype(str) + ">" + variant_set['ALT'].astype(str)
    if not {'PANEL_OF_NORMALS'}.issubset(variant_set.columns):
        variant_set['PANEL_OF_NORMALS'] = False
        
    ## Make sure that specific tags are formatted as integers (not float) during to_csv export
    if {'AMINO_ACID_END','AMINO_ACID_START'}.issubset(variant_set.columns):
        variant_set.loc[variant_set['AMINO_ACID_START'] == ".","AMINO_ACID_START"] = -1
        variant_set.loc[variant_set['AMINO_ACID_END'] == ".","AMINO_ACID_END"] = -1
        variant_set['AMINO_ACID_END'] = variant_set['AMINO_ACID_END'].astype(float).astype(int)
        variant_set['AMINO_ACID_START'] = variant_set['AMINO_ACID_START'].astype(float).astype(int)
    
    for pop in ['AFR','AMR','ASJ','FIN','EAS','NFE','SAS','OTH','']:
        for tag in ['AN','AC','NHOMALT']:
            vcf_info_tag = 'gnomADe_non_cancer_' + str(pop) + '_' + str(tag)
            if vcf_info_tag in variant_set.columns:
                variant_set[vcf_info_tag] = variant_set[vcf_info_tag].astype(str)
                variant_set.loc[variant_set[vcf_info_tag] != ".", vcf_info_tag] = \
                    variant_set.loc[variant_set[vcf_info_tag] != ".", vcf_info_tag].astype(float).astype(int)
    
    for elem in ['NUM_SUBMITTERS','ALLELE_ID','ENTREZGENE','REVIEW_STATUS_STARS','MSID']:
        vcf_info_tag = 'CLINVAR_' + str(elem)
        if vcf_info_tag in variant_set.columns:
            variant_set[vcf_info_tag] = variant_set[vcf_info_tag].astype(str)            
            variant_set.loc[variant_set[vcf_info_tag] != ".", vcf_info_tag] = \
                variant_set.loc[variant_set[vcf_info_tag] != ".", vcf_info_tag].astype(float).astype(int)
    
    for vcf_info_tag in ['ONCOGENE_RANK','TSG_RANK','TCGA_PANCANCER_COUNT','CGC_TIER','DISTANCE',
                         'EXON_AFFECTED','INTRON_POSITION','EXON_POSITION']:
        if vcf_info_tag in variant_set.columns:
            variant_set[vcf_info_tag] = variant_set[vcf_info_tag].astype(str)
            variant_set.loc[variant_set[vcf_info_tag] != ".", vcf_info_tag] = \
                variant_set.loc[variant_set[vcf_info_tag] != ".", vcf_info_tag].astype(float).astype(int)
    
    if germline is True:
        variant_set = set_genotype(variant_set, logger)
    
    return variant_set