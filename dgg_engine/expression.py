#!/usr/bin/env python

import os,re
import csv
import gzip
import pandas as pd
import logging
import numpy as np

from dgg_engine import utils
from dgg_engine import pcgr_vars
from dgg_engine.utils import error_message, warn_message, check_file_exists, remove_file

def parse_expression(expression_fname_tsv: str,
                     sample_id: str,
                     refdata_assembly_dir: str,                                                                    
                     logger: logging.Logger = None):
    """
    Parse the expression data from the user-provided file and verify record identifiers with 
    reference data. Aggregate transcript-level TPM values to generate gene-level TPM values.
    Args:
        expression_fname_tsv (str): Path to the user-provided gene expression file.
        refdata_assembly_dir (str): Path to the build-specific PCGR reference data directory.
        logger (logging.Logger, optional): Logger. Defaults to None.
    """    
    ## Check that the expression file exists
    check_file_exists(expression_fname_tsv, logger)
    ## Read the expression file
    sample_gene_expression = pd.read_csv(expression_fname_tsv, sep = "\t", na_values = ".")
    if not {'TargetID','TPM'}.issubset(sample_gene_expression.columns):
        err_msg = f"Could not find required columns ('TargetID','TPM') in input gene expression file: {expression_fname_tsv} - exiting."
        error_message(err_msg, logger)
    sample_gene_expression  = sample_gene_expression.rename(columns = {'TargetID':'ID'})
    
    ## Remove version number from transcript identifiers
    sample_gene_expression['ID'] = sample_gene_expression['ID'].str.replace(r'(\\.[0-9]{1,})$', '', regex=True)
    sample_gene_expression.sort_values(by=['ID','TPM'], ascending=[True, False], inplace = True)
    dup_ids = len(sample_gene_expression['ID']) - len(sample_gene_expression['ID'].drop_duplicates())
    if dup_ids > 0: 
        logger.warning(f"Found N = {dup_ids} duplicate identifiers - resolving duplicates by keeping the highest TPM value")
        sample_gene_expression = sample_gene_expression.drop_duplicates(subset = ['ID']) 
    
    ## Read the gene identifier index - maps transcript identifiers (Ensembl/Refseq), 
    ## gene symbols/aliases and identifiers (Ensembl) to Entrez gene identifers and
    ## biotype
    gene_index_fname_tsv =  \
        os.path.join(refdata_assembly_dir, "gene", "tsv", "gene_transcript_xref", "gene_index.tsv.gz")
    trans_index_fname_tsv =  \
        os.path.join(refdata_assembly_dir, "gene", "tsv", "gene_transcript_xref", "gene_transcript_xref.tsv.gz")
    check_file_exists(gene_index_fname_tsv, logger)
    gene_index = pd.read_csv(gene_index_fname_tsv, sep = "\t", na_values = ".",  low_memory = False)
    
    sample_identifiers_found = 0
    expression_map = {}
    expression_map['transcript'] = None
    expression_map['gene'] = None
    expression_map['transcript_identifier'] = None
    
    if {'ID','ID_TYPE','ENTREZGENE','ENSEMBL_GENE_ID','SYMBOL','BIOTYPE','AMBIGUOUS_ID'}.issubset(gene_index.columns):
        exp_map = pd.merge(sample_gene_expression, gene_index, how="left", on=["ID"])
        exp_map.loc[:,'ENSEMBL_GENE_ID'] = exp_map.loc[:,'ENSEMBL_GENE_ID'].str.split('|')
        exp_map = exp_map.explode('ENSEMBL_GENE_ID').drop_duplicates().reset_index(drop = True)
        
        if {'ID_TYPE'}.issubset(exp_map.columns):
            ## Filter out records that are not properly verified            
            mask = pd.notna(exp_map['ID_TYPE'])
            exp_map_verified = exp_map[mask]
            
            ## Count number of identifier types present in the input file
            keytype_stats = exp_map_verified['ID_TYPE'].value_counts().reset_index().set_axis(['key','count'], axis=1)
            keytype_stats['pct'] = round((keytype_stats['count'] / sum(keytype_stats['count'])) * 100, 2)
            identifiers_present = {}
            for index, row in keytype_stats.iterrows():
                identifiers_present[row['key']] = row['pct']
                
            if len(identifiers_present.keys()) > 1:
                err_msg = f"More than one identifier type present in input expression file: {expression_fname_tsv} - exiting."
                error_message(err_msg, logger)
            
            identifiers_used_in_input = "Ensembl gene IDs"
            transcript_level = False
            transcript_identifier = None
            if "refseq_transcript_id" in identifiers_present:
                identifiers_used_in_input = "RefSeq transcript IDs"
                transcript_level = True
                transcript_identifier = 'RefSeq'
                
                 ## Create mapping between Refseq transcript identifiers and Ensembl transcript identifiers
                check_file_exists(trans_index_fname_tsv, logger)
                refseq2ensembl = pd.read_csv(trans_index_fname_tsv, sep = "\t", na_values = ".",  low_memory = False)
                if {'ensembl_transcript_id','refseq_transcript_id'}.issubset(refseq2ensembl.columns):
                    refseq2ensembl = refseq2ensembl[['ensembl_transcript_id','refseq_transcript_id']].drop_duplicates().rename(columns=str.upper)
                    refseq2ensembl = refseq2ensembl.loc[~refseq2ensembl['REFSEQ_TRANSCRIPT_ID'].isna(), :]
                    refseq2ensembl.loc[:,'REFSEQ_TRANSCRIPT_ID'] = refseq2ensembl.loc[:,'REFSEQ_TRANSCRIPT_ID'].str.split('&')
                    refseq2ensembl = refseq2ensembl.explode('REFSEQ_TRANSCRIPT_ID').drop_duplicates().reset_index(drop = True)
                    refseq2ensembl = refseq2ensembl.rename(columns = {'REFSEQ_TRANSCRIPT_ID':'ID'}).drop_duplicates(subset = ['ID'])
                            
                    ## Only keep RefSeq transcript identifiers that are mapped to Ensembl transcript identifiers
                    exp_map_verified = pd.merge(exp_map_verified, refseq2ensembl, how="left", on=["ID"])
                    exp_map_verified = exp_map_verified.loc[~exp_map_verified['ENSEMBL_TRANSCRIPT_ID'].isna(), :]
                
            if "ensembl_transcript_id" in identifiers_present:
                identifiers_used_in_input = "Ensembl transcript IDs"
                transcript_level = True
                transcript_identifier = 'Ensembl'
                exp_map_verified = exp_map_verified.assign(ENSEMBL_TRANSCRIPT_ID = exp_map_verified.ID)
                
            if "symbol_alias" in identifiers_present:
                identifiers_used_in_input = "gene symbols"
            
            ## Emit warning if more than 5% of gene/transcript identifiers are not properly verified
            sample_identifiers_found = len(exp_map_verified)
            percent_verified = round((len(exp_map_verified) / len(exp_map)) * 100, 2)
            percent_missing = round(100 - percent_verified, 2)
            if percent_missing > 5:
                logger.warning("Failed to map " + str(percent_missing) +  \
                    "% of gene/transcript identifiers in input TSV file - use proper ENST/RefSeq identifiers")
            logger.info("Verified N = " + str(sample_identifiers_found) + " (" + str(percent_verified) + \
                "%) of gene/transcript identifiers in input gene expression file - using " + str(identifiers_used_in_input))
            
            ## Emit warning if ambiguous gene/transcript identifiers are detected - 
            ## remove them from the analysis (write them to a separate file?)
            n_ambig = len(exp_map_verified[exp_map_verified.AMBIGUOUS_ID == True])
            if n_ambig > 0:
                logger.warning("Detected N = " + str(n_ambig) + " ambiguous gene/transcript identifiers in input gene expression file")
            else:
                logger.info("NO ambiguous gene/transcript identifiers were detected in input gene expression file")
            transcript_expression_map = exp_map_verified[exp_map_verified.AMBIGUOUS_ID == False]
            exp_map_verified_nonzero = transcript_expression_map[transcript_expression_map.TPM > 0]
            
            ## inform the user about the statistics with respect to biotypes (TPM > 0)
            biotype_stats = exp_map_verified_nonzero['BIOTYPE'].value_counts().reset_index().set_axis(['key','count'], axis=1)
            biotype_stats['pct'] = round((biotype_stats['count'] / sum(biotype_stats['count'])) * 100, 2) 
            biotype_stats = biotype_stats.head(10)
            biotype_stats['type_stats_verbose'] = str("'") + biotype_stats['key'].astype(str) + "': " + \
                biotype_stats['pct'].astype(str) + "% (" + biotype_stats['count'].astype(str) + ")"
            logger.info("Top gene/transcript biotype categories (considering records with TPM > 0):")
            j = 0
            while j < len(biotype_stats):
                logger.info(biotype_stats['type_stats_verbose'][j])
                j = j + 1
            
            ## input is provided at the transcript level
            if transcript_level is True:
                expression_map['transcript_identifier'] = transcript_identifier
                expression_map['transcript'] = transcript_expression_map.sort_values(by='TPM', ascending=False)
                expression_map['transcript']['SAMPLE_ID'] = sample_id
                expression_map['transcript'] = \
                    expression_map['transcript'][["SAMPLE_ID","ID", "ID_TYPE","AMBIGUOUS_ID","TPM",
                                                  "ENSEMBL_TRANSCRIPT_ID","ENSEMBL_GENE_ID","SYMBOL",
                                                  "ENTREZGENE","GENENAME","BIOTYPE"]]
    
            ## for transcripts (i.e. provided through VEP) is not supported with transcript expression by
            ## user (missing TPM's), record the minimum TPM value across all other transcripts to be used as proxy
            expression_map['gene_min_tpm'] = transcript_expression_map.groupby(
            ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE']).agg({'TPM':'min'}).reset_index()
            expression_map['gene_min_tpm'].columns = ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE','TPM_MIN']
            expression_map['gene_min_tpm'] = expression_map['gene_min_tpm'].drop_duplicates().sort_values(by='TPM_MIN', ascending=False)
                
            ## make gene level TPM summary  
            expression_map['gene'] = transcript_expression_map.groupby(
                ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE']).agg({'TPM':'sum'}).reset_index()
            ## add log2(TPM + 0.001) for gene-level TPM values (all reference TPM values are in log2(TPM + 0.001))
            expression_map['gene']['TPM_LOG2_GENE'] = np.log2(expression_map['gene']['TPM'] + 0.001)
            expression_map['gene'].columns = ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE','TPM_GENE','TPM_LOG2_GENE']
            expression_map['gene'] = expression_map['gene'].drop_duplicates().sort_values(by='TPM_GENE', ascending=False)
    return(expression_map)


def integrate_variant_expression(variant_set: pd.DataFrame, 
                                 expression_data: dict,
                                 logger: logging.Logger = None) -> pd.DataFrame:
    """
    Integrates annotated variant calls (in variant_set dataframe) with expression data (in expression_data dictionary)
    Both transcript-level TPM values and gene-level TPM values are integrated into the variant set.
    Args:
        variant_set (pd.DataFrame): Pandas dataframe with annotated variants
        expression_data (dict): Dictionary with expression data (transcript and gene-level TPM values)
        logger (logging.Logger, optional): Logger. Defaults to None.
    """    
    if not expression_data is None:
        for s in ['gene','gene_min_tpm']:
            if s in expression_data.keys():
                if not expression_data[s] is None:
                    if s == 'gene':
                        logger.info("Integrating gene-level expression data from tumor into somatic variant set")
                    if expression_data[s].empty:
                        logger.warning('Expression file does not contain any gene-level expression data')
                    else:
                        v = 'TPM_MIN'
                        if s == 'gene':
                            v = 'TPM_GENE'                           
                        if v in variant_set.columns:
                            variant_set.drop(v, inplace=True, axis=1)
                        if {v,'ENSEMBL_GENE_ID'}.issubset(expression_data[s].columns) and \
                            'ENSEMBL_GENE_ID' in variant_set.columns:
                                exp_gene = expression_data[s][[v,'ENSEMBL_GENE_ID']]
                                variant_set = variant_set.merge(exp_gene, on = 'ENSEMBL_GENE_ID', how = 'left')
                    
        if 'transcript' in expression_data.keys():
            if not expression_data['transcript'] is None:
                logger.info("Integrating transcript-level expression data from tumor into somatic variant set")
                if expression_data['transcript'].empty:
                    logger.warning('Expression file does not contain any transcript-level expression data')
                    if 'TPM_GENE' in variant_set.columns:
                        variant_set = variant_set.assign(TPM = variant_set['TPM_GENE'])
                else:
                    exp_data_minimal = expression_data['transcript'][['ENSEMBL_TRANSCRIPT_ID','TPM']]
                    if {'VAR_ID', 'ENSEMBL_TRANSCRIPT_ID'}.issubset(variant_set.columns):
                        var2trans_all = variant_set[['VAR_ID','ENSEMBL_TRANSCRIPT_ID']]
                        mask = pd.notna(var2trans_all['ENSEMBL_TRANSCRIPT_ID'])
                        var2trans = var2trans_all[mask]
                        if not var2trans.empty:                               
                            var2trans = pd.merge(
                                var2trans, exp_data_minimal, on = 'ENSEMBL_TRANSCRIPT_ID', how = 'left')                                
                            var2exp = var2trans.groupby(['VAR_ID','ENSEMBL_TRANSCRIPT_ID']).agg({'TPM':'max'}).reset_index()
                            if 'TPM' in variant_set.columns:
                                variant_set.drop('TPM', inplace=True, axis=1)
                            if {'TPM','VAR_ID','ENSEMBL_TRANSCRIPT_ID'}.issubset(var2exp.columns):
                                var2exp = var2exp[['TPM','VAR_ID','ENSEMBL_TRANSCRIPT_ID']].drop_duplicates()
                                variant_set = variant_set.merge(var2exp, on = ['VAR_ID','ENSEMBL_TRANSCRIPT_ID'], how = 'left')
                                
                        else:
                            variant_set['TPM'] = np.nan
                            logger.warning('Variant file does not contain any entries with valid transcript identifiers')
            else:
                logger.warning('Expression file does not contain any transcript-level expression data')
                if 'TPM_GENE' in variant_set.columns:
                    variant_set = variant_set.assign(TPM = variant_set['TPM_GENE'])
    
    return(variant_set)

def aggregate_tpm_per_cons(variant_set: pd.DataFrame, 
                        expression_data: dict,
                        logger: logging.Logger = None) -> pd.DataFrame:
    
    """
    Integrate transcript-level expression data into somatic variant set
    - aggregate transcript-level expression data per consequence category
    """
 
    cons2exp= pd.DataFrame()
    
    if 'transcript' in expression_data.keys():
        if expression_data['transcript'].empty:
            logger.warning('Expression file does not contain any transcript-level expression data')
            return(cons2exp)
        if variant_set.empty:
            logger.warning('Variant file does not contain any entries with valid transcript identifiers')
            return(cons2exp)
        if {'VAR_ID','VEP_ALL_CSQ'}.issubset(variant_set.columns) and \
            {'ENSEMBL_TRANSCRIPT_ID','TPM'}.issubset(expression_data['transcript'].columns):
            varset = variant_set[['VAR_ID','VEP_ALL_CSQ']]
            
            ## Get all transcript-specific consequences provided by VEP - aggregate expression per gene consequence type
            trans_expression = expression_data['transcript'][['ENSEMBL_TRANSCRIPT_ID','TPM']]
            varset.loc[:,'VEP_ALL_CSQ'] = varset.loc[:,'VEP_ALL_CSQ'].str.split(',')
            varset = varset.explode('VEP_ALL_CSQ').drop_duplicates().reset_index(drop = True)
            varset[['consequence', 'SYMBOL','ENTREZGENE','HGVSC','HGVSP','EXON','FEATURE_TYPE','ENSEMBL_TRANSCRIPT_ID','BIOTYPE']] = \
                varset['VEP_ALL_CSQ'].str.split(':', expand=True)
            varset = varset[varset["ENSEMBL_TRANSCRIPT_ID"].str.contains("ENST")]
            if varset.empty:
                logger.warning('Variant file does not contain any entries with valid transcript identifiers')
                return(cons2exp)
            varset = pd.merge(varset, trans_expression, on = 'ENSEMBL_TRANSCRIPT_ID', how = 'left')
            varset = varset.loc[~varset['TPM'].isna(), :]
            if not varset.empty:
                pd.options.mode.chained_assignment = None
                varset.loc[:, "consTPM_SINGLE"] = varset['consequence'].str.replace('&.+','', regex=True)
                var2trans = varset[['VAR_ID','consTPM_SINGLE','ENTREZGENE','ENSEMBL_TRANSCRIPT_ID']]
                var2trans['ENSEMBL_TRANSCRIPT_ID'] = \
                    var2trans.groupby(['VAR_ID','consTPM_SINGLE','ENTREZGENE'])['ENSEMBL_TRANSCRIPT_ID'].transform(lambda x : ', '.join(x))
                var2trans = var2trans.loc[~var2trans['ENSEMBL_TRANSCRIPT_ID'].isna(), :]
                var2exp = varset.groupby(['VAR_ID','consTPM_SINGLE','ENTREZGENE']).agg({'TPM':'sum'}).reset_index()
                var2exp.columns = ['VAR_ID','consTPM_SINGLE','ENTREZGENE','consTPM']
                cons2exp = pd.merge(var2exp, var2trans, on = ['VAR_ID','consTPM_SINGLE','ENTREZGENE'], how = 'left').drop_duplicates()
                cons2exp = cons2exp.sort_values(by=['VAR_ID','consTPM'], ascending=[True,False])
                cons2exp.columns = ['VAR_ID','consTPM_SINGLE','ENTREZGENE','consTPM','consENSEMBL_TRANSCRIPT_ID']
                cons2exp.drop(['consENSEMBL_TRANSCRIPT_ID'], axis = 1, inplace = True)
                variant_set['consTPM_SINGLE'] = variant_set['CONSEQUENCE'].str.replace('&.+','', regex=True)
                
                if {'VAR_ID','ENTREZGENE','consTPM_SINGLE'}.issubset(variant_set.columns):                    
                    variant_set = variant_set.merge(cons2exp, on = ['VAR_ID','ENTREZGENE','consTPM_SINGLE'], how = 'left')
                    variant_set.drop(['consTPM_SINGLE'], axis = 1, inplace = True)
                  
        
    
    return(variant_set)

def correlate_sample_expression(sample_expression: dict, 
                                yaml_data: dict, 
                                refdata_assembly_dir: str,
                                protein_coding_only: bool = True,
                                logger: logging.Logger = None) -> pd.DataFrame:
    
    exp_sim = {}
    for k in yaml_data['conf']['expression']['similarity_db'].keys():
        exp_sim[k] =  pd.DataFrame()
    sample_id = yaml_data['sample_id']
    drop_columns = ['SYMBOL','BIOTYPE', 'GENENAME','ENTREZGENE','TPM_GENE']
    
    exp_data_sample = sample_expression.copy()
    
    if 'gene' in exp_data_sample.keys():
        if not exp_data_sample['gene'] is None:
            if {'ENSEMBL_GENE_ID','TPM_LOG2_GENE','BIOTYPE'}.issubset(exp_data_sample['gene'].columns):
                if protein_coding_only is True:
                    logger.info("Filtering out non-protein coding genes from expression data")
                    exp_data_sample['gene'] = \
                        exp_data_sample['gene'][exp_data_sample['gene']['BIOTYPE'] == 'protein_coding']
                if len(exp_data_sample['gene']) < 10:
                    logger.warning(
                        'Expression file contains limited protein-coding gene expression records (N = ' + \
                            str(len(exp_data_sample['gene'])) + ') - skipping correlation analysis')
                    return(exp_sim)
                else:
                    for col in drop_columns:
                        if col in exp_data_sample['gene'].columns:
                            exp_data_sample['gene'] = exp_data_sample['gene'].drop(col, axis = 1)   
                    exp_data_sample['gene'] = exp_data_sample['gene'].rename(
                        columns = {'TPM_LOG2_GENE':sample_id})
                
                    if 'tcga' in yaml_data['conf']['expression']['similarity_db'].keys():
                        for cohort in yaml_data['conf']['expression']['similarity_db']['tcga'].keys():                            
                            exp_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                             "tcga", "tcga_" + str(cohort).lower() + "_tpm.tsv.gz")
                            if check_file_exists(exp_fname, strict = False, logger = logger):
                                exp_data_refcohort = pd.read_csv(
                                    exp_fname, 
                                    sep = "\t", na_values = ".", 
                                    low_memory = False)
                                for ann in drop_columns:
                                    if ann in exp_data_refcohort.columns:
                                        exp_data_refcohort = exp_data_refcohort.drop(ann, axis = 1)
                                if 'ENSEMBL_GENE_ID' in exp_data_refcohort.columns:
                                    logger.info('Calculating gene expression correlation between ' + str(sample_id) + \
                                        ' and samples from The Cancer Genome Atlas (TCGA) cohorts - ' + str(cohort))
                                    corr_mat = correlate_samples(exp_data_sample, exp_data_refcohort, sample_id, 'tcga-' + str(cohort).lower(), 
                                                            protein_coding_only,  logger)
                                    ## Make sure sample IDs match patient barcodes (sample metadata is organized per patient)
                                    #corr_mat['EXT_SAMPLE_ID'] = corr_mat['EXT_SAMPLE_ID'].replace(regex=r'-[0-9][0-9][A-Z]$', value='')                                                                      
                                    metadata_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                        'tcga', "tcga_" + str(cohort).lower() + "_sample_metadata.tsv.gz")
                                    if check_file_exists(metadata_fname, strict = False, logger = logger):
                                        sample_metadata = pd.read_csv(
                                            metadata_fname, sep = "\t", na_values = ".", low_memory = False)
                                        corr_mat = corr_mat.merge(
                                            sample_metadata, on = ['EXT_SAMPLE_ID'], how = 'left') 
                                    exp_sim['tcga'] = pd.concat(
                                            [exp_sim['tcga'], corr_mat], ignore_index = True)
                                                                       
                        if not exp_sim['tcga'].empty:
                            if 'CORR' in exp_sim['tcga'].columns:
                                exp_sim['tcga'] = \
                                    exp_sim['tcga'].sort_values(by=['CORR'], ascending=False)
                        
                    for source in ['depmap','treehouse']:
                        source_verbose = 'DepMap cell lines'
                        if source == 'treehouse':
                            source_verbose = 'samples from the Treehouse Childhood Cancer Initiative'
                        if source in yaml_data['conf']['expression']['similarity_db'].keys():
                            exp_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                             source, source + "_tpm.tsv.gz")
                            if check_file_exists(exp_fname, strict = False, logger = logger):
                                exp_data_refcohort = pd.read_csv(
                                    exp_fname, 
                                    sep = "\t", na_values = ".", 
                                    low_memory = False)
                                for ann in drop_columns:
                                    if ann in exp_data_refcohort.columns:
                                        exp_data_refcohort = exp_data_refcohort.drop(ann, axis = 1)
                                logger.info('Calculating gene expression correlation between ' + str(sample_id) + ' and ' + source_verbose)
                                corr = correlate_samples(exp_data_sample, exp_data_refcohort, sample_id, source, 
                                                        protein_coding_only,  logger)
                                exp_sim[source] = pd.concat(
                                    [exp_sim[source], corr], ignore_index = True)
                                
                            if not exp_sim[source].empty:
                                if 'CORR' in exp_sim[source].columns:
                                    exp_sim[source] = \
                                        exp_sim[source].sort_values(by=['CORR'], ascending=False)
                                    metadata_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                        source, source + "_sample_metadata.tsv.gz")
                                    if check_file_exists(metadata_fname, strict = False, logger = logger):
                                        sample_metadata = pd.read_csv(
                                            metadata_fname, sep = "\t", na_values = ".", low_memory = False)
                                        exp_sim[source] = exp_sim[source].merge(
                                            sample_metadata, on = ['EXT_SAMPLE_ID'], how = 'left')
                                
                     
    return(exp_sim)


def find_expression_outliers(sample_expression: dict,
                        yaml_data: dict,
                        refdata_assembly_dir: str,
                        protein_coding_only: bool,
                        logger: logging.Logger) -> pd.DataFrame:
    
    sample_id = yaml_data['sample_id']
    primary_site = yaml_data['conf']['sample_properties']['site']
    required_cols = ['ENSEMBL_GENE_ID','TPM_LOG2_GENE']
    drop_columns = ['SYMBOL','BIOTYPE', 'GENENAME','ENTREZGENE']
    
    exp_data_sample = sample_expression.copy()
    exp_data_refcohort = pd.DataFrame()
    outlier_metrics_valid = pd.DataFrame()
    
    
    ## As of now (April 2024) - match primary site of input sample to TCGA cohorts
    ## Future - allow user to specify TCGA cohort to compare with
    ##
    if primary_site == "Any":
        logger.warning("Primary site not specified in configuration file - skipping expression outlier analysis")
        return(pd.DataFrame())
    else:
        if primary_site in pcgr_vars.SITE_TO_DISEASE.keys():
            comparison_disease_cohort = str(pcgr_vars.SITE_TO_DISEASE[primary_site][0]).lower()
            exp_fname = os.path.join(
                refdata_assembly_dir, "expression", "tsv", 
                "tcga", str(comparison_disease_cohort).lower() + "_tpm.tsv.gz")
            if check_file_exists(exp_fname, strict = False, logger = logger):
                exp_data_refcohort = pd.read_csv(exp_fname, sep = "\t", na_values = ".", low_memory = False)
                
                ## Filter for protein coding genes when doing outlier analysis
                if protein_coding_only:
                    #logger.info("Filtering out non-protein coding genes from expression data")
                    exp_data_refcohort = \
                        exp_data_refcohort[exp_data_refcohort['BIOTYPE'] == 'protein_coding']
                for col in exp_data_refcohort.columns:
                    if col in drop_columns:
                        exp_data_refcohort = exp_data_refcohort.drop(col, axis = 1)
            else:
                logger.warning(f'{exp_fname} not found - skipping expression outlier analysis')
                return(pd.DataFrame())
        else:
            logger.warning("Primary tumor site not specified in configuration file - skipping expression outlier analysis")
            return(pd.DataFrame())
    
    if 'gene' in exp_data_sample.keys():
        if not exp_data_sample['gene'] is None:
            if {'ENSEMBL_GENE_ID','TPM_LOG2_GENE'}.issubset(exp_data_sample['gene'].columns):
                for col in exp_data_sample['gene'].columns:
                    if col not in required_cols:
                        exp_data_sample['gene'] = exp_data_sample['gene'].drop(col, axis = 1)
                exp_data_sample['gene'] = exp_data_sample['gene'].rename(
                    columns = {'TPM_LOG2_GENE':sample_id})
                
                if 'ENSEMBL_GENE_ID' in exp_data_refcohort.columns and \
                    'ENSEMBL_GENE_ID' in exp_data_sample['gene'].columns:
        
                    ref_sample_mat = exp_data_refcohort.merge(
                        exp_data_sample['gene'], on = 'ENSEMBL_GENE_ID', how = 'left')
                    
                    ref_sample_mat = ref_sample_mat.set_index('ENSEMBL_GENE_ID')
                    percentiles = ref_sample_mat.rank(1, pct=True, numeric_only=True).apply(lambda x: round(x * 100, 1))                    
                    sample_percentiles = percentiles[[sample_id]].reset_index().rename(columns={sample_id: 'PERCENTILE'})                    
                    quantiles = pd.concat([
                        round(ref_sample_mat[sample_id], ndigits = 5),
                        round(ref_sample_mat.quantile(0.25, 1), ndigits = 5),
                        round(ref_sample_mat.quantile(0.5, 1), ndigits = 5),
                        round(ref_sample_mat.quantile(0.75, 1), ndigits = 5),
                        round(ref_sample_mat.mean(1), ndigits = 9),
                        round(ref_sample_mat.std(1), ndigits = 9),
                        ], axis=1).rename(
                            columns={0.25: 'Q1', 0.5: 'Q2',0.75: 'Q3', 
                                     sample_id: 'TPM_LOG2_GENE', 0: 'MEAN', 1: 'STD'}).reset_index()                    
                    quantiles['IQR'] = round(quantiles.Q3 - quantiles.Q1, ndigits = 5)
                    quantiles['Z_SCORE'] = round(quantiles[quantiles.STD > 0].apply(
                        lambda row: (row.TPM_LOG2_GENE - row.MEAN) / row.STD, axis=1), ndigits = 5)
                    #ref_sample_mat = exp_data_refcohort.merge(
                    #    exp_data_sample['gene'], on = 'ENSEMBL_GENE_ID', how = 'left')
                    #quantiles['kIQR'] = round(quantiles[quantiles.IQR > 0].apply(
                    #    lambda row: (row.TPM_GENE - row.Q2) / row.IQR, axis=1), ndigits = 5)
                    
                    outlier_metrics = \
                        quantiles.merge(
                            sample_percentiles, 
                            on = 'ENSEMBL_GENE_ID', how = 'left')
                    outlier_metrics['SAMPLE_ID'] = sample_id
                    outlier_metrics['REF_COHORT'] = comparison_disease_cohort
                    outlier_metrics['REF_COHORT_SIZE'] = exp_data_refcohort.shape[1] - 1
                    outlier_metrics = \
                        outlier_metrics[['SAMPLE_ID', 
                                         'REF_COHORT', 
                                         'REF_COHORT_SIZE', 
                                         'ENSEMBL_GENE_ID', 
                                         'TPM_LOG2_GENE', 
                                         'MEAN',
                                         'STD',
                                         'Z_SCORE',
                                         'Q1', 
                                         'Q2', 
                                         'Q3', 
                                         'IQR', 
                                         'PERCENTILE']]
                    
                    mask_valid_ensembl = pd.notna(outlier_metrics['TPM_LOG2_GENE'])
                    outlier_metrics_valid = outlier_metrics[mask_valid_ensembl]                                              
    
    
    return(outlier_metrics_valid)

def correlate_samples(exp_data_sample: dict, 
                      exp_data_refcohort: pd.DataFrame,
                      sample_id: str,
                      db: str,
                      protein_coding_only: bool,                       
                      logger: logging.Logger) -> pd.DataFrame:
    
    corr_similarity = pd.DataFrame()
    if not 'gene' in exp_data_sample:
        logger.warning("No 'gene' entry in expression data dictionary for sample " + sample_id + " - skipping correlation analysis")
        return(corr_similarity)
    if 'ENSEMBL_GENE_ID' in exp_data_refcohort.columns and \
        'ENSEMBL_GENE_ID' in exp_data_sample['gene'].columns:
        
        exp_data_refcohort = exp_data_refcohort.merge(
            exp_data_sample['gene'], on = 'ENSEMBL_GENE_ID', how = 'left')
        exp_data_refcohort = exp_data_refcohort.set_index('ENSEMBL_GENE_ID')
                        
        corr_raw = exp_data_refcohort[[c for c in exp_data_refcohort.columns.tolist() if c != sample_id]].corrwith(
            exp_data_refcohort[sample_id], method='spearman').to_frame()
        
        corr = corr_raw.rename(columns={ corr_raw.columns[0]: "CORR" }).reset_index().rename(
                columns={'index': 'EXT_SAMPLE_ID'})
        corr['SAMPLE_ID'] = sample_id
        corr['EXT_DB'] = db
        corr['PROTEIN_CODING_ONLY'] = protein_coding_only
        #corr['EXT_DB_N'] = exp_data_refcohort.shape[1] - 1
        
        corr_similarity = corr[['SAMPLE_ID','EXT_SAMPLE_ID', 'EXT_DB',
                     'CORR','PROTEIN_CODING_ONLY']]
    
    return(pd.DataFrame(corr_similarity))
            
