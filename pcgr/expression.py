#!/usr/bin/env python

import os,re
import csv
import gzip
import pandas as pd
import logging
import numpy as np

from pcgr import utils
from pcgr.utils import error_message, warn_message, check_file_exists, remove_file

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
    
    ## Read the gene identifier index - maps transcript identifiers (Ensembl/Refseq), 
    ## gene symbols/aliases and identifiers (Ensembl) to Entrez gene identifers and
    ## biotype
    gene_index_fname_tsv =  \
        os.path.join(refdata_assembly_dir, "gene", "tsv", "gene_transcript_xref", "gene_index.tsv.gz")
    check_file_exists(gene_index_fname_tsv, logger)
    gene_index = pd.read_csv(gene_index_fname_tsv, sep = "\t", na_values = ".",  low_memory = False)
    
    sample_identifiers_found = 0
    expression_map = {}
    expression_map['transcript'] = None
    expression_map['gene'] = None
    expression_map['transcript_identifier'] = None
    
    if {'ID','ID_TYPE','ENTREZGENE','ENSEMBL_GENE_ID','SYMBOL','BIOTYPE','AMBIGUOUS_ID'}.issubset(gene_index.columns):
        exp_map = pd.merge(sample_gene_expression, gene_index, how="left", on=["ID"])
        
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
            
            identifiers_used_in_input = "Ensembl gene IDs"
            transcript_level = False
            transcript_identifier = None
            if "refseq_transcript_id" in identifiers_present:
                identifiers_used_in_input = "RefSeq transcript IDs"
                transcript_level = True
                transcript_identifier = 'RefSeq'
            if "ensembl_transcript_id" in identifiers_present:
                identifiers_used_in_input = "Ensembl transcript IDs"
                transcript_level = True
                transcript_identifier = 'Ensembl'
            if "symbol_alias" in identifiers_present:
                identifiers_used_in_input = "gene symbols"
            
            ## Emit warning if more than 5% of gene/transcript identifiers are not properly verified
            sample_identifiers_found = len(exp_map_verified)
            percent_verified = round((len(exp_map_verified) / len(exp_map)) * 100, 2)
            percent_missing = 100 - percent_verified
            if percent_missing > 5:
                logger.warn("Failed to map " + str(percent_missing) +  \
                    "% of gene/transcript identifiers in input TSV file - use proper ENST/RefSeq identifiers")
            logger.info("Verified N = " + str(sample_identifiers_found) + " (" + str(percent_verified) + \
                "%) of gene/transcript identifiers in input gene expression file - using " + str(identifiers_used_in_input))
            
            ## Emit warning if ambiguous gene/transcript identifiers are detected - 
            ## remove them from the analysis (write them to a separate file?)
            n_ambig = len(exp_map_verified[exp_map_verified.AMBIGUOUS_ID == True])
            if n_ambig > 0:
                logger.warn("Detected N = " + str(n_ambig) + " ambiguous gene/transcript identifiers in input gene expression file")
            else:
                logger.info("NO ambiguous gene/transcript identifiers were detected in input gene expression file")
            transcript_expression_map = exp_map_verified[exp_map_verified.AMBIGUOUS_ID == False]
            #print(transcript_expression_map.head(10))
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
                                                  "ENSEMBL_GENE_ID","SYMBOL","ENTREZGENE","GENENAME",
                                                  "BIOTYPE"]]
                
                #print(expression_map['transcript'].head(10))
            
            ## make gene level TPM summary  
            expression_map['gene'] = transcript_expression_map.groupby(
                ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE']).agg({'TPM':'sum'}).reset_index()
            expression_map['gene'].columns = ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE','TPM_GENE']
            expression_map['gene'] = expression_map['gene'].drop_duplicates().sort_values(by='TPM_GENE', ascending=False)
            
            ## for transcripts with missing TPM's, use the minimum TPM value across all transcripts
            expression_map['gene_min_tpm'] = transcript_expression_map.groupby(
                ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE']).agg({'TPM':'min'}).reset_index()
            expression_map['gene_min_tpm'].columns = ['ENSEMBL_GENE_ID','SYMBOL','ENTREZGENE','GENENAME','BIOTYPE','TPM_MIN']
            expression_map['gene_min_tpm'] = expression_map['gene_min_tpm'].drop_duplicates().sort_values(by='TPM_MIN', ascending=False)
                  
    
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
                    logger.info("Integrating gene-level expression data from tumor into somatic variant set")
                    if expression_data[s].empty:
                        logger.warn('Expression file does not contain any gene-level expression data')
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
                    
        if 'transcript' in expression_data.keys() and 'transcript_identifier' in expression_data.keys():
            if not expression_data['transcript'] is None:
                logger.info("Integrating transcript-level expression data from tumor into somatic variant set")
                if expression_data['transcript'].empty:
                    logger.warn('Expression file does not contain any transcript-level expression data')
                else:
                    exp_data_minimal = expression_data['transcript'][['ID','TPM']]
                    transcript_identifier = 'REFSEQ_TRANSCRIPT_ID'
                    if expression_data['transcript_identifier'] == 'Ensembl':
                        transcript_identifier = 'ENSEMBL_TRANSCRIPT_ID'
                    if {'VAR_ID', transcript_identifier}.issubset(variant_set.columns):
                        if expression_data['transcript_identifier'] == 'Ensembl':
                            var2trans_all = variant_set[['VAR_ID','ENSEMBL_TRANSCRIPT_ID']]
                            #var2trans_all.loc[:, 'ID'] = var2trans_all['ENSEMBL_TRANSCRIPT_ID']
                            var2trans_all = var2trans_all.assign(ID = var2trans_all.ENSEMBL_TRANSCRIPT_ID)
                        else:
                            var2trans_all = variant_set[['VAR_ID','ENSEMBL_TRANSCRIPT_ID',transcript_identifier]].rename(
                                columns = {transcript_identifier:'ID'})
                        mask = pd.notna(var2trans_all['ID'])
                        var2trans = var2trans_all[mask]
                        if not var2trans.empty:
                            ## If multiple transcript identifiers are present (RefSeq), take the maximum value
                            var2trans.loc[:,'ID'] = var2trans.loc[:,'ID'].str.split('&')
                            var2trans = var2trans.explode('ID').drop_duplicates().reset_index(drop = True)                                    
                            var2trans = pd.merge(
                                var2trans, exp_data_minimal, on = 'ID', how = 'left')                                
                            var2exp = var2trans.groupby(['VAR_ID','ENSEMBL_TRANSCRIPT_ID']).agg({'TPM':'max'}).reset_index()
                            if 'TPM' in variant_set.columns:
                                variant_set.drop('TPM', inplace=True, axis=1)
                            if {'TPM','VAR_ID','ENSEMBL_TRANSCRIPT_ID'}.issubset(var2exp.columns):
                                var2exp = var2exp[['TPM','VAR_ID','ENSEMBL_TRANSCRIPT_ID']]
                                variant_set = variant_set.merge(var2exp, on = ['VAR_ID','ENSEMBL_TRANSCRIPT_ID'], how = 'left')  
                        else:
                            variant_set['TPM'] = np.nan
                            logger.warn('Variant file does not contain any entries with valid transcript identifiers')
    
    return(variant_set)

def correlate_sample_expression(sample_expression_data: dict, 
                                yaml_data: dict, 
                                refdata_assembly_dir: str,
                                protein_coding_only: bool = True,
                                logger: logging.Logger = None) -> pd.DataFrame:
    
    sample_exp_similarity = {}
    for k in yaml_data['conf']['gene_expression']['similarity_db'].keys():
        sample_exp_similarity[k] =  pd.DataFrame()
    sample_id = yaml_data['sample_id']
    drop_columns = ['SYMBOL','BIOTYPE', 'GENENAME','ENTREZGENE']
    
    if 'gene' in sample_expression_data.keys():
        if not sample_expression_data['gene'] is None:
            if {'ENSEMBL_GENE_ID','TPM_GENE','BIOTYPE'}.issubset(sample_expression_data['gene'].columns):
                if protein_coding_only is True:
                    logger.info("Filtering out non-protein coding genes from expression data")
                    #print('Pre-filter: ' + str(len(sample_expression_data['gene'])))
                    sample_expression_data['gene'] = \
                        sample_expression_data['gene'][sample_expression_data['gene']['BIOTYPE'] == 'protein_coding']
                    #print('Post-filter: ' + str(len(sample_expression_data['gene'])))
                if len(sample_expression_data['gene']) < 10:
                    logger.warn(
                        'Expression file contains limited protein-coding gene expression records (N = ' + \
                            str(len(sample_expression_data['gene'])) + ') - skipping correlation analysis')
                    return(sample_exp_similarity)
                else:
                    for col in drop_columns:
                        if col in sample_expression_data['gene'].columns:
                            sample_expression_data['gene'] = sample_expression_data['gene'].drop(col, axis = 1)   
                    sample_expression_data['gene'] = sample_expression_data['gene'].rename(
                        columns = {'TPM_GENE':sample_id})
                
                    if 'tcga' in yaml_data['conf']['gene_expression']['similarity_db'].keys():
                        for cohort in yaml_data['conf']['gene_expression']['similarity_db']['tcga'].keys():                            
                            exp_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                             "tcga", "tcga_" + str(cohort).lower() + "_tpm.tsv.gz")
                            if check_file_exists(exp_fname, strict = False, logger = logger):
                                sample_comp_exp_mat = pd.read_csv(
                                    exp_fname, 
                                    sep = "\t", na_values = ".", 
                                    low_memory = False)
                                for ann in drop_columns:
                                    if ann in sample_comp_exp_mat.columns:
                                        sample_comp_exp_mat = sample_comp_exp_mat.drop(ann, axis = 1)
                                if 'ENSEMBL_GENE_ID' in sample_comp_exp_mat.columns:
                                    logger.info('Calculating gene expression correlation between ' + str(sample_id) + \
                                        ' and samples from The Cancer Genome Atlas (TCGA) cohorts - ' + str(cohort))
                                    corr_mat = correlate_samples(sample_expression_data, sample_id, 'tcga-' + str(cohort).lower(), 
                                                            protein_coding_only, sample_comp_exp_mat, logger)
                                    ## Make sure sample IDs match patient barcodes (sample metadata is organized per patient)
                                    #corr_mat['EXT_SAMPLE_ID'] = corr_mat['EXT_SAMPLE_ID'].replace(regex=r'-[0-9][0-9][A-Z]$', value='')                                                                      
                                    metadata_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                        'tcga', "tcga_" + str(cohort).lower() + "_sample_metadata.tsv.gz")
                                    if check_file_exists(metadata_fname, strict = False, logger = logger):
                                        sample_metadata = pd.read_csv(
                                            metadata_fname, sep = "\t", na_values = ".", low_memory = False)
                                        corr_mat = corr_mat.merge(
                                            sample_metadata, on = ['EXT_SAMPLE_ID'], how = 'left') 
                                    sample_exp_similarity['tcga'] = pd.concat(
                                            [sample_exp_similarity['tcga'], corr_mat], ignore_index = True)
                                                                       
                        if not sample_exp_similarity['tcga'].empty:
                            if 'CORRELATION' in sample_exp_similarity['tcga'].columns:
                                sample_exp_similarity['tcga'] = \
                                    sample_exp_similarity['tcga'].sort_values(by=['CORRELATION'], ascending=False)
                        
                    for source in ['depmap','treehouse']:
                        source_verbose = 'DepMap cell lines'
                        if source == 'treehouse':
                            source_verbose = 'samples from the Treehouse Childhood Cancer Initiative'
                        if source in yaml_data['conf']['gene_expression']['similarity_db'].keys():
                            exp_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                             source, source + "_tpm.tsv.gz")
                            if check_file_exists(exp_fname, strict = False, logger = logger):
                                sample_comp_exp_mat = pd.read_csv(
                                    exp_fname, 
                                    sep = "\t", na_values = ".", 
                                    low_memory = False)
                                for ann in drop_columns:
                                    if ann in sample_comp_exp_mat.columns:
                                        sample_comp_exp_mat = sample_comp_exp_mat.drop(ann, axis = 1)
                                logger.info('Calculating gene expression correlation between ' + str(sample_id) + ' and ' + source_verbose)
                                corr = correlate_samples(sample_expression_data, sample_id, source, 
                                                        protein_coding_only, sample_comp_exp_mat, logger)
                                sample_exp_similarity[source] = pd.concat(
                                    [sample_exp_similarity[source], corr], ignore_index = True)
                                
                            if not sample_exp_similarity[source].empty:
                                if 'CORRELATION' in sample_exp_similarity[source].columns:
                                    sample_exp_similarity[source] = \
                                        sample_exp_similarity[source].sort_values(by=['CORRELATION'], ascending=False)
                                    metadata_fname = os.path.join(refdata_assembly_dir, "expression", "tsv", 
                                        source, source + "_sample_metadata.tsv.gz")
                                    if check_file_exists(metadata_fname, strict = False, logger = logger):
                                        sample_metadata = pd.read_csv(
                                            metadata_fname, sep = "\t", na_values = ".", low_memory = False)
                                        sample_exp_similarity[source] = sample_exp_similarity[source].merge(
                                            sample_metadata, on = ['EXT_SAMPLE_ID'], how = 'left')
                                
    #else:
        #print("BALLE")                    
    return(sample_exp_similarity)


def correlate_samples(sample_expression_data: dict, 
                      sample_id: str,
                      db: str,
                      protein_coding_only: bool, 
                      sample_comp_exp_mat: pd.DataFrame,
                      logger: logging.Logger) -> pd.DataFrame:
    
    corr_similarity = pd.DataFrame()
    if not 'gene' in sample_expression_data:
        logger.warn("No 'gene' entry in expression data dictionary for sample " + sample_id + " - skipping correlation analysis")
        return(corr_similarity)
    if 'ENSEMBL_GENE_ID' in sample_comp_exp_mat.columns and \
        'ENSEMBL_GENE_ID' in sample_expression_data['gene'].columns:
        
        sample_comp_exp_mat = sample_comp_exp_mat.merge(
            sample_expression_data['gene'], on = 'ENSEMBL_GENE_ID', how = 'left')
        sample_comp_exp_mat = sample_comp_exp_mat.set_index('ENSEMBL_GENE_ID')
                        
        corr_raw = sample_comp_exp_mat[[c for c in sample_comp_exp_mat.columns.tolist() if c != sample_id]].corrwith(
            sample_comp_exp_mat[sample_id], method='spearman').to_frame()
        
        corr = corr_raw.rename(columns={ corr_raw.columns[0]: "CORRELATION" }).reset_index().rename(
                columns={'index': 'EXT_SAMPLE_ID'})
        corr['SAMPLE_ID'] = sample_id
        corr['EXT_DB'] = db
        corr['PROTEIN_CODING_ONLY'] = protein_coding_only
        corr['EXT_DB_N'] = sample_comp_exp_mat.shape[1] - 1
        
        corr_similarity = corr[['SAMPLE_ID','EXT_SAMPLE_ID', 'EXT_DB',
                     'CORRELATION','PROTEIN_CODING_ONLY','EXT_DB_N']]
    
    return(pd.DataFrame(corr_similarity))
            
