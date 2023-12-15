#!/usr/bin/env python

import os,re
import csv
import gzip
import pybedtools
import pandas as pd
import logging

from pcgr import utils
from pybedtools import BedTool
from pcgr.annoutils import nuclear_chromosomes
from pcgr.utils import error_message, warn_message, check_file_exists
from pcgr.biomarker import load_biomarkers

def annotate_cna_segments(output_fname: str, 
                          output_dir: str, 
                          cna_segment_file: str, 
                          pcgr_build_db_dir: str, 
                          build: str, 
                          sample_id: str,
                          n_copy_amplifications: int = 5,
                          overlap_fraction: float = 0.5,
                          logger: logging.Logger = None) -> int:
    """
    Annotate copy number aberrations in a given segment file.
    Args:
        output_fname (str): File name of the annotated output file.
        output_dir (str): Directory to save the annotated file.
        cna_segment_file (str): Path to the user-provided copy number aberrations segment file.
        pcgr_build_db_dir (str): Path to the build-specific PCGR database directory.
        build (str): Genome assembly build of input data.
        sample_id( (str): Sample identifier
        output_dir (str): Directory to save the annotated file.
        n_copy_amplifications (int, optional): Number of copies to consider as gains/amplifications. Defaults to 5.
        overlap_fraction (float, optional): Fraction of overlap required for annotation. Defaults to 0.5.
    Returns:
        int: 0 if successful.
    """
    
    #logger = getlogger("pcgr-annotate-cna-segments")
    pybedtools.set_tempdir(output_dir)
    
    temp_files = []
    
    # Read user-defined copy number aberrations segment file
    check_file_exists(cna_segment_file, logger)
    cna_query_segment_df = pd.read_csv(cna_segment_file, sep="\t", na_values=".")
    
    ## Check that required columns are present
    if not {'Chromosome','Start','End','nMajor','nMinor'}.issubset(cna_query_segment_df.columns):
        err_msg = f"Could not find required columns in CNA segment file: {cna_segment_file} - exiting."
        error_message(err_msg, logger)
    
    ## Remove 'chr' prefix from chromosome names
    for elem in ['Chromosome','nMajor','nMinor']:
        cna_query_segment_df = cna_query_segment_df.astype({elem:'string'})
    cna_query_segment_df['Chromosome'] = \
        cna_query_segment_df['Chromosome'].str.replace("^(C|c)hr", "", regex = True)
    
    ## Only consider chromosome names that are present in the nuclear chromosomes list
    cna_query_segment_df = cna_query_segment_df[cna_query_segment_df['Chromosome'].isin(nuclear_chromosomes)]
    if cna_query_segment_df.empty is True:
        warn_msg = f"Could not find any CNA query segments listed on nuclear chromosomes: {nuclear_chromosomes} - returning."
        warn_message(warn_msg, logger)
        return 0
    
    ## Create segment identifier column
    cna_query_segment_df['segment_id'] = \
        'chr' + cna_query_segment_df['Chromosome'].str.cat(
            cna_query_segment_df['Start'].astype(str), sep = ":").str.cat(
                cna_query_segment_df['End'].astype(str), sep='-')
    
    ## Create Name column of BED file
    cna_query_segment_df["Name"] = \
        cna_query_segment_df["segment_id"].str.cat(
            cna_query_segment_df["nMajor"], sep="|").str.cat(
                cna_query_segment_df["nMinor"], sep = "|")
    cna_query_segment_df = cna_query_segment_df[['Chromosome','Start','End','Name']]
    
    
    ## Check that 'End' of segments do not exceed chromosome lengths
    chromsizes_fname = \
        os.path.join(pcgr_build_db_dir, 'chromsize.' + build + '.tsv')
    
    check_file_exists(chromsizes_fname, logger)
    chromsizes = pd.read_csv(chromsizes_fname, sep="\t", header=None, names=['Chromosome', 'ChromLength'])
    cna_query_segment_df = cna_query_segment_df.merge(
        chromsizes, left_on=["Chromosome"], right_on=["Chromosome"], how="left")
    segments_beyond_chromlength = \
        cna_query_segment_df[cna_query_segment_df['End'] > cna_query_segment_df['ChromLength']]
    
    ## Issue warning if segments exceed chromosome lengths
    if not segments_beyond_chromlength.empty is True:
        warn_msg = f"Ignoring n = {len(segments_beyond_chromlength)} copy number segments that " + \
            f"exceed the chromosomal lengths of {build}"
        warn_message(warn_msg, logger)
    cna_query_segment_df = \
        cna_query_segment_df[cna_query_segment_df['End'] <= cna_query_segment_df['ChromLength']]
    
    cna_query_segment_df = cna_query_segment_df[['Chromosome','Start','End','Name']]
    
    ## transform cna segments to pybedtools object
    cna_query_segment_bed = pybedtools.BedTool.from_dataframe(cna_query_segment_df)
    temp_files.append(cna_query_segment_bed.fn)
    
    ## annotate segments with cytobands
    cna_query_segment_df = annotate_cytoband(cna_query_segment_bed, output_dir, pcgr_build_db_dir, logger)

    ## annotate with protein-coding transcripts
    cna_query_segment_bed = pybedtools.BedTool.from_dataframe(cna_query_segment_df)
    temp_files.append(cna_query_segment_bed.fn)
    
    cna_query_segment_df = annotate_transcripts(
       cna_query_segment_bed, output_dir, pcgr_build_db_dir, overlap_fraction=overlap_fraction, logger=logger)

    cna_query_segment_df['segment_length_mb'] = \
        ((cna_query_segment_df['segment_end'] - cna_query_segment_df['segment_start']) / 1e6).astype(float).round(5)
    
    ## load copy-number biomarker evidence
   
    biomarkers = {}
    cna_actionable_dict = {}
    
    for db in ['cgi','civic']:
        variant_fname = os.path.join(pcgr_build_db_dir, 'biomarker','tsv', f"{db}.variant.tsv.gz")
        clinical_fname = os.path.join(pcgr_build_db_dir, 'biomarker','tsv', f"{db}.clinical.tsv.gz")
        logger.info(f"Loading copy-number biomarker evidence from {db} ..")
        biomarkers[db] = load_biomarkers(
            logger, variant_fname, clinical_fname, biomarker_vartype = 'CNA')
        
        for key in biomarkers[db]['other']:
            biomarker_data = biomarkers[db]['other'][key]
            biomarker_item = str(db) + '|' + str(biomarker_data[0]['variant_id']) + \
                    '|' + str(biomarker_data[0]['clinical_evidence_items']) + '|by_cna_segment'
            if not key in cna_actionable_dict:               
                cna_actionable_dict[key] = biomarker_item
            else:
                cna_actionable_dict[key] = cna_actionable_dict[key] + ',' + biomarker_item
            
    cna_actionable_df = pd.DataFrame(cna_actionable_dict.items(), columns=['aberration_key', 'biomarker_match'])
    
    ## Mark copy number amplifications (threshold defined by user) in input
    cna_query_segment_df['aberration_key'] = 'nan'
    cna_query_segment_df['loss_cond'] = True
    cna_query_segment_df.loc[cna_query_segment_df['n_major'] + cna_query_segment_df['n_minor'] < n_copy_amplifications,"amp_cond"] = False
    cna_query_segment_df.loc[cna_query_segment_df['n_major'] + cna_query_segment_df['n_minor'] >= n_copy_amplifications,"amp_cond"] = True
    
    cna_query_segment_df.loc[cna_query_segment_df.amp_cond, 'aberration_key'] =  \
        cna_query_segment_df.loc[cna_query_segment_df.amp_cond, 'entrezgene'].astype(str) + '_amplification'
    
    ## Mark homozygous deletions in input
    cna_query_segment_df['amp_cond'] = True
    cna_query_segment_df.loc[cna_query_segment_df['n_major'] + cna_query_segment_df['n_minor'] > 0,"loss_cond"] = False
    cna_query_segment_df.loc[cna_query_segment_df['n_major'] + cna_query_segment_df['n_minor'] == 0,"loss_cond"] = True
    
    cna_query_segment_df.loc[cna_query_segment_df.loss_cond, 'aberration_key'] =  \
        cna_query_segment_df.loc[cna_query_segment_df.loss_cond, 'entrezgene'].astype(str) + '_ablation'

    ## Append actionability evidence to input amplifications (column 'biomarker_match')
    cna_query_segment_df = cna_query_segment_df.merge(
        cna_actionable_df, left_on=["aberration_key"], right_on=["aberration_key"], how="left")
    cna_query_segment_df.drop(['amp_cond', 'loss_cond', 'aberration_key'], axis=1, inplace=True)    
    cna_query_segment_df.loc[cna_query_segment_df['biomarker_match'].isnull(),"biomarker_match"] = '.'
    
    ## remove all temporary files
    for fname in temp_files:
        utils.remove(fname)
        
    cna_query_segment_df.columns = map(str.upper, cna_query_segment_df.columns)
    cna_query_segment_df.rename(columns = {'CHROMOSOME':'CHROM','SEGMENT_ID':'VAR_ID'}, inplace = True)
    hgname = "hg38"
    if build == "grch37":
        hgname = "hg19"
    ucsc_browser_prefix = \
        f"http://genome.ucsc.edu/cgi-bin/hgTracks?db={hgname}&position="
        
    cna_query_segment_df['SEGMENT_LINK'] = \
        "<a href='" + ucsc_browser_prefix + cna_query_segment_df['VAR_ID'].astype(str) + \
            "' target='_blank'>" + cna_query_segment_df['VAR_ID'].astype(str) + "</a>"
    cna_query_segment_df['VAR_ID'] = \
        cna_query_segment_df['VAR_ID'].str.cat(
            cna_query_segment_df['N_MAJOR'].astype(str), sep=":").str.cat(
                cna_query_segment_df['N_MINOR'].astype(str), sep=":")
    
    cna_query_segment_df['SAMPLE_ID'] = sample_id
    cna_query_segment_df.to_csv(output_fname, sep="\t", header=True, index=False)
                
    return 0


def annotate_cytoband(cna_segments_bt: BedTool, output_dir: str, pcgr_build_db_dir: str, logger: logging.Logger) -> pd.DataFrame:
    
    pybedtools.set_tempdir(output_dir)    
    temp_files = []
    
    # BED file with cytoband annotations
    cytoband_bed_fname = \
        os.path.join(pcgr_build_db_dir, 'misc','bed','cytoband', 'cytoband.bed.gz')
    
    cytoband_annotated_segments = pd.DataFrame()
    
    check_file_exists(cytoband_bed_fname, logger)
    cytoband_bed = pybedtools.BedTool(cytoband_bed_fname)
        
    tmp = cna_segments_bt.intersect(cytoband_bed, loj=True)        
    cytoband_intersect = pd.read_csv(
        tmp.fn, sep="\t", header=None, 
        names=['chromosome','segment_start','segment_end','segment_name',
                'cytoband_chrom','cytoband_start','cytoband_end','cytoband_name'])
    temp_files.append(tmp.fn)
    tmp1 = cytoband_intersect.groupby(['chromosome','segment_start','segment_end']).first()
    segments_first = tmp1.index.to_frame()
    tmp2 = cytoband_intersect.groupby(['chromosome','segment_start','segment_end']).last()
    segments_last = tmp2.index.to_frame()
        
    cytoband_first = pd.concat([segments_first, tmp1], axis = 1)
    cytoband_first = cytoband_first.reset_index(drop = True)
    cytoband_first = \
        cytoband_first[['chromosome','segment_start','segment_end','segment_name','cytoband_name']]
    cytoband_first.rename({'cytoband_name': 'cytoband_first'}, axis=1, inplace=True)
        
    cytoband_last = pd.concat([segments_last, tmp2], axis = 1)
    cytoband_last = cytoband_last.reset_index(drop = True)
    cytoband_last = \
        cytoband_last[['chromosome','segment_start','segment_end','segment_name','cytoband_name']]
    cytoband_last.rename({'cytoband_name': 'cytoband_last'}, axis=1, inplace=True)
        
    segments_cytoband = cytoband_first.merge(
        cytoband_last, left_on=["chromosome","segment_start","segment_end","segment_name"], 
        right_on=["chromosome","segment_start","segment_end","segment_name"], how="left")
        
    cytoband_first_annotations = segments_cytoband.cytoband_first.str.split('\\|', expand = True)
    cytoband_first_annotations.columns = ['first_cytoband','first_arm','first_arm_length','first_focal_threshold']
    cytoband_first_annotations = cytoband_first_annotations.astype({'first_focal_threshold':'int'})
    cytoband_last_annotations = segments_cytoband.cytoband_last.str.split('\\|', expand = True)
    cytoband_last_annotations.columns = ['last_cytoband','last_arm','last_arm_length','last_focal_threshold']
        
    cytoband_all = pd.concat([segments_cytoband, cytoband_first_annotations, cytoband_last_annotations], axis = 1)
    cytoband_all['segment_start'] = cytoband_all['segment_start'].astype(int)
    cytoband_all['segment_end'] = cytoband_all['segment_end'].astype(int)       
    cytoband_all.loc[:,'segment_length'] = cytoband_all['segment_end'] - cytoband_all['segment_start']
    cytoband_all['event_type'] =  'broad'
    cytoband_all.loc[cytoband_all['segment_length'] < cytoband_all['first_focal_threshold'], "event_type"] = "focal"
        
    cytoband_all['cytoband'] = cytoband_all['first_cytoband'].where(
    cytoband_all['first_cytoband'] == cytoband_all['last_cytoband'], 
    cytoband_all['first_cytoband'] + "-" + cytoband_all['last_cytoband'])
        
    cytoband_all['segment_name'] = cytoband_all['segment_name'].str.cat(cytoband_all['first_arm'], sep = "|").str.cat(
        cytoband_all['cytoband'],sep="|").str.cat(cytoband_all['event_type'], sep="|")    
    cytoband_annotated_segments = cytoband_all[['chromosome','segment_start','segment_end','segment_name']]
    
    ## remove all temporary files
    for fname in temp_files:
        utils.remove(fname)
            
    return cytoband_annotated_segments


def annotate_transcripts(cna_segments_bt: BedTool, output_dir: str,
                         pcgr_build_db_dir: str, overlap_fraction: float,
                         logger: logging.Logger) -> pd.DataFrame:
    
    """
    Annotate the CNA segments with gene transcripts and return the annotated segments as a DataFrame.
    
    Parameters:
    - cna_segments_bt: A BedTool object representing the CNA segments.
    - output_dir: A string specifying the output directory (for temporary files generation).
    - pcgr_build_db_dir: A string specifying the directory of the build-specific PCGR data bundle.
    - overlap_fraction: A float representing the fraction of overlap required for annotation.
    
    Returns:
    - cna_segments_annotated: A DataFrame containing the annotated CNA segments.
    """
    
    pybedtools.set_tempdir(output_dir)
    temp_files = []
        
    # BED file with protein-coding transcripts
    gene_transcript_bed_fname = \
        os.path.join(pcgr_build_db_dir, 'gene','bed','gene_transcript_xref', 'gene_transcript_xref_pc_nopad.bed.gz')
    gene_xref_tsv_fname = \
        os.path.join(pcgr_build_db_dir, "gene", "tsv", "gene_transcript_xref", "gene_transcript_xref.tsv.gz")
        
    cna_segments_annotated = pd.DataFrame()

    if os.path.exists(gene_transcript_bed_fname):
        gene_transcript_bed = pybedtools.BedTool(gene_transcript_bed_fname)               
        
        ## Intersect query segments with annotated gene transcripts
        transcript_cna_annotations = gene_transcript_bed.intersect(cna_segments_bt, wao=True, f=overlap_fraction)
        temp_files.append(transcript_cna_annotations.fn)
        
        if os.path.exists(transcript_cna_annotations.fn) and os.path.getsize(transcript_cna_annotations.fn) > 0:
            colnames = ['chromosome','transcript_start','transcript_end','transcript_annotations',
                        'chromosome2','segment_start','segment_end','segment_name','bp_overlap']
            cna_transcript_annotations = \
                pd.read_csv(transcript_cna_annotations.fn, sep="\t", na_values=".", 
                            names=colnames, header=None, low_memory=False)
            
            ## ignore transcripts that do not overlap with any copy number segment
            cna_transcript_annotations = \
                cna_transcript_annotations[cna_transcript_annotations['segment_start'] != -1]
            
            ## calculate fraction of overlap between segments and transcripts
            cna_transcript_annotations['transcript_overlap_percent'] = \
                cna_transcript_annotations['transcript_end'] - cna_transcript_annotations['transcript_start']
            cna_transcript_annotations['transcript_overlap_percent'] = ((cna_transcript_annotations['bp_overlap'] / \
                cna_transcript_annotations['transcript_overlap_percent']) * 100).astype(float).round(2)
                               
            ## pull out segment identifier, major and minor copy number from the name column, and cytoband information
            major_minor_annotations_df = cna_transcript_annotations.segment_name.str.split('\\|', expand = True)
            major_minor_annotations_df.columns = ['segment_id','n_major','n_minor','chromosome_arm','cytoband','event_type']
            cna_transcript_annotations = pd.concat([cna_transcript_annotations, major_minor_annotations_df], axis = 1)
            
            ## select only relevant columns
            core_segment_annotations = \
                cna_transcript_annotations[['chromosome','segment_start','segment_end','segment_id',
                                            'n_major','n_minor','chromosome_arm','cytoband','event_type',
                                            'transcript_overlap_percent',
                                            'transcript_start','transcript_end']]
            core_segment_annotations = core_segment_annotations.astype({'n_major':'int','n_minor':'int'})
            ## pull out gene cross-reference annotations
            gene_annotations = cna_transcript_annotations.transcript_annotations.str.split('\\|', expand = True)
            gene_annotations.columns = ["ensembl_transcript_id", "ensembl_gene_id","symbol","entrezgene",
                                      "refseq_transcript_id","actionable_gene", "tsg", "tsg_support","tsg_rank",
                                      "oncogene","oncogene_support","oncogene_rank"]
            for elem in ['actionable_gene','oncogene','tsg']:
                gene_annotations = gene_annotations.astype({elem:'string'})
                gene_annotations.loc[gene_annotations[elem] == "", elem] = "FALSE"
                gene_annotations.loc[gene_annotations[elem] == "1", elem] = "TRUE"
            
            for elem in ['symbol','ensembl_gene_id','entrezgene','refseq_transcript_id','tsg_support',
                         'tsg_rank','oncogene_support','oncogene_rank']:
                gene_annotations = gene_annotations.astype({elem:'string'})
                gene_annotations.loc[gene_annotations[elem] == "", elem] = "."
            
            cna_segments_annotated = pd.concat([core_segment_annotations, gene_annotations], axis = 1)
            
            if {'genename'}.issubset(cna_segments_annotated.columns):
                cna_segments_annotated.drop('genename', inplace=True, axis=1)
            
            if os.path.exists(gene_xref_tsv_fname):
                gene_xref_df = pd.read_csv(gene_xref_tsv_fname, sep="\t", na_values=".", usecols=["entrezgene","name"])
                gene_xref_df = gene_xref_df[gene_xref_df['entrezgene'].notnull()].drop_duplicates()
                gene_xref_df["entrezgene"] = gene_xref_df["entrezgene"].astype("int64").astype("string")
                gene_xref_df.rename(columns = {'name':'genename'}, inplace = True)                                        
                cna_segments_annotated = cna_segments_annotated.merge(
                    gene_xref_df, left_on=["entrezgene"], right_on=["entrezgene"], how="left")
                cna_segments_annotated["entrezgene"] = cna_segments_annotated['entrezgene'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                cna_segments_annotated = cna_segments_annotated.fillna('.')
            else:
                err_msg = f"Could not find {gene_xref_tsv_fname} needed for gene name annotation - exiting"
                error_message(err_msg, logger)                
        else:
            err_msg = f"{transcript_cna_annotations.fn} is empty or not found - exiting"
            error_message(err_msg, logger)
    else:
        err_msg = f"{gene_transcript_bed_fname} is empty or not found - exiting"
        error_message(err_msg, logger)
    
    ## remove all temporary files
    for fname in temp_files:
        utils.remove(fname)
    
    return(cna_segments_annotated)
        


def is_valid_cna(cna_segment_file, logger):
    """
    Function that checks whether the CNA segment file adheres to the correct format
    """
    cna_reader = csv.DictReader(open(cna_segment_file,'r'), delimiter='\t')
    ## check that required columns are present
    if not ('Chromosome' in cna_reader.fieldnames and 
            'Segment_Mean' in cna_reader.fieldnames and 
            'Start' in cna_reader.fieldnames and 
            'End' in cna_reader.fieldnames):
        err_msg = "Copy number segment file (" + str(cna_segment_file) + \
            ") is missing required column(s): 'Chromosome', 'Start', 'End', or  'Segment_Mean'\n. " + \
                "Column names present in file: " + str(cna_reader.fieldnames)
        return error_message(err_msg, logger)

    cna_dataframe = pd.read_csv(cna_segment_file, sep="\t")
    if cna_dataframe.empty is True:
        err_msg = 'Copy number segment file is empty - contains NO segments'
        return error_message(err_msg, logger)
    if not cna_dataframe['Start'].dtype.kind in 'i': ## check that 'Start' is of type integer
        err_msg = '\'Start\' column of copy number segment file contains non-integer values'
        return error_message(err_msg, logger)
    if not cna_dataframe['End'].dtype.kind in 'i': ## check that 'End' is of type integer
        err_msg = '\'End\' column of copy number segment file contains non-integer values'
        return error_message(err_msg, logger)

    if not cna_dataframe['Segment_Mean'].dtype.kind in 'if': ## check that 'Segment_Mean' is of type integer/float
        err_msg = '\'Segment_Mean\' column of copy number segment file contains non-numerical values'
        return error_message(err_msg, logger)

    for rec in cna_reader:
        if int(rec['End']) < int(rec['Start']): ## check that 'End' is always greather than 'Start'
            err_msg = 'Detected wrongly formatted chromosomal segment - \'Start\' is greater than \'End\' (' + \
                str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
            return error_message(err_msg, logger)
        if int(rec['End']) < 1 or int(rec['Start']) < 1: ## check that 'Start' and 'End' is always non-negative
            err_msg = 'Detected wrongly formatted chromosomal segment - \'Start\' or \'End\' is less than or equal to zero (' + \
                str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
            return error_message(err_msg, logger)
    logger.info(f'Copy number segment file ({cna_segment_file}) adheres to the correct format')
    return 0

def is_valid_cna2(cna_segment_file, logger):
    """
    Function that checks whether the CNA segment file adheres to the correct format
    """
    cna_reader = csv.DictReader(open(cna_segment_file,'r'), delimiter='\t')
    ## check that required columns are present
    if not ('Chromosome' in cna_reader.fieldnames and 
            'nMajor' in cna_reader.fieldnames and 
            'nMinor' in cna_reader.fieldnames and 
            'Start' in cna_reader.fieldnames and 
            'End' in cna_reader.fieldnames):
        err_msg = "Copy number segment file (" + str(cna_segment_file) + ") is missing one or more required column(s): " + \
            "'Chromosome', 'Start', 'End', 'nMinor', or 'nMajor'\n. Column names present in file: " + str(cna_reader.fieldnames)
        return error_message(err_msg, logger)

    cna_dataframe = pd.read_csv(cna_segment_file, sep="\t")
    if cna_dataframe.empty is True:
        err_msg = 'Copy number segment file is empty - contains NO segments'
        return error_message(err_msg, logger)
    
    for elem in ['Start','End','nMajor','nMinor']:
        if not cna_dataframe[elem].dtype.kind in 'i':
            err_msg = 'Copy number segment file contains non-integer values for column: "' + elem + '"'
            return error_message(err_msg, logger)

    for rec in cna_reader:
        if int(rec['End']) < int(rec['Start']): ## check that 'End' is always greather than 'Start'
            err_msg = 'Detected wrongly formatted chromosomal segment - \'Start\' is greater than \'End\' (' + \
                str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
            return error_message(err_msg, logger)
        if int(rec['End']) < 1 or int(rec['Start']) < 1: ## check that 'Start' and 'End' is always non-negative
            err_msg = 'Detected wrongly formatted chromosomal segment - \'Start\' or \'End\' is less than or equal to zero (' + \
                str(rec['Chromosome']) + ':' + str(rec['Start']) + '-' + str(rec['End']) + ')'
            return error_message(err_msg, logger)
    logger.info(f'Copy number segment file ({cna_segment_file}) adheres to the correct format')
    return 0
