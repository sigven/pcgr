#!/usr/bin/env python

import os
import pandas as pd

from pcgr.utils import check_file_exists, remove_file

def update_maf(maf_tmp_fname: str, 
               maf_fname: str, 
               allelic_support_tags: dict,
               logger = None,
               update_allelic_support = False,
               debug = False):

    """
    Update MAF file from vcf2maf.pl with allelic support data (t_depth, t_ref_count, t_alt_count etc).
    
    Args:
        maf_tmp_fname (str): File name of the temporary MAF file.
        maf_fname (str): File name of the final MAF file.
        allelic_support_tags (dict): Dictionary of allelic support tags (encoded in VCF INFO field and retained in MAF).
        logger: Logger object for logging messages.
        update_allelic_support (bool): Flag indicating whether to update allelic support.
    Returns:
        int: 0 if successful.
    """
    
    # Read MAF file generated with vcf2maf.pl
    check_file_exists(maf_tmp_fname, logger)
    
    header_line = "#version 2.4"
    with open(maf_tmp_fname) as f:
        header_line = f.readline().strip('\n')
    f.close()
        
    raw_maf_data = pd.read_csv(maf_tmp_fname, sep="\t", header=1, dtype='string',na_values=['.'], low_memory=False)
    if update_allelic_support is False:
        # write to file
        os.rename(maf_tmp_fname, maf_fname)
    else:
    
        if 'tumor_dp_tag' in allelic_support_tags:
            if allelic_support_tags['tumor_dp_tag'] != "_NA_":                               
                if {allelic_support_tags['tumor_dp_tag']}.issubset(raw_maf_data.columns):
                    if raw_maf_data[raw_maf_data[allelic_support_tags['tumor_dp_tag']].isna()].empty:
                        raw_maf_data.loc[:,"t_depth"] = raw_maf_data.loc[:,allelic_support_tags['tumor_dp_tag']]
                        
                        if 'tumor_af_tag' in allelic_support_tags: 
                            if allelic_support_tags['tumor_af_tag'] != "_NA_":    
                                if {allelic_support_tags['tumor_af_tag']}.issubset(raw_maf_data.columns):                                   
                                    if raw_maf_data[raw_maf_data[allelic_support_tags['tumor_af_tag']].isna()].empty:
                                        raw_maf_data['t_alt_count'] = None
                                        raw_maf_data.loc[:,"t_alt_count"] = \
                                            raw_maf_data.loc[:,allelic_support_tags['tumor_af_tag']].astype(float) * raw_maf_data.loc[:,"t_depth"].astype(int)
                                    
                                        raw_maf_data.loc[:,"t_alt_count"] = round(raw_maf_data.loc[:,"t_alt_count"].astype(float),0).astype(int)
                                        raw_maf_data['t_ref_count'] = None
                                        raw_maf_data.loc[:,"t_ref_count"] = \
                                            raw_maf_data.loc[:,"t_depth"].astype(int) - raw_maf_data.loc[:,"t_alt_count"]
        
        if 'control_dp_tag' in allelic_support_tags:
            if allelic_support_tags['control_dp_tag'] != "_NA_":                    
                if {allelic_support_tags['control_dp_tag']}.issubset(raw_maf_data.columns):
                    if raw_maf_data[raw_maf_data[allelic_support_tags['control_dp_tag']].isna()].empty:
                        raw_maf_data.loc[:,"n_depth"] = raw_maf_data.loc[:,allelic_support_tags['control_dp_tag']]
                        
                        if 'control_af_tag' in allelic_support_tags: 
                            if allelic_support_tags['control_af_tag'] != "_NA_":    
                                if {allelic_support_tags['control_af_tag']}.issubset(raw_maf_data.columns):
                                    if raw_maf_data[raw_maf_data[allelic_support_tags['control_af_tag']].isna()].empty:
                                        raw_maf_data['n_alt_count'] = None
                                        raw_maf_data.loc[:,"n_alt_count"] = \
                                            raw_maf_data.loc[:,allelic_support_tags['control_af_tag']].astype(float) * raw_maf_data.loc[:,"n_depth"].astype(int)
                                    
                                        raw_maf_data.loc[:,"n_alt_count"] = round(raw_maf_data.loc[:,"n_alt_count"].astype(float),0).astype(int)
                                        raw_maf_data['n_ref_count'] = None
                                        raw_maf_data.loc[:,"n_ref_count"] = \
                                            raw_maf_data.loc[:,"n_depth"].astype(int) - raw_maf_data.loc[:,"n_alt_count"]
        
        #raw_maf_data = raw_maf_data.fillna("")
        with open(maf_fname, 'w') as f:
            f.write(f'{header_line}\n')
        f.close()
        raw_maf_data.to_csv(maf_fname, sep="\t", index=False, mode='a')
        if not debug:
            remove_file(maf_tmp_fname)                     
                
                
        
            
    return 0