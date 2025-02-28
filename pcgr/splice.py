#!/usr/bin/env python

import os,re,sys
import csv
from logging import Logger


def load_splice_effects(mutsplice_db_fname: str, logger: Logger):
    """
    Load splice variants with functional effects (MutSpliceDB) from a file.
    """
   
    splice_variants = {}
    if not os.path.exists(mutsplice_db_fname):
        logger.info(f"ERROR: File '{mutsplice_db_fname}' does not exist - exiting")
        exit(1)

    with open(mutsplice_db_fname, mode='rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            key = str(row['ENTREZGENE']) + '_' + str(row['REFSEQ_TRANSCRIPT_ID']) + '_' + str(row['HGVSC'])
            splice_variants[key] = \
                str(row['ENTREZGENE']) + '|' + \
                str(row['SYMBOL']) + '|' + \
                str(row['REFSEQ_TRANSCRIPT_ID']) + '|' + \
                str(row['HGVSC']) + '|' + \
                str(row['SPLICE_EFFECT']) + '|' + \
                str(row['AFFECTED_EXONS']) + '|' + \
                str(row['SOURCE'])
         
    return splice_variants