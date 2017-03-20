#!/usr/bin/env python

import os,re,sys
import logging
import csv
import gzip
from bx.intervals.intersection import IntervalTree

csv.field_size_limit(500 * 1024 * 1024)


def index_gene(gene_annotations_file_path, index = 'ensGene_transcript'):
   
   gene_xref = {}
   with gzip.open(gene_annotations_file_path, 'rb') as tsvfile:
      greader = csv.DictReader(tsvfile, delimiter='\t', quotechar='|')
      for rec in greader:
         if index != 'symbol':
            transcript_id = rec['refseq_mrna']
            if index == 'ensGene_transcript':
               transcript_id = rec['ensembl_transcript_id']
            if index == 'refGene_transcript' and not (transcript_id.startswith('NM_') or transcript_id.startswith('NR_')):
               continue
            if index == 'ensGene_transcript' and not transcript_id.startswith('ENST'):
               continue
            if not gene_xref.has_key(transcript_id):
               gene_xref[transcript_id] = []
            gene_xref[transcript_id].append(rec)
         else:
            if not rec['gene_biotype'] == 'LRG_gene':
               if not gene_xref.has_key(rec['symbol']):
                  gene_xref[rec['symbol']] = []
               gene_xref[rec['symbol']].append(rec)
               if not gene_xref.has_key(rec['ensembl_gene_id']):
                  gene_xref[rec['ensembl_gene_id']] = []
               gene_xref[rec['ensembl_gene_id']].append(rec)
   return gene_xref
