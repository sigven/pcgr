#!/usr/bin/env python

import os
import re
import sys
import csv
import cyvcf2
import gzip

from cyvcf2 import VCF, Writer
from pcgr import utils
from pcgr.utils import check_subprocess, error_message

csv.field_size_limit(500 * 1024 * 1024)
threeLettertoOneLetterAA = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Gly': 'G', 'His': 'H',
                            'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': 'X'}

nuclear_chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']

def read_infotag_file(vcf_info_tags_tsv, scope = "vep"):
    """
    Function that reads a file that lists VCF INFO tags produced by PCGR/CPSR/gvanno.
    An example line of the VCF info tag file is the following:

    tag	number	type	description category
    Consequence	.	String	"Impact modifier for the consequence type (picked by VEP's --flag_pick_allele option)."   vep

    A dictionary is returned, with the tag as the key, and the full dictionary record as the value
    """
    info_tag_xref = {}  # dictionary returned
    if not os.path.exists(vcf_info_tags_tsv):
        ## issue warning
        return info_tag_xref
    tsvfile = open(vcf_info_tags_tsv, 'r')
    reader = csv.DictReader(tsvfile, delimiter='\t')
    for row in reader:
        if not row['tag'] in info_tag_xref:
            if scope in row['category']:
                info_tag_xref[row['tag']] = row

    return info_tag_xref

def read_vcfanno_tag_file(vcfanno_tag_file, logger):

    infotag_results = {}

    if not os.path.exists(vcfanno_tag_file):
        err_msg = f"vcfanno configuration file ({vcfanno_tag_file}) does not exist"
        error_message(err_msg, logger)
    
    f = open(vcfanno_tag_file, "r")
    info_elements = f.readlines()
    f.close()


    for i in info_elements:
        elem_info_content = i.rstrip().split(',')
        row = {}
        for elem in elem_info_content:
            if elem.startswith('Type'):
                row['type'] = elem.replace('Type=', '')
            if elem.startswith('##INFO'):
                row['tag'] = elem.replace('##INFO=<ID=', '')
            if elem.startswith('Number'):
                row['number'] = elem.replace('Number=', '')
            if elem.startswith('Description'):
                row['description'] = re.sub(r'">|Description="', '', elem)        
            
        if 'tag' in row and not row['tag'] in infotag_results:
            infotag_results[row['tag']] = row
    
    return(infotag_results)
        
    
def read_genexref_namemap(gene_xref_namemap_tsv, logger):
    """
    Function that reads a file that lists names of tags in GENE_TRANSCRIPT_XREF annotation.
    """

    namemap_xref = {}  # dictionary returned
    if not os.path.exists(gene_xref_namemap_tsv):
        err_msg = f"gene_transcript_xref BED mapping file ({gene_xref_namemap_tsv}) does not exist"
        error_message(err_msg, logger)
    
    with gzip.open(gene_xref_namemap_tsv, mode='rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            namemap_xref[row['name']] = int(row['index'])
    f.close()

    return namemap_xref


def write_pass_vcf(annotated_vcf, logger):
    """
    Function prints all PASS variants from a VCF file. New VCF file appends '.pass.' to the filename.
    """
    #out_vcf = re.sub(r'\.annotated\.vcf\.gz$','.annotated.pass.vcf',annotated_vcf)
    out_vcf = re.sub(r'\.vcf\.gz$', '.pass.vcf', annotated_vcf)
    vcf = VCF(annotated_vcf)
    w = Writer(out_vcf, vcf)

    num_rejected = 0
    num_pass = 0
    for rec in vcf:
        if rec.FILTER is None or rec.FILTER == 'None':
            w.write_record(rec)
            num_pass += 1
        else:
            num_rejected += 1

    vcf.close()
    w.close()

    logger.info('Number of non-PASS/REJECTED variant calls: ' +
                str(num_rejected))
    logger.info('Number of PASSed variant calls: ' + str(num_pass))
    if num_pass == 0:
        logger.warning(
            'There are zero variants with a \'PASS\' filter in the VCF file')
        os.system('bgzip -dc ' + str(annotated_vcf) +
                  ' egrep \'^#\' > ' + str(out_vcf))
    # else:
    os.system('bgzip -f ' + str(out_vcf))
    os.system('tabix -f -p vcf ' + str(out_vcf) + '.gz')

    return


def map_regulatory_variant_annotations(vep_csq_records):
    """
    Function that considers an array of VEP CSQ records and appends all regulatory variant consequent annotations (open chromatin, TF_binding_site,
    CTCF_binding_site, promoter (flanks), enhancers ) into a single comma-separated string. Each individual regulatory annotation is formatted as:
    <Consequence>|<Feature_type>|<Feature>|<BIOTYPE>|<MOTIF_NAME>|<MOTIF_POS>|<HIGH_INF_POS>|<MOTIF_SCORE_CHANGE>|<TRANSCRIPTION_FACTORS>
    """

    regulatory_annotation = '.'
    if len(vep_csq_records) == 1:
        return regulatory_annotation

    j = 0
    regulatory_annotations = []
    while j < len(vep_csq_records):
        missing_key_annotation = False
        for k in ['Feature_type', 'Consequence', 'Feature']:
            if not k in vep_csq_records[j].keys():
                missing_key_annotation = True

        if missing_key_annotation is False:

            # RegulatoryFeature annotations - open chromatin, promoters (flanks), enhancers, CTCF binding sites
            if vep_csq_records[j]['Feature_type'] == 'RegulatoryFeature':
                biotype = ""
                if re.match(r"^(enhancer|promoter|open|CTCF|TF_)", vep_csq_records[j]['BIOTYPE']):
                    biotype = vep_csq_records[j]['BIOTYPE']

                annotation = str(vep_csq_records[j]['Consequence']) + '|' + \
                    str(vep_csq_records[j]['Feature_type']) + '|' + \
                    str(vep_csq_records[j]['Feature']) + '|' + \
                    str(biotype) + '|||||'

                regulatory_annotations.append(annotation)

            # MotifFeature annotations (TF)
            if vep_csq_records[j]['Feature_type'] == 'MotifFeature':
                missing_motif_annotation = False
                for annotation in ['MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'TRANSCRIPTION_FACTORS']:
                    if not annotation in vep_csq_records[j].keys():
                        missing_motif_annotation = True

                if missing_motif_annotation is False:
                    annotation = str(vep_csq_records[j]['Consequence']) + '|' + \
                        str(vep_csq_records[j]['Feature_type']) + '|' + \
                        str(vep_csq_records[j]['Feature']) + '|TF_binding_site|' + \
                        str(vep_csq_records[j]['MOTIF_NAME']) + '|' + \
                        str(vep_csq_records[j]['MOTIF_POS']) + '|' + \
                        str(vep_csq_records[j]['HIGH_INF_POS']) + '|' + \
                        str(vep_csq_records[j]['MOTIF_SCORE_CHANGE']) + '|' + \
                        str(vep_csq_records[j]['TRANSCRIPTION_FACTORS'])

                    regulatory_annotations.append(annotation)

        j = j + 1

    if len(regulatory_annotations) > 0:
        regulatory_annotation = ','.join(regulatory_annotations)

    return (regulatory_annotation)


def get_correct_cpg_transcript(vep_csq_records):
    """
    Function that considers an array of VEP CSQ records and picks most relevant consequence (and gene) from
    neighbouring genes/transcripts of relevance for cancer predisposition (cpg = cancer predisposition gene)
    """

    csq_idx = 0
    if len(vep_csq_records) == 1:
        return csq_idx

    # some variants iare assigned multiple transcript consequences
    # if cancer predisposition genes are in the vicinity of other genes, choose the cancer predisposition gene
    # if there are neighbouring cancer-predispositon genes, choose custom gene, preferring coding change (see below, KLLN/PTEN, XPC/TMEM43, NTHL1/TSC2)
    csq_idx_dict = {}
    for g in ['KLLN', 'PTEN', 'XPC', 'TMEM43', 'NTHL1', 'TSC2']:
        csq_idx_dict[g] = {}
        csq_idx_dict[g]['idx'] = -1
        csq_idx_dict[g]['coding'] = False

    j = 0
    while j < len(vep_csq_records):
        if 'CANCER_PREDISPOSITION_SOURCE' in vep_csq_records[j].keys() or 'GE_PANEL_ID' in vep_csq_records[j].keys():
            csq_idx = j
            if 'SYMBOL' in vep_csq_records[j].keys():
                if vep_csq_records[j]['SYMBOL'] in csq_idx_dict.keys():
                    csq_idx_dict[str(vep_csq_records[j]['SYMBOL'])]['idx'] = j
                    if vep_csq_records[j]['CODING_STATUS'] == 'coding':
                        csq_idx = j  # prefer coding on over anything else
                        csq_idx_dict[str(vep_csq_records[j]
                                         ['SYMBOL'])]['coding'] = True
        j = j + 1

    if csq_idx_dict['KLLN']['idx'] != -1 and csq_idx_dict['PTEN']['idx'] != -1:
        csq_idx = csq_idx_dict['PTEN']['idx']
        if csq_idx_dict['KLLN']['coding'] is True:
            csq_idx = csq_idx_dict['KLLN']['idx']

    if csq_idx_dict['XPC']['idx'] != -1 and csq_idx_dict['TMEM43']['idx'] != -1:
        csq_idx = csq_idx_dict['XPC']['idx']
        if csq_idx_dict['TMEM43']['coding'] is True:
            csq_idx = csq_idx_dict['TMEM43']['idx']

    if csq_idx_dict['TSC2']['idx'] != -1 and csq_idx_dict['NTHL1']['idx'] != -1:
        csq_idx = csq_idx_dict['TSC2']['idx']
        if csq_idx_dict['NTHL1']['coding'] is True:
            csq_idx = csq_idx_dict['NTHL1']['idx']

    if csq_idx is None:
        csq_idx = 0
    return csq_idx


def threeToOneAA(aa_change):

    for three_letter_aa in threeLettertoOneLetterAA.keys():
        aa_change = aa_change.replace(
            three_letter_aa, threeLettertoOneLetterAA[three_letter_aa])

    aa_change = re.sub(r'[A-Z]{1}fsX([0-9]{1,}|\?)', 'fs', aa_change)
    return aa_change


def assign_cds_exon_intron_annotations(csq_record):

    csq_record['CODING_STATUS'] = 'noncoding'
    csq_record['EXONIC_STATUS'] = 'nonexonic'
    csq_record['SPLICE_DONOR_RELEVANT'] = False
    csq_record['NULL_VARIANT'] = False
    csq_record['INTRON_POSITION'] = 0
    csq_record['EXON_POSITION'] = 0
    csq_record['CDS_CHANGE'] = '.'
    csq_record['HGVSp_short'] = '.'
    csq_record['PROTEIN_CHANGE'] = '.'
    csq_record['EXON_AFFECTED'] = '.'
    csq_record['LOSS_OF_FUNCTION'] = False

    coding_csq_pattern = r"(stop_(lost|gained)|start_lost|frameshift_|missense_|splice_(donor|acceptor)|protein_altering|inframe_)"
    wes_csq_pattern = r"(stop_(lost|gained)|start_lost|frameshift_|missense_|splice_(donor|acceptor)|protein_altering|inframe_|synonymous|(start|stop)_retained)"
    null_pattern = r"(stop_gained|frameshift_)"
    splice_region_pattern = r"(splice_|intron_variant)"
    splice_donor_pattern = r"(splice_region_variant|splice_donor_variant|splice_donor_region_variant|splice_donor_5th_base_variant)"
    if re.search(coding_csq_pattern, str(csq_record['Consequence'])) is not None:
        csq_record['CODING_STATUS'] = 'coding'

    if re.search(wes_csq_pattern, str(csq_record['Consequence'])) is not None:
        csq_record['EXONIC_STATUS'] = 'exonic'

    if 'LoF' in csq_record:
        if csq_record['LoF'] != '.':
            csq_record['LOSS_OF_FUNCTION'] = True
            
        ## Don't list LoF as True if consequence is assigned as missense
        if re.search(r"^missense", csq_record['Consequence']) is not None:
            if csq_record['LOSS_OF_FUNCTION'] is True:
                csq_record['LOSS_OF_FUNCTION'] = False
        

    if re.search(null_pattern, str(csq_record['Consequence'])) is not None:
        csq_record['NULL_VARIANT'] = True
    
    if re.search(splice_donor_pattern, str(csq_record['Consequence'])) is not None and re.search(r'(\+3(A|G)>|\+4A>|\+5G>)', str(csq_record['HGVSc'])) is not None:
        csq_record['SPLICE_DONOR_RELEVANT'] = True

    if re.search(splice_region_pattern, str(csq_record['Consequence'])) is not None:
        match = re.search(
            r"((-|\+)[0-9]{1,}(dup|del|inv|((ins|del|dup|inv|delins)(A|G|C|T){1,})|(A|C|T|G){1,}>(A|G|C|T){1,}))$", str(csq_record['HGVSc']))
        if match is not None:
            pos = re.sub(
                r"(\+|dup|del|delins|ins|inv|(A|G|C|T){1,}|>)", "", match.group(0))
            if utils.is_integer(pos):
                csq_record['INTRON_POSITION'] = int(pos)

    if 'NearestExonJB' in csq_record.keys():
        if not csq_record['NearestExonJB'] is None:
            if re.search(r"synonymous_|missense_|stop_|inframe_|start_", str(csq_record['Consequence'])) is not None and str(csq_record['NearestExonJB']) != "":
                exon_pos_info = csq_record['NearestExonJB'].split("+")
                if len(exon_pos_info) == 4:
                    if utils.is_integer(exon_pos_info[1]) and str(exon_pos_info[2]) == "end":
                        csq_record['EXON_POSITION'] = -int(exon_pos_info[1])
                    if utils.is_integer(exon_pos_info[1]) and str(exon_pos_info[2]) == "start":
                        csq_record['EXON_POSITION'] = int(exon_pos_info[1])

    if not csq_record['HGVSc'] is None:
        if csq_record['HGVSc'] != '.':
            if 'splice_acceptor_variant' in csq_record['Consequence'] or 'splice_donor_variant' in csq_record['Consequence'] \
                or 'splice_donor_5th_base_variant' in csq_record['Consequence'] or 'splice_region_variant' in csq_record['Consequence'] \
                    or 'splice_polypyrimidine_tract_variant' in csq_record['Consequence']:
                key = str(csq_record['Consequence']) + \
                    ':' + str(csq_record['HGVSc'])
                csq_record['CDS_CHANGE'] = key

    if csq_record['Amino_acids'] is None or csq_record['Protein_position'] is None or csq_record['Consequence'] is None:
        return
    
    if not csq_record['Protein_position'] is None:
        if csq_record['Protein_position'].startswith('-'):
            return

    protein_change = '.'
    if '/' in csq_record['Protein_position']:
        protein_position = str(csq_record['Protein_position'].split('/')[0])
        if '-' in protein_position:
            if protein_position.split('-')[0].isdigit():
                csq_record['AMINO_ACID_START'] = protein_position.split('-')[0]
            if protein_position.split('-')[1].isdigit():
                csq_record['AMINO_ACID_END'] = protein_position.split('-')[1]
        else:
            if protein_position.isdigit():
                csq_record['AMINO_ACID_START'] = protein_position
                csq_record['AMINO_ACID_END'] = protein_position

    if not csq_record['HGVSp'] is None:
        if csq_record['HGVSp'] != '.':
            if ':' in csq_record['HGVSp']:
                protein_identifier = str(csq_record['HGVSp'].split(':')[0])
                if protein_identifier.startswith('ENSP'):
                    protein_change_VEP = str(csq_record['HGVSp'].split(':')[1])
                    protein_change = threeToOneAA(protein_change_VEP)
                    csq_record['PROTEIN_CHANGE'] = protein_change_VEP

    if 'synonymous_variant' in csq_record['Consequence']:
        protein_change = 'p.' + \
            str(csq_record['Amino_acids']) + \
            str(protein_position) + str(csq_record['Amino_acids'])
        if 'stop_lost' in str(csq_record['Consequence']) and '/' in str(csq_record['Amino_acids']):
            protein_change = 'p.X' + \
                str(protein_position) + \
                str(csq_record['Amino_acids']).split('/')[1]

    csq_record['HGVSp_short'] = protein_change
    exon_number = 'NA'
    if not csq_record['EXON'] is None:
        if csq_record['EXON'] != '.':
            if '/' in csq_record['EXON']:
                exon_number = str(csq_record['EXON']).split('/')[0]
                csq_record['AFFECTED_EXON'] = exon_number
                if '-' in exon_number:
                    csq_record['AFFECTED_EXON'] = exon_number.split('-')[0]
                num_exons = str(csq_record['EXON']).split('/')[1]
                if exon_number == num_exons:
                    csq_record['LAST_EXON'] = True

    if not csq_record['INTRON'] is None:
        if csq_record['INTRON'] != '.':
            if '/' in csq_record['INTRON']:
                intron_number = str(csq_record['INTRON']).split('/')[0]
                num_introns = str(csq_record['INTRON']).split('/')[1]
                if intron_number == num_introns:
                    csq_record['LAST_INTRON'] = True

    if not csq_record['HGVSc'] is None:
        if csq_record['HGVSc'] != '.':
            if protein_change != '.':
                key = str(csq_record['Consequence']) + ':' + str(csq_record['HGVSc']) + ':exon' + str(exon_number) + ':' + str(protein_change)
                csq_record['CDS_CHANGE'] = key

    return


def make_transcript_xref_map(rec, fieldmap, xref_tag='GENE_TRANSCRIPT_XREF'):
    transcript_xref_map = {}
    if not rec.INFO.get(xref_tag) is None:
        for tref in rec.INFO.get(xref_tag).split(','):
            xrefs = tref.split('|')
            ensembl_transcript_id = str(xrefs[0])
            transcript_xref_map[ensembl_transcript_id] = {}
            for annotation in fieldmap.keys():
                annotation_index = fieldmap[annotation]
                if annotation_index > (len(xrefs) - 1):
                    continue
                if xrefs[annotation_index] != '':
                    transcript_xref_map[ensembl_transcript_id][annotation] = xrefs[annotation_index]

    return (transcript_xref_map)
