#!/usr/bin/env python

import os,re
import csv
import gzip

from pcgr.annoutils import assign_cds_exon_intron_annotations
from pcgr import pcgr_vars
from pcgr.utils import getlogger, check_file_exists, get_perl_exports, get_maxentscan_dir
from importlib.resources import files

def get_vep_command(file_paths, conf_options, input_vcf, output_vcf, debug = False):

    output_vcf_gz = f'{output_vcf}.gz'
    genome_assembly = conf_options['genome_assembly']

    vep_dir = file_paths['vep_dir']
    fasta_assembly = os.path.join(
        file_paths['refdata_assembly_dir'],
        'misc','fasta','assembly',
        f'Homo_sapiens.{pcgr_vars.VEP_ASSEMBLY[genome_assembly]}.dna.primary_assembly.fa.gz')

    logger = getlogger('check-fasta-files')
    check_file_exists(fasta_assembly, logger = logger)

    # for logging only
    plugins_in_use = "NearestExonJB, MaxEntScan"

    # List all VEP flags used when calling VEP
    vep_flags = (
        f'--hgvs --af_gnomad --variant_class --domains --symbol --protein --ccds --mane '
        f'--uniprot --appris --biotype --tsl --canonical --format vcf --cache --numbers '
        f'--total_length --allele_number --failed 1 --no_stats --no_escape --xref_refseq --vcf '
        f'--check_ref --dont_skip --flag_pick_allele_gene --plugin NearestExonJB,max_range=50000 '
        f'--plugin MaxEntScan,{get_maxentscan_dir()} '
        f'--force_overwrite --species homo_sapiens --offline')
    vep_options = (
        f'--dir {vep_dir} --assembly {pcgr_vars.VEP_ASSEMBLY[genome_assembly]} --cache_version {pcgr_vars.VEP_VERSION} '
        f'--fasta {fasta_assembly} '
        f'--pick_order {conf_options["conf"]["vep"]["vep_pick_order"]} '
        f'--buffer_size {conf_options["conf"]["vep"]["vep_buffer_size"]} '
        f'--fork {conf_options["conf"]["vep"]["vep_n_forks"]} '
        f'{vep_flags} '
        f'{"--verbose" if debug else "--quiet"}')

    gencode_set_in_use = "GENCODE - all transcripts"
    if conf_options['conf']['vep']['vep_no_intergenic'] == 1:
        vep_options += ' --no_intergenic'
    if conf_options['conf']['vep']['vep_regulatory'] == 1:
        vep_options += ' --regulatory'
    if conf_options['conf']['vep']['vep_gencode_basic'] == 1:
        vep_options += ' --gencode_basic'
        gencode_set_in_use = "GENCODE - basic transcript set (--gencode_basic)"

    # Compose full VEP command
    vep_main_command = f'{get_perl_exports()} && vep --input_file {input_vcf} --output_file {output_vcf} {vep_options}'
    vep_bgzip_command = f'bgzip -f -c {output_vcf} > {output_vcf_gz}'
    vep_tabix_command = f'tabix -f -p vcf {output_vcf_gz}'
    if debug:
        print(vep_main_command)

    vep_cmd = {}
    vep_cmd['main'] = vep_main_command
    vep_cmd['bgzip'] = vep_bgzip_command
    vep_cmd['tabix'] = vep_tabix_command
    vep_cmd['gencode_set_in_use'] = gencode_set_in_use
    vep_cmd['plugins_in_use'] = plugins_in_use
    vep_cmd['fasta_assembly'] = fasta_assembly

    return(vep_cmd)


def get_csq_record_annotations(csq_fields, varkey, logger, vep_csq_fields_map, transcript_xref_map, grantham_scores = None):
    """
    Generates a dictionary object containing the annotations of a CSQ record.

    Parameters:
    - csq_fields (list): A list of CSQ fields.
    - varkey (str): The VARKEY value.
    - logger (Logger): The logger object.
    - vep_csq_fields_map (dict): A dictionary mapping VEP CSQ fields to their indices.
    - transcript_xref_map (dict): A dictionary mapping Ensembl transcript IDs to their annotations.

    Returns:
    - csq_record (dict): A dictionary object containing the annotations of a CSQ record.
    """

    j = 0
    csq_record = {}

    csq_record['VARKEY'] = varkey
    ensembl_transcript_id = '.'

    # loop over block annotation elements (separated with '|'), and assign their values to the csq_record dictionary object
    while (j < len(csq_fields)):
        if j in vep_csq_fields_map['index2field']:

            ## consider non-empty CSQ fields
            if csq_fields[j] != '':
                csq_record[vep_csq_fields_map['index2field'][j]] = str(csq_fields[j])
                if vep_csq_fields_map['index2field'][j] == 'Feature':

                    ensembl_transcript_id = str(csq_fields[j])
                    if ensembl_transcript_id in transcript_xref_map:
                        for annotation in transcript_xref_map[ensembl_transcript_id].keys():
                            if annotation != 'SYMBOL':
                                ## assign additional (non-VEP provided) gene/transcript annotations from the custom
                                ## transcript_xref_map as key,value pairs in the csq_record object
                                csq_record[annotation] = transcript_xref_map[ensembl_transcript_id][annotation]
                    else:
                        if re.match(r'ENST', ensembl_transcript_id):
                            logger.warning(
                                'Could not find transcript xrefs for ' + str(ensembl_transcript_id))

                # Specifically assign PFAM protein domain as a csq_record key
                if vep_csq_fields_map['index2field'][j] == 'DOMAINS':
                    domain_identifiers = str(csq_fields[j]).split('&')
                    for v in domain_identifiers:
                        if v.startswith('Pfam'):
                            csq_record['PFAM_DOMAIN'] = str(re.sub(r'\.[0-9]{1,}$', '', re.sub(r'Pfam:', '', v)))

                # Assign COSMIC/DBSNP mutation ID's as individual key,value pairs in the csq_record object
                if vep_csq_fields_map['index2field'][j] == 'Existing_variation':
                    var_identifiers = str(csq_fields[j]).split('&')
                    parsed_identifiers = {'COSMIC_ID':[], 'DBSNP_RSID':[]}
                    for v in var_identifiers:
                        if v.startswith('COSV') or v.startswith('COSM'):
                            parsed_identifiers['COSMIC_ID'].append(v)
                        if v.startswith('rs'):
                            parsed_identifiers['DBSNP_RSID'].append(v)
                    for db in parsed_identifiers.keys():
                        if len(parsed_identifiers[db]) > 0:
                            csq_record[db] = '&'.join(parsed_identifiers[db])

                ## Sort (potentially multiple) variant consequence elements from VEP (they appear unsorted in some cases)
                ## Example: intron_variant&splice_region_variant
                if vep_csq_fields_map['index2field'][j] == "Consequence":
                    consequence_elements = sorted(str(csq_fields[j]).split('&'))
                    csq_record['Consequence'] = '&'.join(consequence_elements)
                
                if vep_csq_fields_map['index2field'][j] == "MaxEntScan_diff":
                    csq_record['MaxEntScan_diff'] = float(csq_fields[j])

            else:
                csq_record[vep_csq_fields_map['index2field'][j]] = None
        j = j + 1

    ## if VEP/Ensembl does not provide a symbol, use symbol provided by PCGR/CPSR gene_transcript_xref map
    if csq_record['SYMBOL'] is None and ensembl_transcript_id != ".":
        if ensembl_transcript_id in transcript_xref_map:
            csq_record['SYMBOL'] = transcript_xref_map[ensembl_transcript_id]['SYMBOL']

    # Assign coding status, protein change, grantham distance, coding sequence change, last exon/intron status etc as VCF info tags
    csq_record = assign_cds_exon_intron_annotations(csq_record, grantham_scores, logger)

    return(csq_record)


def pick_single_gene_csq(vep_csq_results,
                         pick_criteria_ordered = "mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length",
                         logger = None,
                         debug = 0):


    csq_candidates = []

    for csq_elem in vep_csq_results['picked_gene_csq']:
        if csq_elem is None:
            continue
        csq_candidate = {}

        ## default values (undefined properties)
        csq_candidate['mane_select'] = 1
        csq_candidate['mane_plus_clinical'] = 1
        csq_candidate['canonical'] = 1
        csq_candidate['appris'] = 18
        csq_candidate['biotype'] = 1
        csq_candidate['tsl'] = 6
        csq_candidate['ccds'] = 1
        csq_candidate['rank'] = 42

        ## set to picked as default
        csq_candidate['PICKED'] = True
        csq_candidate['varkey'] = csq_elem['VARKEY']
        csq_candidate['conskey'] = str(csq_elem['SYMBOL']) + ':' + str(csq_elem['Consequence'])

        ## MANE select status - lower value prioritized
        if csq_elem['MANE_SELECT'] is not None:
            csq_candidate['mane_select'] = 0

        ## MANE PLUS clnical status - lower value prioritized
        if csq_elem['MANE_PLUS_CLINICAL'] is not None:
            csq_candidate['mane_plus_clinical'] = 0

        ## CANONICAL status - lower value prioritized
        if csq_elem['CANONICAL'] is not None:
            if csq_elem['CANONICAL'] == 'YES':
                csq_candidate['canonical'] = 0

        ## APPRIS level - lower value prioritized
        if csq_elem['APPRIS'] is not None:
            if 'ALTERNATIVE' not in csq_elem['APPRIS']:
                csq_candidate['appris'] = int(re.sub(r'[A-Z]{1,}:?', '', csq_elem['APPRIS']))
            else:
                csq_candidate['appris'] = int(re.sub(r'ALTERNATIVE:','', csq_elem['APPRIS'])) + 10

        ## Biotype - lower value prioritized
        if csq_elem['BIOTYPE'] is not None:
            if csq_elem['BIOTYPE'] == 'protein_coding':
                csq_candidate['biotype'] = 0

        ## CCDS - lower value prioritized
        if csq_elem['CCDS'] is not None:
            csq_candidate['ccds'] = 0

        ## Consequence rank - lower value prioritized
        if csq_elem['Consequence'] is not None:
            main_cons = csq_elem['Consequence'].split('&')[0]
            if main_cons in pcgr_vars.VEP_consequence_rank:
                csq_candidate['rank'] = int(pcgr_vars.VEP_consequence_rank[main_cons])
            else:
                warn_msg = f"Missing Consequence in pcgr_vars.VEP_consequence_rank: {csq_elem['Consequence']} -  '{main_cons}'"
                if logger is not None:
                    logger.warn(warn_msg)

        ## TSL - lower value prioritized
        if csq_elem['TSL'] is not None:
            csq_candidate['tsl'] = int(csq_elem['TSL'])

        csq_candidates.append(csq_candidate)

    ## Go through pick criteria in pre-defined order
    ## - set 'PICKED' = False for all csq elements with a score above the minimum value for a given criterion
    ## - when there is only one element with 'PICKED' equal to TRUE, break out of the loop, and report the chosen transcript CSQ element
    chosen_csq_index = 0
    for rank_criterion in pick_criteria_ordered.split(','):
        if rank_criterion != 'length':
            lowest_score = 100
            i = 0
            for candidate in csq_candidates:
                if candidate[rank_criterion] <= lowest_score:
                    lowest_score = candidate[rank_criterion]
                i = i + 1

            for candidate in csq_candidates:
                if candidate[rank_criterion] > lowest_score:
                    candidate['PICKED'] = False

            j = 0
            num_picked = 0
            for candidate in csq_candidates:
                if candidate['PICKED'] is True:
                    num_picked += 1
                    chosen_csq_index = j
                j = j + 1

            if num_picked == 1:
                # if debug:
                #     print()
                #     for c in csq_candidates:
                #         all_rank_criterions = []
                #         all_rank_criterions.append('PICKED:' + str(c['PICKED']))
                #         all_rank_criterions.append(c['varkey'])
                #         all_rank_criterions.append(c['conskey'])
                #         for rank_criterion in pick_criteria_ordered.split(','):
                #             if rank_criterion in c:
                #                 all_rank_criterions.append(rank_criterion + ':' + str(c[rank_criterion]))

                #         rank_str = ' - '.join(map(str, all_rank_criterions))
                #         print(rank_str)
                #     print()
                break

    return(chosen_csq_index)

def parse_vep_csq(rec, transcript_xref_map, vep_csq_fields_map, grantham_scores, vep_pick_order, logger, pick_only=True,
                  csq_identifier='CSQ', debug = 0, targets_entrez_gene = None):

    """
    Function that parses the comma-separated CSQ elements found in the rec.INFO object (VCF)
    - creates an individual CSQ record for all transcript-specific elements provided as comma-separated elements in the CSQ tag
    - each individual record is gathered as a dictionary of properties (defined by vep_csq_field_map), i.e.
    - 'CSQ=A|missense_variant|KRAS++' in the VCF INFO element gives csq_record['Consequence'] = 'missense_variant',
       csq_record['SYMBOL'] = 'KRAS' etc.
    - if argument 'pick_only' is TRUE, only elements with 'PICK' == 1' is chosen
    """

    all_csq = []
    primary_csq_pick = []
    secondary_csq_pick = []
    all_transcript_consequences = []


    varkey = str(rec.CHROM) + '_' + str(rec.POS) + '_' + str(rec.REF) + '_' + str(','.join(rec.ALT))

    #var_gwas_info = rec.INFO.get('GWAS_HIT')
    #gwas_hit = False
    #if var_gwas_info is not None:
    #    gwas_hit = True

    #found_in_target = 0

    ## Retrieve the INFO element provided by VEP (default 'CSQ') in the VCF object, and
    ## loop through all transcript-specific consequence blocks provided, e.g.
    #  CSQ=A|intron_variant|||.., A|splice_region_variant|||, and so on.
    for csq in rec.INFO.get(csq_identifier).split(','):
        csq_fields = csq.split('|')
        entrezgene = '.'

        ## Entrez gene identifier is not provided by VEP, pull out this from 'transcript_xref_map' for a given
        ## transcript-specific CSQ block
        ##  - used for 'consequence_entry' object that are added to 'vep_all_csq' array
        k = 0
        while(k < len(csq_fields)):
            if k in vep_csq_fields_map['index2field']:
                if vep_csq_fields_map['index2field'][k] == 'Feature':
                    ensembl_transcript_id = csq_fields[k]
                    if ensembl_transcript_id != '' and ensembl_transcript_id.startswith('ENST'):
                        if ensembl_transcript_id in transcript_xref_map.keys():
                            if 'ENTREZGENE' in transcript_xref_map[ensembl_transcript_id].keys():
                                entrezgene = transcript_xref_map[ensembl_transcript_id]['ENTREZGENE']
            k = k + 1


        ## CPSR - consider all picked gene-specific consequences (considering that a variant may occasionally
        ## overlap other, non-CPSR targets)
        if pick_only is False:
            csq_record = get_csq_record_annotations(csq_fields, varkey, logger, vep_csq_fields_map, transcript_xref_map, grantham_scores)
            all_csq.append(csq_record)
            if csq_record['PICK'] == '1':
                if 'Feature_type' in csq_record:
                    if csq_record['Feature_type'] == 'RegulatoryFeature':
                        primary_csq_pick.append(csq_record)
                        secondary_csq_pick.append(csq_record)
                    if csq_record['Feature_type'] == 'Transcript':
                        ## Entrez gene identifier is provided by VEP
                        if 'ENTREZGENE' in csq_record.keys():
                            ## Consequence in CPSR target - append to primary_csq_pick
                            if csq_record['ENTREZGENE'] in targets_entrez_gene.keys():
                                primary_csq_pick.append(csq_record)
                                ## Consequence not in CPSR targets, nor with Entrez identifier
                                # (pseudogenes etc) - append to secondary_csq_pick
                            else:
                                secondary_csq_pick.append(csq_record)

                        ## Entrez gene identifier is not provided by VEP
                        else:
                            if entrezgene != "." and entrezgene in targets_entrez_gene.keys():
                                primary_csq_pick.append(csq_record)
                            else:
                                ## Consequence not in CPSR targets, nor with Entrez identifier
                                # (pseudogenes etc) - append to secondary_csq_pick
                                warning = 1
                                secondary_csq_pick.append(csq_record)
                    else:
                        ## intergenic
                        if csq_record['Feature_type'] is None:
                            secondary_csq_pick.append(csq_record)

        ## PCGR - consider picked transcript consequence only
        else:
            # loop over VEP consequence blocks PICK'ed according to VEP's ranking scheme
            # only consider the primary/picked consequence when expanding with annotation tags

            if csq_fields[vep_csq_fields_map['field2index']['PICK']] == "1":
                csq_record = get_csq_record_annotations(
                    csq_fields, varkey, logger, vep_csq_fields_map, transcript_xref_map, grantham_scores)
                # Append transcript consequence to primary_csq_pick
                primary_csq_pick.append(csq_record)
        symbol = "."
        hgvsc = "."
        hgvsp = "."
        exon = "."
        if csq_fields[vep_csq_fields_map['field2index']['EXON']] != "":
            if "/" in csq_fields[vep_csq_fields_map['field2index']['EXON']]:
                exon = str(csq_fields[vep_csq_fields_map['field2index']['EXON']].split('/')[0])
        if csq_fields[vep_csq_fields_map['field2index']['SYMBOL']] != "":
            symbol = str(csq_fields[vep_csq_fields_map['field2index']['SYMBOL']])
        if csq_fields[vep_csq_fields_map['field2index']['HGVSc']] != "":
            hgvsc = str(csq_fields[vep_csq_fields_map['field2index']['HGVSc']].split(':')[1])
        else:
            if len(primary_csq_pick) == 1:
                if 'HGVSc' in primary_csq_pick[0]:
                    if primary_csq_pick[0]['HGVSc'] is not None:
                        if ':' in primary_csq_pick[0]['HGVSc']:
                            hgvsc = str(primary_csq_pick[0]['HGVSc'].split(':')[1])
        if csq_fields[vep_csq_fields_map['field2index']['HGVSp']] != "":
            hgvsp = str(csq_fields[vep_csq_fields_map['field2index']['HGVSp']].split(':')[1])
        consequence_entry = (str(csq_fields[vep_csq_fields_map['field2index']['Consequence']]) + ":" +
            str(symbol) + ":" +
            str(entrezgene) + ":" +
            str(hgvsc) + ":" +
            str(hgvsp) + ":" +
            str(exon) + ":" +
            str(csq_fields[vep_csq_fields_map['field2index']['Feature_type']]) + ":" +
            str(csq_fields[vep_csq_fields_map['field2index']['Feature']]) + ":" +
            str(csq_fields[vep_csq_fields_map['field2index']['BIOTYPE']]))
        all_transcript_consequences.append(consequence_entry)


    ## CPSR - consider all picked VEP blocks
    if pick_only is False:
        if len(primary_csq_pick) == 0:
            # variant does not overlap virtual panel genes, but found with
            # other targeted regions (GWAS hits etc) - choose secondary csq pick (if non-empty)
            if len(secondary_csq_pick) > 0:
                primary_csq_pick = secondary_csq_pick
            else:
                ## No primary picked VEP block, no secondary picked VEP block -
                # consider all VEP blocks (does it ever happen?)
                primary_csq_pick = all_csq
                logger.warning(f"NO VEP BLOCK PICKED for {varkey}")
        else:
            ## don't consider regulatory consequence as single chosen consequence
            if len(primary_csq_pick) == 1 and len(secondary_csq_pick) > 0:
                if 'Feature_type' in primary_csq_pick[0]:
                    if primary_csq_pick[0]['Feature_type'] == 'RegulatoryFeature':
                        primary_csq_pick = secondary_csq_pick

    vep_csq_results = {}
    vep_csq_results['picked_gene_csq'] = primary_csq_pick
    vep_csq_results['all_csq'] = all_transcript_consequences
    vep_csq_results['picked_csq'] = None
    vep_chosen_csq_idx = 0

    ## If multiple transcript-specific variant consequences highlighted by --pick_allele_gene,
    ## prioritize/choose block of consequence according to 'vep_pick_order'
    if len(vep_csq_results['picked_gene_csq']) > 1:
        vep_chosen_csq_idx = pick_single_gene_csq(
            vep_csq_results, pick_criteria_ordered = vep_pick_order, logger = logger, debug = debug)
        vep_csq_results['picked_csq'] = vep_csq_results['picked_gene_csq'][vep_chosen_csq_idx]
    else:
        ## check that size if 1, otherwise prompt error below
        #logger.info('ERROR: No VEP block chosen by --pick_allele_gene')
        if len(vep_csq_results['picked_gene_csq']) == 1:
            vep_csq_results['picked_csq'] = vep_csq_results['picked_gene_csq'][vep_chosen_csq_idx]
        else:
            logger.error('ERROR: No VEP block chosen by --pick_allele_gene')



    return (vep_csq_results)

