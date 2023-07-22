#!/usr/bin/env python

import argparse
import random
import re, os
import glob

from pcgr import utils
from pcgr.utils import check_subprocess
from pcgr import annoutils
from pcgr.annoutils import read_vcfanno_tag_file, get_vcf_info_tags, print_vcf_header


def __main__():
    parser = argparse.ArgumentParser(description='Run brentp/vcfanno - annotate a VCF file against multiple VCF files in parallel',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'query_vcf', help='Bgzipped input VCF file with query variants (SNVs/InDels)')
    parser.add_argument(
        'out_vcf', help='Output VCF file with appended annotations from multiple VCF files')
    parser.add_argument(
        'pcgr_db_dir', help='PCGR assembly-specific data directory')
    parser.add_argument(
        '--num_processes', help="Number of processes vcfanno can use during annotation", default=4)
    parser.add_argument("--clinvar", action="store_true",
                        help="Annotate VCF with annotations from ClinVar")
    # parser.add_argument("--ncer", action="store_true",
    #                     help="Annotate VCF with ranking of variant deleteriousness in non-coding regions (ncER)")
    parser.add_argument('--dbmts', action="store_true",
                        help="Annotate VCF file with variants predicted to cause loss/gain of miRNA target sites in 3'UTR regions")
    parser.add_argument('--gerp', action="store_true",
                        help="Annotate VCF file with GERP RS scores (cancer predisposition genes/actionable/secondary findings/GWAS loci only)")
    parser.add_argument("--dbnsfp", action="store_true",
                        help="Annotate VCF with annotations from database of non-synonymous functional predictions")
    parser.add_argument("--tcga", action="store_true",
                        help="Annotate VCF with variant frequencies from the The Cancer Genome Atlas")
    parser.add_argument("--gene_transcript_xref", action="store_true",
                        help="Annotate VCF with transcript annotations from PCGR (drug targets, actionable genes, cancer gene roles, etc)")
    parser.add_argument("--gwas", action="store_true",
                        help="Annotate VCF against moderate-to-low cancer risk variants, as identified from genome-wide association studies (GWAS)")
    parser.add_argument("--rmsk", action="store_true",
                        help="Annotate VCF against known sequence repeats, as identified by RepeatMasker (rmsk)")
    parser.add_argument("--simplerepeat", action="store_true",
                        help="Annotate VCF against known sequence repeats, as identified by Tandem Repeats Finder (simplerepeats)")
    parser.add_argument("--winmsk", action="store_true",
                        help="Annotate VCF against known sequence repeats, as identified by Windowmasker (winmsk)")
    parser.add_argument("--gnomad_non_cancer", action="store_true",
                        help="Annotate VCF with population-specific allele frequencies in gnomAD non-cancer subjects (cancer predisposition genes/actionable/secondary findings/GWAS loci only)")
    parser.add_argument("--pon_vcf", dest="pon_vcf",
                        help="Annotate VCF with calls from panel of normals (PON VCF)")
    parser.add_argument("--keep_logs", action="store_true")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="Print full commands to log, default: %(default)s")

    args = parser.parse_args()

    logger = utils.getlogger('pcgr-vcfanno')

    query_info_tags = get_vcf_info_tags(args.query_vcf)
    vcfheader_file = args.out_vcf + '.tmp.' + \
        str(random.randrange(0, 10000000)) + '.header.txt'
    conf_fname = args.out_vcf + '.tmp.conf.toml'
    print_vcf_header(args.query_vcf, vcfheader_file,
                     logger, chromline_only=False)
    
    vcfanno_tracks = {}
    ## BED
    vcfanno_tracks['rmsk'] = args.rmsk
    vcfanno_tracks['winmsk'] = args.winmsk
    vcfanno_tracks['simplerepeat'] = args.simplerepeat
    vcfanno_tracks['gene_transcript_xref'] = args.gene_transcript_xref

    ## VCF
    vcfanno_tracks['tcga'] = args.tcga
    vcfanno_tracks['gwas'] = args.gwas
    vcfanno_tracks['dbmts'] = args.dbmts
    vcfanno_tracks['dbnsfp'] = args.dbnsfp
    vcfanno_tracks['gerp'] = args.gerp
    vcfanno_tracks['clinvar'] = args.clinvar
    vcfanno_tracks['gnomad_non_cancer'] = args.gnomad_non_cancer

    run_vcfanno(args.num_processes, args.query_vcf, vcfanno_tracks, query_info_tags, vcfheader_file,
                args.pcgr_db_dir, conf_fname, args.pon_vcf, args.out_vcf, args.keep_logs, args.debug, logger)


def run_vcfanno(num_processes, query_vcf, vcfanno_tracks, query_info_tags, vcfheader_file, pcgr_db_dir, conf_fname,
                pon_vcf, output_vcf, keep_logs, debug, logger):

    """
    Function that annotates a VCF file with vcfanno against a user-defined set of germline and somatic VCF files
    """

    ## Collect metadata (VCF INFO tags) for annotations populated with vcfanno
    metadata_vcf_infotags = {}
    infotags = {}

    track_file_info = {}

    ## INFO tags used for vcfanno annotation
    track_file_info['tags_fname'] = {}

    ## Source file (VCF/BED) used for vcfanno annotation
    track_file_info['track_fname'] = {}
    
    for variant_track in ['clinvar','tcga','gwas','dbmts','dbnsfp','gnomad_non_cancer']:
        track_file_info['tags_fname'][variant_track] = os.path.join(pcgr_db_dir,'variant','vcf', variant_track, f'{variant_track}.vcfanno.vcf_info_tags.txt')
        track_file_info['track_fname'][variant_track] = os.path.join(pcgr_db_dir,'variant','vcf', variant_track, f'{variant_track}.vcf.gz')

    for bed_track in ['simplerepeat','winmsk','rmsk','gerp']:
        track_file_info['tags_fname'][bed_track] = os.path.join(pcgr_db_dir,'misc','bed', bed_track, f'{bed_track}.vcfanno.vcf_info_tags.txt')
        track_file_info['track_fname'][bed_track] = os.path.join(pcgr_db_dir,'misc','bed', bed_track, f'{bed_track}.bed.gz')

    track_file_info['tags_fname']['gene_transcript_xref'] = os.path.join(pcgr_db_dir,'gene','bed', 'gene_transcript_xref', 'gene_transcript_xref.vcfanno.vcf_info_tags.txt')
    track_file_info['track_fname']['gene_transcript_xref'] = os.path.join(pcgr_db_dir,'gene','bed', 'gene_transcript_xref', 'gene_transcript_xref.bed.gz')
    
    for track in track_file_info['tags_fname']:

        if not vcfanno_tracks[track] is True:
            continue

        infotags_vcfanno = read_vcfanno_tag_file(track_file_info['tags_fname'][track], logger)
        infotags[track] = infotags_vcfanno.keys()
        for tag in infotags_vcfanno:
            if tag in query_info_tags:
                logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the ' + str(
                    track) + ' VCF/BED annotation file. This tag will be overwritten if not renamed in the query VCF')
            metadata_vcf_infotags[tag] = infotags_vcfanno[tag]
        
        ## append track to vcfanno configuration file
        append_to_conf_file(track, infotags[track], track_file_info['track_fname'][track], conf_fname)

        ## Append VCF INFO tags to VCF header file
        check_subprocess(logger, f'cat {track_file_info["tags_fname"][track]} >> {vcfheader_file}', debug=False)
    
    
   
    panel_normal_tags = ["PANEL_OF_NORMALS"]
   
    if not pon_vcf is None:
        if "PANEL_OF_NORMALS" in query_info_tags:
            logger.warning(
                "Query VCF has INFO tag \"PANEL_OF_NORMALS\" - this is also present in the panel of normal VCF file. This tag will be overwritten if not renamed in the query VCF")
        
        vcf_info_tags_file = os.path.join(pcgr_db_dir,'variant','vcf', 'panel_of_normals', 'panel_of_normals.vcfanno.vcf_info_tags.txt')
        check_subprocess(
            logger, f'cat {vcf_info_tags_file} >> {vcfheader_file}', debug=False)
        fh = open(conf_fname, 'a')
        fh.write('[[annotation]]\n')
        fh.write('file="' + str(pon_vcf) + '"\n')
        fields_string = 'fields = ["' + '","'.join(panel_normal_tags) + '"]'
        ops = ['self'] * len(panel_normal_tags)
        ops_string = 'ops=["' + '","'.join(ops) + '"]'
        fh.write(fields_string + '\n')
        fh.write(ops_string + '\n\n')
        fh.close()

    out_vcf_vcfanno_unsorted1 = output_vcf + '.tmp.unsorted.1'
    query_prefix = re.sub(r'\.vcf.gz$', '', query_vcf)
    print_vcf_header(query_vcf, vcfheader_file, logger, chromline_only=True)
    
    vcfanno_command = f"vcfanno -p={num_processes} {conf_fname} {query_vcf} > {out_vcf_vcfanno_unsorted1} 2> {query_prefix}.vcfanno.log"
    if debug:
        logger.info(f"vcfanno command: {vcfanno_command}")
    check_subprocess(logger, vcfanno_command, debug)

    check_subprocess(
        logger, f'cat {vcfheader_file} > {output_vcf}', debug=False)
    check_subprocess(
        logger, f'cat {out_vcf_vcfanno_unsorted1} | grep -v \'^#\' >> {output_vcf}', debug=False)
    check_subprocess(logger, f'bgzip -f {output_vcf}', debug)
    check_subprocess(logger, f'tabix -f -p vcf {output_vcf}.gz', debug)
    if not keep_logs:
        for tmpf in glob.glob(f"{output_vcf}.tmp*"):
            utils.remove(tmpf)
        utils.remove(f"{query_prefix}.vcfanno.log")


def append_to_conf_file(datasource, datasource_info_tags, datasource_track_fname, conf_fname):
    """
    Function that appends data to a vcfanno conf file ('conf_fname') according to user-defined ('datasource').
    The datasource defines the set of tags that will be appended during annotation
    """
    fh = open(conf_fname, 'a')
    fh.write('[[annotation]]\n')
    fh.write('file="' + str(datasource_track_fname) + '"\n')
    if datasource == 'ncer' or datasource == 'gerp':        
        fh.write('columns=[4]\n')
        names_string = 'names=["' + '","'.join(datasource_info_tags) + '"]'
        fh.write(names_string + '\n')
        fh.write('ops=["mean"]\n\n')
    elif datasource == 'gene_transcript_xref' or datasource == 'rmsk':       
        fh.write('columns=[4]\n')
        names_string = 'names=["' + '","'.join(datasource_info_tags) + '"]'
        fh.write(names_string + '\n')
        fh.write('ops=["concat"]\n\n')
    elif datasource == 'simplerepeat' or datasource == 'winmsk':        
        fh.write('columns=[4]\n')
        names_string = 'names=["' + '","'.join(datasource_info_tags) + '"]'
        fh.write(names_string + '\n')
        fh.write('ops=["flag"]\n\n')
    else:        
        fields_string = 'fields = ["' + '","'.join(datasource_info_tags) + '"]'
        ops = ['concat'] * len(datasource_info_tags)
        ops_string = 'ops=["' + '","'.join(ops) + '"]'
        fh.write(fields_string + '\n')
        fh.write(ops_string + '\n\n')
    fh.close()
    return



if __name__ == "__main__":
    __main__()
