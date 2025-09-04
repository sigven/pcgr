#!/usr/bin/env python

import csv
import re
import argparse
import os
import sys
#import pandas as pd
import gzip

from cyvcf2 import VCF
from pcgr.utils import error_message, check_subprocess, getlogger, random_id_generator, sort_bed, check_file_exists, remove_file
from pcgr.vcf import check_existing_vcf_info_tags, check_retained_vcf_info_tags
from pcgr.annoutils import read_infotag_file,read_vcfanno_tag_file


def __main__():

    parser = argparse.ArgumentParser(description='Verify input data for CPSR')
    parser.add_argument('refdata_assembly_dir',help='Assembly-specific directory containing the reference data files for PCGR/CPSR')
    parser.add_argument('input_vcf', help='VCF input file with query variants (SNVs/InDels)')
    parser.add_argument('validated_vcf',help="Name of VCF file with validated, decomposed query variants that overlap target genes (SNVs/InDels)")
    parser.add_argument('custom_target_tsv',help='Custom text/TSV file indicating user-defined target genes from panel 0 for screening and reporting')
    parser.add_argument('custom_target_bed',help='Name of BED file populated with regions associated with custom target genes defined by user')
    parser.add_argument('retained_info_tags',help='Comma-separated string of VCF INFO tags in query VCF to be retained for output')
    parser.add_argument('sample_id',help='CPSR sample_name')
    parser.add_argument('virtual_panel_id',type=str,help='virtual panel identifier(s)')
    parser.add_argument('diagnostic_grade_only', type=int, default=0, choices=[0,1], help="Green virtual panels only (Genomics England PanelApp)")
    parser.add_argument('gwas_findings', type=int, default=0, choices=[0,1], help="Include GWAS findings")
    parser.add_argument('secondary_findings', type=int, default=0, choices=[0,1], help="Include secondary findings")
    parser.add_argument('pgx_findings', type=int, default=0, choices=[0,1], help="Include classification of variants related to chemotherapy toxicity")
    parser.add_argument('--output_dir', dest='output_dir', help='Output directory')
    parser.add_argument("--debug", action="store_true", help="Print full commands to log")
    args = parser.parse_args()

    ret = validate_cpsr_input(args.refdata_assembly_dir,
                              args.input_vcf,
                              args.validated_vcf,
                              args.custom_target_tsv,
                              args.custom_target_bed,
                              args.retained_info_tags,
                              args.sample_id,
                              args.virtual_panel_id,
                              args.diagnostic_grade_only,
                              args.gwas_findings,
                              args.secondary_findings,
                              args.pgx_findings,
                              args.output_dir,
                              args.debug)
    if ret != 0:
       sys.exit(1)


def get_valid_custom_genelist(genelist_fname, genelist_bed_fname, refdata_assembly_dir,  
                              gwas_findings, secondary_findings, pgx_findings, logger, debug):
    """
    Function that checks whether the custom genelist contains valid entries from the complete exploratory track
    """
    
    random_id = random_id_generator(15)
    
    genelist_reader = csv.DictReader(open(genelist_fname,'r'), delimiter='\n', fieldnames=['ensembl_gene_id'])
    virtualpanel_track_bed = os.path.join(
        refdata_assembly_dir, "gene","bed","gene_virtual_panel", "0.bed.gz")
    virtualpanel_track_tsv = os.path.join(
        refdata_assembly_dir, "gene","tsv","gene_virtual_panel", "gene_virtual_panel.tsv.gz")
    genelist_bed_fname_unsorted = f'{genelist_bed_fname}.{random_id}.unsorted.bed'
    
    check_file_exists(virtualpanel_track_bed, logger)
    check_file_exists(virtualpanel_track_tsv, logger)

    customlist_identifiers = {}
    superpanel_track = []
    superpanel_identifiers_all = {}
    valid_custom_identifiers = []
    valid_custom_symbols = []

    for row in genelist_reader:
        if not re.match(r'^ENSG[0-9]{1,}$', str(row['ensembl_gene_id']).rstrip()):
            err_msg = "Custom list of genes from CPSR superpanel (panel 0) should be provided as Ensembl " + \
                "gene identifiers, '" + str(row['ensembl_gene_id']) + "' is not a valid identifier"
            return error_message(err_msg, logger)
        else:
            customlist_identifiers[str(row['ensembl_gene_id']).strip()] = 1

    with gzip.open(virtualpanel_track_tsv, mode='rt') as f:
        virtualpanel_reader = csv.DictReader(f, delimiter = '\t')
        for row in virtualpanel_reader:
            if row['id'] == '0':
                superpanel_track.append(dict(row))

    i = 0
    while i < len(superpanel_track):
        superpanel_identifiers_all[superpanel_track[i]['ensembl_gene_id']] = superpanel_track[i]['symbol']
        i = i + 1

    for g in customlist_identifiers.keys():
        if g in superpanel_identifiers_all.keys():
            valid_custom_identifiers.append(g)
            valid_custom_symbols.append(superpanel_identifiers_all[g])
        else:
            logger.warning("Ignoring custom-provided gene identifier (" + str(g) + ") NOT found in CPSR superpanel (panel 0)")
            logger.warning("Choose only Ensembl gene identifiers from panel 0 in this file : " + str(virtualpanel_track_tsv))
    all_valid_custom_geneset = ', '.join(sorted(valid_custom_symbols))

    logger.info('Detected n = ' + str(len(valid_custom_identifiers)) + ' valid targets in custom-provided gene list file (--custom_list)):')
    logger.info(all_valid_custom_geneset)

    if len(valid_custom_identifiers) == 0:
        logger.info('')
        logger.info("NO valid gene identifiers from panel 0 in custom-provided genelist - exiting")
        logger.info('')
        exit(1)

    ## Add custom set of genes to target BED
    logger.info('Creating BED file with custom target genes: ' + str(genelist_bed_fname))
    
    ## Create BED file with custom target genes, potentially filtering out GWAS tag SNPs, ACMG Secondary Findings, and CPIC PGx Oncology overlaps
    target_gene_pattern = '|'.join([f"{g}" for g in valid_custom_identifiers])
    filter_bed_command = get_bed_filtering_command(
        virtualpanel_track_bed, target_gene_pattern, gwas_findings, secondary_findings, pgx_findings)
    check_subprocess(logger, f'{filter_bed_command} > {genelist_bed_fname_unsorted}', debug)
    
    ## Sort regions in target BED
    sort_bed(genelist_bed_fname_unsorted, genelist_bed_fname, debug, logger)

    return 0

def get_bed_filtering_command(target_bed_gz, target_gene_pattern, gwas_findings, secondary_findings, pgx_findings):
    
    ## All targets - genes of interest + GWAS tag SNPS + ACMG Secondary Findings + CPIC PGx Oncology (chemo toxicity, DPYD)
    all_targets = target_gene_pattern + "|(\\|tag\\|)|ACMG_SF|CPIC_PGX_ONCOLOGY"
    
    ## Patterns to filter out GWAS tag SNPs, ACMG Secondary Findings, and CPIC PGx Oncology
    ignore_acmg_sf_targets = "awk 'BEGIN{FS=\"\\t\"}{if($4 !~ /ACMG_SF/ || ($4 ~ /ACMG_SF/ && $4 ~ /" + str(target_gene_pattern) + "/))print;}'"
    ignore_gwas_targets = "egrep -v '(\\|tag\\|)'"
    ignore_pgx_targets = "egrep -v 'CPIC_PGX_ONCOLOGY'"
    
    filter_bed_command = f"bgzip -dc {target_bed_gz} | egrep '{all_targets}'"
    if gwas_findings == 0:
        filter_bed_command += f' | {ignore_gwas_targets}'
    if secondary_findings == 0:
        filter_bed_command += f' | {ignore_acmg_sf_targets}'
    if pgx_findings == 0:
        filter_bed_command += f' | {ignore_pgx_targets}'
    return filter_bed_command

def simplify_vcf(input_vcf, validated_vcf, vcf, custom_bed, refdata_assembly_dir, virtual_panel_id, 
                 sample_id, diagnostic_grade_only, gwas_findings, secondary_findings, pgx_findings, output_dir, logger, debug):

    """
    Function that performs four separate checks/filters on the validated input VCF:
    1. Remove/Strip off any genotype data (not needed for annotation)
    2. If VCF have variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'), 
        these are decomposed into variants with a single alternative allele
    3. Filters against predisposition loci (virtual panel id or custom target) - includes secondary finding targets/GWAS-loci if set by user
    4. Final VCF file is sorted and indexed (bgzip + tabix)
    """

    random_str = random_id_generator(10) 

    temp_files = {}
    temp_files['vcf_1'] = \
        os.path.join(output_dir, f'{sample_id}.cpsr_validate.bcftools.{random_str}_1.vcf')
    temp_files['vcf_2'] = \
        os.path.join(output_dir, f'{sample_id}.cpsr_validate.bcftools.{random_str}_2.vcf.gz')
    temp_files['vcf_3'] = \
        os.path.join(output_dir, f'{sample_id}.cpsr_validate.bftools.{random_str}_3.vcf')
    temp_files['vcf_4'] = \
        os.path.join(output_dir, f'{sample_id}.cpsr_validate.decomp.{random_str}_4.vcf')
    bcftools_simplify_log = \
        os.path.join(output_dir, f'{sample_id}.cpsr_validate.bcftools.{random_str}.log')
    vt_decompose_log = \
        os.path.join(output_dir, f'{sample_id}.cpsr_validate.vt_decompose.{random_str}.log')
    virtual_panels_tmp_bed = \
        os.path.join(output_dir, f'{sample_id}.cpsr_virtual_panels_all.{random_str}.tmp.bed')
    virtual_panels_bed = \
        os.path.join(output_dir, f'{sample_id}.cpsr_virtual_panels_all.{random_str}.bed')

    multiallelic_list = list()
    for rec in vcf:
        POS = rec.start + 1
        alt = ",".join(str(n) for n in rec.ALT)
        if len(rec.ALT) > 1:
            variant_id = f"{rec.CHROM}:{POS}_{rec.REF}->{alt}"
            multiallelic_list.append(variant_id)

    # bgzip + tabix required for sorting
    cmd_vcf1 = f'bcftools view {input_vcf} | bgzip -cf > {temp_files["vcf_2"]} && tabix -p vcf {temp_files["vcf_2"]} && ' + \
        f'bcftools sort --temp-dir {output_dir} -Oz {temp_files["vcf_2"]} > {temp_files["vcf_3"]} 2> {bcftools_simplify_log}' + \
        f' && tabix -p vcf {temp_files["vcf_3"]}'
    logger.info('Extracting variants on autosomal/sex/mito chromosomes only (1-22,X,Y) with bcftools')
    # Keep only autosomal/sex chromosomes, sub chr prefix
    # Note: M/MT variants are skipped - requires additional cache/handling from VEP, 
    # see e.g. https://github.com/Ensembl/ensembl-vep/issues/464
    chrom_to_keep = [str(x) for x in [*range(1,23), 'X', 'Y',]]
    chrom_to_keep = ','.join([*['chr' + chrom for chrom in chrom_to_keep], *[chrom for chrom in chrom_to_keep]])
    cmd_vcf2 = f'bcftools view --regions {chrom_to_keep} {temp_files["vcf_3"]} | sed \'s/^chr//\' > {temp_files["vcf_1"]}'

    check_subprocess(logger, cmd_vcf1, debug)
    check_subprocess(logger, cmd_vcf2, debug)

    if multiallelic_list:
        logger.warning(f"There were {len(multiallelic_list)} multiallelic sites detected. Showing (up to) the first 100:")
        print('----')
        print(', '.join(multiallelic_list[:100]))
        print('----')
        logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
        command_decompose = f'vt decompose -s {temp_files["vcf_1"]} > {temp_files["vcf_4"]} 2> {vt_decompose_log}'
        check_subprocess(logger, command_decompose, debug)
    else:
        logger.info('All sites seem to be decomposed - skipping decomposition of multiallelic sites')
        command_copy = f'cp {temp_files["vcf_1"]} {temp_files["vcf_4"]}'
        check_subprocess(logger, command_copy, debug)


    if not custom_bed == 'None':
        logger.info('Limiting variant set to user-defined screening loci (custom list from panel 0)')
        if check_file_exists(custom_bed):
            target_variants_intersect_cmd = \
                "bedtools intersect -wa -u -header -a " + str(temp_files['vcf_4']) + \
                " -b " + str(custom_bed) + " > " + str(validated_vcf)
            check_subprocess(logger, target_variants_intersect_cmd, debug)
    else:
        logger.info("Limiting variant set to cancer predisposition loci (virtual panel id(s): '" + str(virtual_panel_id) + "')")      

        ## Concatenate all panel BEDs to one big virtual panel BED, sort and make unique
        panel_ids = str(virtual_panel_id).split(',')
        for pid in set(panel_ids):
            ge_panel_identifier = "GE_PANEL_" + str(pid)
            target_bed_gz = os.path.join(
                refdata_assembly_dir, 'gene','bed','gene_virtual_panel', str(pid) + ".bed.gz")
            if diagnostic_grade_only == 1 and virtual_panel_id != "0":
                logger.info('Considering diagnostic-grade only genes in panel ' + str(pid) + ' - (GREEN status in Genomics England PanelApp)')
                target_bed_gz = os.path.join(
                    refdata_assembly_dir, 'gene','bed','gene_virtual_panel', str(pid) + ".GREEN.bed.gz")            
            check_file_exists(target_bed_gz, logger)
            
            ## Append the virtual panel BED file with specific panel target genes, potentially filtering out GWAS hits, 
            ## ACMG Secondary Findings, and CPIC PGx Oncology overlapp (i.e. based on user parameters)
            target_gene_pattern = str(ge_panel_identifier) + ":"
            if pid == '0':
                target_gene_pattern = "ENSG[0-9]{1,}"
            filter_bed_command = get_bed_filtering_command(
                target_bed_gz, target_gene_pattern, gwas_findings, secondary_findings, pgx_findings)
            check_subprocess(logger, f'{filter_bed_command} >> {virtual_panels_tmp_bed}', debug)

        ## sort the collection of virtual panels
        sort_bed(virtual_panels_tmp_bed, virtual_panels_bed, debug, logger)

        if check_file_exists(virtual_panels_bed):
            target_variants_intersect_cmd = f'bedtools intersect -wa -u -header -a {temp_files["vcf_4"]} -b ' + \
                f'{virtual_panels_bed} > {validated_vcf}'
            check_subprocess(logger, target_variants_intersect_cmd, debug)


    check_subprocess(logger, f'bgzip -cf {validated_vcf} > {validated_vcf}.gz', debug)
    check_subprocess(logger, f'tabix -p vcf {validated_vcf}.gz', debug)
    if not debug:
        for fn in [virtual_panels_bed, 
                   temp_files["vcf_1"],
                   temp_files["vcf_2"],
                   temp_files["vcf_3"],
                   temp_files["vcf_4"],
                   bcftools_simplify_log, 
                   vt_decompose_log]:
            #print(f"Deleting {fn}")
            remove_file(fn)
        
        remove_file(temp_files["vcf_2"] + str('.tbi'))
        remove_file(temp_files["vcf_3"] + str('.tbi'))

    if check_file_exists(f'{validated_vcf}.gz'):
        vcf = VCF(validated_vcf + '.gz')
        i = 0
        for rec in vcf:
            i = i + 1
        if len(vcf.seqnames) == 0 or i == 0:
            logger.info('')
            logger.info("Query VCF contains NO variants within the selected cancer predisposition geneset (or "\
                "GWAS loci/secondary findings) - quitting workflow")
            logger.info('')
            exit(1)

def validate_cpsr_input(refdata_assembly_dir, 
                        input_vcf, 
                        validated_vcf,
                        custom_list_fname, 
                        custom_list_bed_fname,
                        retained_info_tags, 
                        sample_id, 
                        virtual_panel_id, 
                        diagnostic_grade_only, 
                        gwas_findings, 
                        secondary_findings, 
                        pgx_findings,
                        output_dir, 
                        debug):
    """
    Function that reads the input files to CPSR (VCF file + custom gene list) and performs the following checks:
    0. If custom gene list (panel) is provided, checks the validity of this list
    1. Check that no INFO annotation tags in the query VCF coincides with those generated by CPSR
    2. Check that custom VCF INFO tags set by user as retained for output is found in query VCF
    3. Check that if VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
    4. The resulting VCF file is sorted and indexed (bgzip + tabix)
    """
    logger = getlogger('cpsr-validate-input-arguments')

    custom_target_fname = {}
    custom_target_fname['tsv'] = custom_list_fname
    custom_target_fname['bed'] = 'None'
    
    if not custom_target_fname['tsv'] == 'None':
        logger.info('Establishing BED track with custom list of genes from panel 0')
        custom_target_fname['bed'] = custom_list_bed_fname
        get_valid_custom_genelist(custom_target_fname['tsv'], 
                                  custom_target_fname['bed'], 
                                  refdata_assembly_dir, 
                                  gwas_findings, 
                                  secondary_findings, 
                                  pgx_findings,
                                  logger, 
                                  debug)
    
    if not input_vcf == 'None':
        
        vcf_object = VCF(input_vcf)
        
        ## Check that VCF does not already contain INFO tags that will be appended through CPSR annotation
        ## - First add info tags generated by VEP, and those generated by CPSR
        populated_infotags_other_fname = os.path.join(refdata_assembly_dir, 'vcf_infotags_other.tsv')
        populated_infotags_vep_fname = os.path.join(refdata_assembly_dir, 'vcf_infotags_vep.tsv')
        tags_cpsr = read_infotag_file(populated_infotags_other_fname, scope = "cpsr")
        tags_vep = read_infotag_file(populated_infotags_vep_fname, scope = "vep")
        tags_cpsr.update(tags_vep)
        
        
        ## - Next, add INFO tags generated through vcfanno annotation
        track_file_info = {}
        track_file_info['tags_fname'] = {}
        for variant_track in ['clinvar','tcga','gwas','dbmts','dbnsfp','gnomad_non_cancer']:
            track_file_info['tags_fname'][variant_track] = os.path.join(
                refdata_assembly_dir,'variant','vcf', variant_track, f'{variant_track}.vcfanno.vcf_info_tags.txt')

        for bed_track in ['simplerepeat','winmsk','rmsk','gerp']:
            track_file_info['tags_fname'][bed_track] = os.path.join(
                refdata_assembly_dir,'misc','bed', bed_track, f'{bed_track}.vcfanno.vcf_info_tags.txt')

        track_file_info['tags_fname']['gene_transcript_xref'] = os.path.join(
            refdata_assembly_dir,'gene','bed', 'gene_transcript_xref', 'gene_transcript_xref.vcfanno.vcf_info_tags.txt')
    
        for track in track_file_info['tags_fname']:
            infotags_vcfanno = read_vcfanno_tag_file(track_file_info['tags_fname'][track], logger)
            tags_cpsr.update(infotags_vcfanno)
        
        ## Check that no INFO annotation tags in the query VCF coincides with those generated by CPSR
        tag_check = check_existing_vcf_info_tags(vcf_object, tags_cpsr, logger)
        if tag_check == -1:
            return -1        

        ## Check that retained VCF INFO tags requested by user are present in VCF
        if retained_info_tags != "None":
            custom_check = check_retained_vcf_info_tags(vcf_object, retained_info_tags, logger)
            if custom_check == -1:
                return -1

        samples = vcf_object.samples
        if len(samples) > 1:
            err_msg = "Query VCF contains more than one sample column (" + ', '.join(samples) + ") - " + \
                "CPSR expects a germline VCF with a single sample column - exiting"
            return error_message(err_msg, logger)
        
        simplify_vcf(input_vcf, 
                     validated_vcf,
                     vcf_object, 
                     custom_target_fname['bed'], 
                     refdata_assembly_dir, 
                     virtual_panel_id, 
                     sample_id, 
                     diagnostic_grade_only, 
                     gwas_findings, 
                     secondary_findings, 
                     pgx_findings,
                     output_dir, 
                     logger, 
                     debug)

    return 0

if __name__=="__main__":
    __main__()
