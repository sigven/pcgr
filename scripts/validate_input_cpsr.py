#!/usr/bin/env python

import csv
import re
import argparse
import os
import sys
import gzip

from cyvcf2 import VCF
from pcgr.utils import getlogger, error_message, check_subprocess, random_id_generator, sort_bed, check_file_exists, remove_file
from pcgr.vcf import check_existing_vcf_info_tags, check_retained_vcf_info_tags
from pcgr.validate import bgzip_sort_filter_vcf, detect_multiallelic_sites, collect_workflow_info_tags


def __main__():

    parser = argparse.ArgumentParser(description='Verify input data for CPSR')
    parser.add_argument('refdata_assembly_dir',help='Assembly-specific directory containing the reference data files for PCGR/CPSR')
    parser.add_argument('input_vcf', help='VCF input file with query variants (SNVs/InDels)')
    parser.add_argument('validated_vcf', help="Name of VCF file with validated, decomposed query variants that overlap target genes (SNVs/InDels)")
    parser.add_argument('input_sv_vcf', help='VCF input file with germline structural variants (use "None" if absent)')
    parser.add_argument('validated_sv_vcf', help='Name of VCF file for validated germline SV variants (use "None" if absent)')
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
                              args.input_sv_vcf,
                              args.validated_sv_vcf,
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

def filter_hom_ref_variants(input_vcf, output_vcf, logger, debug):
    """
    Remove homozygous-reference (0/0) records from a VCF using bcftools.
    Returns the number of records removed.
    """
    import subprocess

    def _count_records(vcf_path):
        result = subprocess.run(
            f'bcftools view -H {vcf_path} 2>/dev/null | wc -l',
            shell=True, capture_output=True, text=True)
        return int(result.stdout.strip())

    n_total = _count_records(input_vcf)
    # GT="ref" matches only hom-ref (0/0); 2>/dev/null suppresses ##contig header warnings
    check_subprocess(logger, f'bcftools view -e \'GT="ref"\' {input_vcf} 2>/dev/null > {output_vcf}', debug)
    n_after = _count_records(output_vcf)
    n_removed = n_total - n_after
    if n_removed > 0:
        logger.warning(
            f"Filtered out {n_removed} of {n_total} variant(s) with a homozygous-reference genotype (0/0) "
        )
    return n_removed


def simplify_vcf(input_vcf, validated_vcf, vcf_obj, custom_bed, refdata_assembly_dir, virtual_panel_id,
                 sample_id, diagnostic_grade_only, gwas_findings, secondary_findings, pgx_findings, output_dir, logger, debug):

    """
    Function that performs four separate checks/filters on the validated input VCF:
    1. Remove/Strip off any genotype data (not needed for annotation)
    2. If VCF have variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'),
        these are decomposed into variants with a single alternative allele
    3. Filters against predisposition loci (virtual panel id or custom target) - includes secondary finding targets/GWAS-loci if set by user
    4. Homozygous-reference (0/0) variants are removed
    5. Final VCF file is sorted and indexed (bgzip + tabix)
    """

    random_str = random_id_generator(10)

    vcf_filtered         = os.path.join(output_dir, f'{sample_id}.cpsr_validate.filtered.{random_str}.vcf')
    vcf_decomposed       = os.path.join(output_dir, f'{sample_id}.cpsr_validate.decomp.{random_str}.vcf')
    vcf_hom_ref_filtered = os.path.join(output_dir, f'{sample_id}.cpsr_validate.hom_ref_filtered.{random_str}.vcf')
    vt_decompose_log     = os.path.join(output_dir, f'{sample_id}.cpsr_validate.vt_decompose.{random_str}.log')
    virtual_panels_tmp_bed = os.path.join(output_dir, f'{sample_id}.cpsr_virtual_panels_all.{random_str}.tmp.bed')
    virtual_panels_bed   = os.path.join(output_dir, f'{sample_id}.cpsr_virtual_panels_all.{random_str}.bed')

    multiallelic_list = detect_multiallelic_sites(vcf_obj)

    # Sort, chromosome-filter, keep genotype columns — temp files managed inside
    bgzip_sort_filter_vcf(input_vcf, vcf_filtered, output_dir, sample_id, 'cpsr', logger, debug)

    if multiallelic_list:
        logger.warning(f"There were {len(multiallelic_list)} multiallelic sites detected. Showing (up to) the first 100:")
        print('----')
        print(', '.join(multiallelic_list[:100]))
        print('----')
        logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
        command_decompose = f'vt decompose -s {vcf_filtered} > {vcf_decomposed} 2> {vt_decompose_log}'
        check_subprocess(logger, command_decompose, debug)
    else:
        logger.info('All sites seem to be decomposed - skipping decomposition of multiallelic sites')
        check_subprocess(logger, f'cp {vcf_filtered} {vcf_decomposed}', debug)

    if not custom_bed == 'None':
        logger.info('Limiting variant set to user-defined screening loci (custom list from panel 0)')
        if check_file_exists(custom_bed):
            target_variants_intersect_cmd = \
                f'bedtools intersect -wa -u -header -a {vcf_decomposed} -b {custom_bed} > {validated_vcf}'
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
            ## ACMG Secondary Findings, and CPIC PGx Oncology overlap (i.e. based on user parameters)
            target_gene_pattern = str(ge_panel_identifier) + ":"
            if pid == '0':
                target_gene_pattern = "ENSG[0-9]{1,}"
            filter_bed_command = get_bed_filtering_command(
                target_bed_gz, target_gene_pattern, gwas_findings, secondary_findings, pgx_findings)
            check_subprocess(logger, f'{filter_bed_command} >> {virtual_panels_tmp_bed}', debug)

        ## sort the collection of virtual panels
        sort_bed(virtual_panels_tmp_bed, virtual_panels_bed, debug, logger)

        if check_file_exists(virtual_panels_bed):
            target_variants_intersect_cmd = (
                f'bedtools intersect -wa -u -header -a {vcf_decomposed} '
                f'-b {virtual_panels_bed} > {validated_vcf}')
            check_subprocess(logger, target_variants_intersect_cmd, debug)

    ## Filter out homozygous-reference (0/0) records — not informative for germline classification
    if len(vcf_obj.samples) > 0:
        filter_hom_ref_variants(validated_vcf, vcf_hom_ref_filtered, logger, debug)
        check_subprocess(logger, f'mv {vcf_hom_ref_filtered} {validated_vcf}', debug)

    check_subprocess(logger, f'bgzip -cf {validated_vcf} > {validated_vcf}.gz', debug)
    check_subprocess(logger, f'tabix -p vcf {validated_vcf}.gz', debug)
    if not debug:
        for fn in [virtual_panels_bed, vcf_filtered, vcf_decomposed, vt_decompose_log, vcf_hom_ref_filtered]:
            remove_file(fn)

    if check_file_exists(f'{validated_vcf}.gz'):
        vcf = VCF(validated_vcf + '.gz')
        i = 0
        n_hom_ref = 0
        has_gt = len(vcf.samples) > 0
        for rec in vcf:
            i += 1
            if has_gt:
                gt = rec.genotypes
                if gt and len(gt) > 0:
                    alleles = gt[0][:2]
                    if alleles == [0, 0]:
                        n_hom_ref += 1
        if len(vcf.seqnames) == 0 or i == 0:
            logger.info('')
            logger.info("Query VCF contains NO variants within the selected cancer predisposition geneset (or "\
                "GWAS loci/secondary findings) - quitting workflow")
            logger.info('')
            exit(1)
        if has_gt and n_hom_ref == i and i > 0:
            logger.error('')
            logger.error(f"All {i} variant(s) in the query VCF have a homozygous-reference genotype (i.e. '0/0'). "
                         "CPSR expects germline variants with heterozygous/homozygous genotypes ('0/1' or '1/1'). "                          
                         "Please verify your input VCF - exiting.")
            logger.error('')
            exit(1)

def validate_cpsr_sv_input(input_sv_vcf, validated_sv_vcf, sample_id, output_dir, logger, debug):
    """
    Validate and clean a germline structural variant VCF for CPSR.

    Checks performed:
    1. SVTYPE INFO tag is present in the VCF header and all records carry a recognised value
       (DEL, DUP, INS, INV, BND)
    2. END or SVLEN INFO tag is present for non-BND record types
    3. A sample (genotype) column is present — expected for germline SV callers
    4. No symbolic-allele-free records (i.e. no plain SNV/indel records slipping through)
    5. Chromosome filtering (1-22, X, Y), bgzip compression, and tabix indexing
       — multiallelic decomposition via vt is deliberately skipped (not valid for SVs)
    """
    valid_svtypes = {'DEL', 'DUP', 'INS', 'INV', 'BND'}
    random_str = random_id_generator(10)

    vcf_object = VCF(input_sv_vcf)

    ## Check 1: SVTYPE INFO tag in header
    header_info_ids = [h['ID'] for h in vcf_object.header_iter() if h.type == 'INFO']
    if 'SVTYPE' not in header_info_ids:
        err_msg = "SV VCF is missing the SVTYPE INFO tag in the header — required for structural variant input"
        return error_message(err_msg, logger)

    ## Check 2: Sample/genotype column present (germline SV callers always emit one)
    if len(vcf_object.samples) == 0:
        err_msg = "SV VCF contains no sample column — expecting a single-sample VCF"
        return error_message(err_msg, logger)
    if len(vcf_object.samples) > 1:
        err_msg = (f"SV VCF contains more than one sample column ({', '.join(vcf_object.samples)}) — "
                   f"Expecting a single-sample VCF")
        return error_message(err_msg, logger)

    ## Check 3 & 4: Per-record SVTYPE validity, END/SVLEN presence, no plain SNV/indel ALTs
    n_records = 0
    invalid_svtype = []
    missing_span = []
    plain_allele = []
    for rec in vcf_object:
        n_records += 1
        svtype = rec.INFO.get('SVTYPE')
        if svtype is None or svtype not in valid_svtypes:
            invalid_svtype.append(f"{rec.CHROM}:{rec.POS} SVTYPE={svtype}")
        # Flag records that have neither a recognised SVTYPE nor symbolic/BND ALT notation —
        # these are likely accidental SNV/indel records. Records with a valid SVTYPE are
        # allowed to carry a literal sequence ALT (e.g. Manta writes the inserted sequence
        # explicitly for SVTYPE=INS rather than using <INS>).
        if svtype not in valid_svtypes:
            alt_alleles = [str(a) for a in rec.ALT]
            for alt in alt_alleles:
                if not (alt.startswith('<') or '[' in alt or ']' in alt):
                    plain_allele.append(f"{rec.CHROM}:{rec.POS} ALT={alt[:60]}{'...' if len(alt) > 60 else ''}")
        if svtype != 'BND':
            if rec.INFO.get('END') is None and rec.INFO.get('SVLEN') is None:
                missing_span.append(f"{rec.CHROM}:{rec.POS} SVTYPE={svtype}")

    if invalid_svtype:
        logger.warning(
            f"{len(invalid_svtype)} SV record(s) carry an unrecognised or missing SVTYPE "
            f"(expected one of: {', '.join(sorted(valid_svtypes))}). "
            f"First offenders: {'; '.join(invalid_svtype[:5])}"
        )
    if missing_span:
        logger.warning(
            f"{len(missing_span)} non-BND SV record(s) lack both END and SVLEN INFO tags. "
            f"First offenders: {'; '.join(missing_span[:5])}"
        )
    if plain_allele:
        err_msg = (f"SV VCF contains {len(plain_allele)} record(s) with no recognised SVTYPE and "
                   f"no symbolic/BND ALT allele — these appear to be SNV/indel records. "
                   f"First offenders: {'; '.join(plain_allele[:5])}")
        return error_message(err_msg, logger)

    if n_records == 0:
        err_msg = "SV VCF contains no variant records — exiting"
        return error_message(err_msg, logger)

    ## Filter to autosomal/sex chromosomes, strip chr prefix, bgzip + tabix
    ## Note: vt decompose is intentionally omitted — not valid for symbolic SV alleles
    sv_filtered = os.path.join(output_dir, f'{sample_id}.cpsr_sv_validate.filtered.{random_str}.vcf')
    bgzip_sort_filter_vcf(input_sv_vcf, sv_filtered, output_dir, sample_id, 'cpsr', logger, debug)

    check_subprocess(logger, f'bgzip -cf {sv_filtered} > {validated_sv_vcf}.gz', debug)
    check_subprocess(logger, f'tabix -p vcf {validated_sv_vcf}.gz', debug)

    if not debug:
        remove_file(sv_filtered)

    logger.info(f'SV VCF validation complete — {n_records} record(s) processed')
    return 0


def validate_cpsr_input(refdata_assembly_dir,
                        input_vcf,
                        validated_vcf,
                        input_sv_vcf,
                        validated_sv_vcf,
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
        tags_cpsr = collect_workflow_info_tags(refdata_assembly_dir, 'cpsr', logger)

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

    if not input_sv_vcf == 'None':
        sv_ret = validate_cpsr_sv_input(
            input_sv_vcf, validated_sv_vcf, sample_id, output_dir, logger, debug)
        if sv_ret != 0:
            return -1

    return 0

if __name__=="__main__":
    __main__()
