#!/usr/bin/env python

from pcgr import pcgr_vars


def create_config(arg_dict):
    config_options = {}
    config_options['tumor_purity'] = 'NA'
    config_options['tumor_ploidy'] = 'NA'
    if not arg_dict['tumor_purity'] is None:
        config_options['tumor_purity'] = float(arg_dict['tumor_purity'])
    if not arg_dict['tumor_ploidy'] is None:
        config_options['tumor_ploidy'] = float(arg_dict['tumor_ploidy'])

    config_options['assay'] = arg_dict['assay']

    config_options['other'] = {}
    config_options['other']['vcfanno_n_proc'] = int(arg_dict['vcfanno_n_proc'])
    config_options['other']['vep_buffer_size'] = int(arg_dict['vep_buffer_size'])
    config_options['other']['vep_pick_order'] = str(arg_dict['vep_pick_order'])
    config_options['other']['vep_n_forks'] = int(arg_dict['vep_n_forks'])
    config_options['other']['visual_theme'] = str(arg_dict['report_theme'])
    config_options['other']['nonfloating_toc'] = int(arg_dict['report_nonfloating_toc'])
    config_options['other']['vep_no_intergenic'] = int(arg_dict['vep_no_intergenic'])
    config_options['other']['vep_regulatory'] = int(arg_dict['vep_regulatory'])
    config_options['other']['vcf2maf'] = int(arg_dict['vcf2maf'])
    config_options['other']['basic'] = int(arg_dict['basic'])
    config_options['other']['preserved_info_tags'] = str(arg_dict['preserved_info_tags'])
    config_options['other']['list_noncoding'] = int(arg_dict['show_noncoding'])
    config_options['other']['force_overwrite'] = int(arg_dict['force_overwrite'])
    config_options['other']['no_vcf_validate'] = int(arg_dict['no_vcf_validate'])

    config_options['rna'] = {}
    config_options['clinicaltrials'] = {}
    config_options['clinicaltrials']['run'] = int(arg_dict['include_trials'])

    config_options['msi'] = {}
    config_options['msi']['run'] = int(arg_dict['estimate_msi_status'])
    config_options['msigs'] = {}
    config_options['msigs']['run'] = int(arg_dict['estimate_signatures'])
    config_options['msigs']['mutation_limit'] = int(arg_dict['min_mutations_signatures'])
    config_options['msigs']['all_reference_signatures'] = int(arg_dict['all_reference_signatures'])
    config_options['msigs']['include_artefact_signatures'] = int(arg_dict['include_artefact_signatures'])

    config_options['tmb'] = {}
    config_options['tmb']['run'] = int(arg_dict['estimate_tmb'])
    config_options['tmb']['target_size_mb'] = arg_dict['target_size_mb']
    config_options['tmb']['algorithm'] = arg_dict['tmb_algorithm']

    config_options['cna'] = {}
    config_options['cna']['logR_homdel'] = float(arg_dict['logr_homdel'])
    config_options['cna']['logR_gain'] = float(arg_dict['logr_gain'])
    config_options['cna']['cna_overlap_pct'] = float(arg_dict['cna_overlap_pct'])

    config_options['allelic_support'] = {}
    config_options['allelic_support']['tumor_dp_min'] = int(arg_dict['tumor_dp_min'])
    config_options['allelic_support']['control_dp_min'] = int(arg_dict['control_dp_min'])
    config_options['allelic_support']['tumor_af_min'] = float(arg_dict['tumor_af_min'])
    config_options['allelic_support']['control_af_max'] = float(arg_dict['control_af_max'])
    config_options['allelic_support']['control_dp_tag'] = str(arg_dict['control_dp_tag'])
    config_options['allelic_support']['control_af_tag'] = str(arg_dict['control_af_tag'])
    config_options['allelic_support']['tumor_dp_tag'] = str(arg_dict['tumor_dp_tag'])
    config_options['allelic_support']['tumor_af_tag'] = str(arg_dict['tumor_af_tag'])
    config_options['allelic_support']['call_conf_tag'] = str(arg_dict['call_conf_tag'])

    config_options['tumor_only'] = {}
    config_options['tumor_only']['tumor_only'] = int(arg_dict['tumor_only'])
    config_options['tumor_only']['cell_line'] = int(arg_dict['cell_line'])
    config_options['tumor_only']['exclude_pon'] = int(arg_dict['exclude_pon'])
    config_options['tumor_only']['exclude_likely_hom_germline'] = int(arg_dict['exclude_likely_hom_germline'])
    config_options['tumor_only']['exclude_likely_het_germline'] = int(arg_dict['exclude_likely_het_germline'])
    config_options['tumor_only']['exclude_dbsnp_nonsomatic'] = int(arg_dict['exclude_dbsnp_nonsomatic'])
    config_options['tumor_only']['exclude_nonexonic'] = int(arg_dict['exclude_nonexonic'])

    for pop in ['eur', 'afr', 'amr', 'eas', 'sas', 'global']:
        tag = 'maf_onekg_' + str(pop)
        if arg_dict[tag]:
            config_options['tumor_only'][tag] = float(arg_dict[tag])

    for pop in ['nfe', 'fin', 'amr', 'eas', 'sas', 'asj', 'oth', 'afr', 'global']:
        tag = 'maf_gnomad_' + str(pop)
        if arg_dict[tag]:
            config_options['tumor_only'][tag] = float(arg_dict[tag])

    config_options['tumor_type'] = {}
    config_options['tumor_type']['type'] = str(pcgr_vars.tsites[arg_dict['tsite']])

    return config_options
