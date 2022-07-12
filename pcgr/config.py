#!/usr/bin/env python

from pcgr import pcgr_vars


def create_config(arg_dict):
    config_options = {
      'sample_id': arg_dict['sample_id'],
      'genome_assembly': arg_dict['genome_assembly'],
      'debug': arg_dict['debug'],
      'tumor_purity': 'NA',
      'tumor_ploidy': 'NA',
      'tumor_type': {
        'type': str(pcgr_vars.tsites[arg_dict['tsite']]),
      },
      'assay': arg_dict['assay'],
      'rna': {},
      'clinicaltrials': {
        'run': int(arg_dict['include_trials']),
      },
      'msi': {
        'run': int(arg_dict['estimate_msi_status'])
      }
    }

    if not arg_dict['tumor_purity'] is None:
        config_options['tumor_purity'] = float(arg_dict['tumor_purity'])
    if not arg_dict['tumor_ploidy'] is None:
        config_options['tumor_ploidy'] = float(arg_dict['tumor_ploidy'])

    config_options['other'] = {
      'vcfanno_n_proc': int(arg_dict['vcfanno_n_proc']),
      'vep_buffer_size': int(arg_dict['vep_buffer_size']),
      'vep_pick_order': str(arg_dict['vep_pick_order']),
      'vep_n_forks': int(arg_dict['vep_n_forks']),
      'visual_theme': str(arg_dict['report_theme']),
      'nonfloating_toc': int(arg_dict['report_nonfloating_toc']),
      'vep_no_intergenic': int(arg_dict['vep_no_intergenic']),
      'vep_regulatory': int(arg_dict['vep_regulatory']),
      'vep_gencode_all': int(arg_dict['vep_gencode_all']),
      'vcf2maf': int(arg_dict['vcf2maf']),
      'basic': int(arg_dict['basic']),
      'preserved_info_tags': str(arg_dict['preserved_info_tags']),
      'list_noncoding': int(arg_dict['show_noncoding']),
      'force_overwrite': int(arg_dict['force_overwrite']),
      'no_vcf_validate': int(arg_dict['no_vcf_validate']),
    }

    config_options['msigs'] = {
      'run': int(arg_dict['estimate_signatures']),
      'mutation_limit': int(arg_dict['min_mutations_signatures']),
      'all_reference_signatures': int(arg_dict['all_reference_signatures']),
      'include_artefact_signatures': int(arg_dict['include_artefact_signatures']),
      'prevalence_reference_signatures': int(arg_dict['prevalence_reference_signatures']),
    }

    config_options['tmb'] = {
      'run': int(arg_dict['estimate_tmb']),
      'target_size_mb': arg_dict['target_size_mb'],
      'algorithm': arg_dict['tmb_algorithm'],
    }

    config_options['cna'] = {
      'logR_homdel': float(arg_dict['logr_homdel']),
      'logR_gain': float(arg_dict['logr_gain']),
      'cna_overlap_pct': float(arg_dict['cna_overlap_pct']),
    }

    config_options['allelic_support'] = {
      'tumor_dp_min': int(arg_dict['tumor_dp_min']),
      'control_dp_min': int(arg_dict['control_dp_min']),
      'tumor_af_min': float(arg_dict['tumor_af_min']),
      'control_af_max': float(arg_dict['control_af_max']),
      'control_dp_tag': str(arg_dict['control_dp_tag']),
      'control_af_tag': str(arg_dict['control_af_tag']),
      'tumor_dp_tag': str(arg_dict['tumor_dp_tag']),
      'tumor_af_tag': str(arg_dict['tumor_af_tag']),
      'call_conf_tag': str(arg_dict['call_conf_tag']),
    }

    config_options['tumor_only'] = {
      'tumor_only': int(arg_dict['tumor_only']),
      'cell_line': int(arg_dict['cell_line']),
      'exclude_pon': int(arg_dict['exclude_pon']),
      'exclude_likely_hom_germline': int(arg_dict['exclude_likely_hom_germline']),
      'exclude_likely_het_germline': int(arg_dict['exclude_likely_het_germline']),
      'exclude_dbsnp_nonsomatic': int(arg_dict['exclude_dbsnp_nonsomatic']),
      'exclude_nonexonic': int(arg_dict['exclude_nonexonic']),
    }

    for pop in ['eur', 'afr', 'amr', 'eas', 'sas', 'global']:
        tag = 'maf_onekg_' + str(pop)
        if arg_dict[tag]:
            config_options['tumor_only'][tag] = float(arg_dict[tag])

    for pop in ['nfe', 'fin', 'amr', 'eas', 'sas', 'asj', 'oth', 'afr', 'global']:
        tag = 'maf_gnomad_' + str(pop)
        if arg_dict[tag]:
            config_options['tumor_only'][tag] = float(arg_dict[tag])

    return config_options
