#!/usr/bin/env python
"""
Import of homologous recombination deficiency (HRD) results from external HRD
assays/tools.

An external HRD result is read from a single tab-separated file with a header
and one data row (one sample), see parse_hrd_generic(), and normalised to a
common representation (see pcgr_vars.HRD_EXTERNAL_COLUMNS).

'HRD' covers measurements of different kinds - genomic scar scores (e.g. the
genomic instability score, HRD-LOH + LST + TAI, or genome-wide LOH), classifiers
based on mutational signatures (e.g. CHORD), and models that combine both (e.g.
HRDetect). Scores of different kinds are on different scales, and are not
comparable. The kind of score is therefore not inferred by PCGR, but stated by
the user in the file (HRD_SCORE_TYPE, see pcgr_vars.HRD_SCORE_TYPES).

The HRD status (HR-deficient/HR-proficient) is either reported in the file
(HRD_STATUS), or derived by PCGR from the score with a user-defined threshold
(--input_hrd_threshold). Thresholds are specific to the assay/tool (and often to
the tumor type), hence there is no default threshold.

File names are not used for any decision - the content (header) of the file is.
"""

import logging
import os
from typing import Optional

from pcgr import pcgr_vars
from pcgr.external_input import (
    ExternalInputError, is_missing, match_synonym, optional_text, read_single_row_tsv,
    to_count, to_float, write_external_result)
from pcgr.utils import error_message, warn_message

## accepted spellings of the HRD status (compared after upper-casing and mapping '_'/' ' to '-')
_STATUS_SYNONYMS = {
    pcgr_vars.HRD_STATUS_DEFICIENT: {'HR-DEFICIENT', 'HRDEFICIENT', 'DEFICIENT', 'HRD', 'HR-D',
                                     'HRD-POSITIVE', 'HRD+', 'POSITIVE'},
    pcgr_vars.HRD_STATUS_PROFICIENT: {'HR-PROFICIENT', 'HRPROFICIENT', 'PROFICIENT', 'HRP', 'HR-P',
                                      'HRD-NEGATIVE', 'HRD-', 'NEGATIVE'},
    pcgr_vars.HRD_STATUS_INDETERMINATE: {'INDETERMINATE', 'UNDETERMINED', 'EQUIVOCAL', 'INCONCLUSIVE',
                                         'UNKNOWN', 'BORDERLINE'},
}

_SCORE_COMPONENTS = ['HRD_LOH', 'HRD_LST', 'HRD_TAI']


class HrdInputError(ExternalInputError):
    """Raised when an external HRD result file cannot be read or has invalid content."""


def _normalise_status(value) -> str:
    status = match_synonym(value, _STATUS_SYNONYMS)
    if status is None:
        raise HrdInputError(
            f"Unrecognised HRD status ('{value}') - must be one of: "
            f"{pcgr_vars.HRD_STATUS_DEFICIENT}, {pcgr_vars.HRD_STATUS_PROFICIENT}, "
            f"{pcgr_vars.HRD_STATUS_INDETERMINATE} "
            "(or a common synonym, e.g. 'HRD'/'HRP', 'Positive'/'Negative')")
    return status


def _normalise_score_type(value) -> str:
    key = str(value).strip().lower().replace('-', '_').replace(' ', '_')
    if key not in pcgr_vars.HRD_SCORE_TYPES:
        raise HrdInputError(
            f"Unrecognised HRD score type ('{value}') in 'HRD_SCORE_TYPE' - must be one of: "
            f"{', '.join(pcgr_vars.HRD_SCORE_TYPES)}")
    return key


def parse_hrd_generic(fname: str) -> dict:
    """
    Read an external HRD result: a tab-separated file with a header and exactly
    one data row (one sample). Columns:

      HRD_CALLER        required - name of the HRD assay/tool
      HRD_SCORE_TYPE    required - kind of score, one of: genomic_scar,
                        mutational_signature, composite, other (see pcgr_vars.HRD_SCORE_TYPES)
      HRD_STATUS        optional* - HR-deficient, HR-proficient or Indeterminate
                        (common synonyms such as 'HRD'/'HRP' or 'Positive'/'Negative'
                        are accepted, case-insensitive)
      HRD_SCORE         optional* - numeric score as reported by the tool (higher
                        values indicating HR deficiency)
      CALLER_VERSION    optional
      HRD_SCORE_DESCRIPTION
                        optional - description of the score, e.g. 'GIS (HRD-LOH + LST + TAI)'
      HRD_SCORE_THRESHOLD_USED
                        optional - threshold on the score applied by the tool for HRD_STATUS
      HRD_LOH, HRD_LST, HRD_TAI
                        optional - genomic scar components (non-negative integers)

      * at least one of HRD_STATUS and HRD_SCORE must be given

    Additional columns are ignored. Missing optional values may be left empty
    or set to '.' / 'NA'. The status is not resolved here (see resolve_hrd_status()).

    Raises HrdInputError if the file is invalid.
    """
    row = read_single_row_tsv(fname, 'HRD', ['HRD_CALLER', 'HRD_SCORE_TYPE'], HrdInputError)
    if 'HRD_STATUS' not in row and 'HRD_SCORE' not in row:
        raise HrdInputError(f"HRD file ({os.path.basename(fname)}) must have an 'HRD_STATUS' and/or an 'HRD_SCORE' column")

    if is_missing(row['HRD_CALLER']):
        raise HrdInputError("'HRD_CALLER' cannot be empty")
    if is_missing(row['HRD_SCORE_TYPE']):
        raise HrdInputError("'HRD_SCORE_TYPE' cannot be empty")

    result = {
        'HRD_STATUS': None if is_missing(row.get('HRD_STATUS')) else _normalise_status(row['HRD_STATUS']),
        'HRD_STATUS_BASIS': None,
        'HRD_CALLER': row['HRD_CALLER'].strip(),
        'CALLER_VERSION': optional_text(row.get('CALLER_VERSION')),
        'HRD_SCORE_TYPE': _normalise_score_type(row['HRD_SCORE_TYPE']),
        'HRD_SCORE_DESCRIPTION': optional_text(row.get('HRD_SCORE_DESCRIPTION')),
        'HRD_SCORE': to_float(row.get('HRD_SCORE'), 'HRD_SCORE', HrdInputError),
        'HRD_SCORE_THRESHOLD_USED': to_float(
            row.get('HRD_SCORE_THRESHOLD_USED'), 'HRD_SCORE_THRESHOLD_USED', HrdInputError),
    }
    for component in _SCORE_COMPONENTS:
        result[component] = to_count(row.get(component), component, HrdInputError)
    if result['HRD_STATUS'] is None and result['HRD_SCORE'] is None:
        raise HrdInputError(f"HRD file ({os.path.basename(fname)}): 'HRD_STATUS' and 'HRD_SCORE' cannot both be empty")
    return result


def resolve_hrd_status(
        result: dict,
        threshold: Optional[float] = None,
        logger: Optional[logging.Logger] = None) -> dict:
    """
    Resolve the HRD status of a parsed HRD result (see parse_hrd_generic()):

      - with a threshold (--input_hrd_threshold): the status is derived from the
        score - HR-deficient if score >= threshold, otherwise HR-proficient. A status
        given in the file is then not used (a warning is given if it differs)
      - without a threshold: the status given in the file is used

    Raises HrdInputError if the status cannot be resolved (threshold, but no score;
    no threshold, and no status).
    """
    result = dict(result)
    if threshold is not None:
        if result['HRD_SCORE'] is None:
            raise HrdInputError(
                "Option '--input_hrd_threshold' requires a score ('HRD_SCORE') in the HRD file")
        derived = pcgr_vars.HRD_STATUS_DEFICIENT if result['HRD_SCORE'] >= threshold else pcgr_vars.HRD_STATUS_PROFICIENT
        if result['HRD_STATUS'] is not None and result['HRD_STATUS'] != derived and logger is not None:
            warn_message(
                f"HRD status in file ('{result['HRD_STATUS']}') differs from the status derived with "
                f"'--input_hrd_threshold' ({threshold:g}: '{derived}') - using the derived status", logger)
        result['HRD_STATUS'] = derived
        result['HRD_STATUS_BASIS'] = 'derived'
        result['HRD_SCORE_THRESHOLD_USED'] = float(threshold)
    else:
        if result['HRD_STATUS'] is None:
            raise HrdInputError(
                "HRD file has no HRD status ('HRD_STATUS') - please provide a threshold on the score "
                "with '--input_hrd_threshold' (thresholds are specific to the HRD assay/tool, hence "
                "there is no default)")
        result['HRD_STATUS_BASIS'] = 'caller'
    return result


def read_external_hrd(
        fname: str,
        threshold: Optional[float] = None,
        logger: Optional[logging.Logger] = None,
        quiet: bool = False) -> dict:
    """
    Read and validate an external HRD result file, and resolve the HRD status.
    Terminates with an error message (as for other PCGR input validation) if the
    file is invalid. Returns a dict with the normalised fields (see pcgr_vars.HRD_EXTERNAL_COLUMNS,
    except SAMPLE_ID). With quiet = True, only errors are reported (used for early
    validation of the input file; the result is logged when the file is read again
    in the HRD section of the workflow).

    threshold: user-defined threshold on the HRD score for HR-deficient
    (--input_hrd_threshold), see resolve_hrd_status().
    """
    if logger is None:
        logger = logging.getLogger("pcgr-hrd-external")
    try:
        result = resolve_hrd_status(
            parse_hrd_generic(fname), threshold=threshold, logger=None if quiet else logger)
    except HrdInputError as e:
        error_message(str(e), logger)

    if quiet:
        return result

    score_info = f", score = {result['HRD_SCORE']:g}" if result['HRD_SCORE'] is not None else ''
    basis = 'derived by PCGR' if result['HRD_STATUS_BASIS'] == 'derived' else 'reported in file'
    logger.info(
        f"External HRD result ({result['HRD_CALLER']}, {pcgr_vars.HRD_SCORE_TYPES[result['HRD_SCORE_TYPE']]}): "
        f"{result['HRD_STATUS']} ({basis}){score_info}")
    if result['HRD_STATUS_BASIS'] == 'derived':
        logger.info(
            f"HRD status derived with a user-defined threshold on the HRD score: "
            f"HR-deficient if score >= {result['HRD_SCORE_THRESHOLD_USED']:g}")
    return result


def write_external_hrd(result: dict, sample_id: str, fname: str) -> None:
    """Write a normalised external HRD result (one row) as a tab-separated file; missing values as '.'."""
    write_external_result(result, sample_id, fname, pcgr_vars.HRD_EXTERNAL_COLUMNS)


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(
        description="Read an external HRD result and print the normalised result")
    parser.add_argument("hrd_file", help="HRD result file")
    parser.add_argument("--sample_id", default="SAMPLE")
    parser.add_argument("--threshold", type=float, default=None,
                        help="threshold on the HRD score for HR-deficient (score >= threshold)")
    args = parser.parse_args()
    res = read_external_hrd(args.hrd_file, threshold=args.threshold)
    res['SAMPLE_ID'] = args.sample_id
    for col in pcgr_vars.HRD_EXTERNAL_COLUMNS:
        print(f"{col}\t{res.get(col) if res.get(col) is not None else '.'}")
