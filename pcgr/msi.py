#!/usr/bin/env python
"""
Import of microsatellite instability (MSI) results from external MSI callers.

An external MSI result is read from a single summary file and normalised to a
common representation (one row, see NORMALISED_COLUMNS):

  - 'generic'       tab-separated file with one row and the columns described in
                    parse_msi_generic() (required: MSI_STATUS, MSI_CALLER)
  - 'msisensor-pro' native summary file of MSIsensor-pro (v1.x; the file with columns
                    'Total_Number_of_Sites', 'Number_of_Unstable_Sites' and '%',
                    written by both the 'pro' (tumor-only) and 'msi' (tumor-normal)
                    commands). The tool reports no MSI status, so the status is
                    derived from the percentage of unstable sites (see
                    MSI_CALLER_DEFAULT_THRESHOLDS; the threshold can be set by the user)

Callers differ in what they report: MSIsensor-pro gives counts of evaluated and
unstable sites but no status, the status of other callers is given directly, and
the scores are on different scales (hence MSI_SCORE_TYPE). Optional fields that
a caller does not report are left empty.

File names are not used for any decision - the content (header) of the file is.
"""

import logging
import os
from typing import Optional

import pandas as pd

from pcgr.utils import error_message, pd_to_csv

MSI_CALLERS = ['generic', 'msisensor-pro']

MSI_STATUS_H = 'MSI-H'
MSI_STATUS_MSS = 'MSS'
MSI_STATUS_INDETERMINATE = 'Indeterminate'

## Default threshold on the MSI score for callers that report a score but no MSI status
## (the status is then derived by PCGR). For MSIsensor-pro: percentage of unstable sites at
## or above which a result is called MSI-H. MSIsensor-pro itself does not report a status
## (nor define a cutoff) - this is the cutoff used in practice by users of the tool. Practice
## varies (assay, tumor-only versus tumor-normal), hence the threshold can be set by the user
## (--input_msi_threshold). Callers that report their own status (e.g. MANTIS, generic format)
## have no entry here - their status is used as is.
MSI_CALLER_DEFAULT_THRESHOLDS = {'msisensor-pro': 20.0}
MSISENSOR_PRO_MSIH_CUTOFF = MSI_CALLER_DEFAULT_THRESHOLDS['msisensor-pro']

MSISENSOR_PRO_SUMMARY_HEADER = ['Total_Number_of_Sites', 'Number_of_Unstable_Sites', '%']

NORMALISED_COLUMNS = [
    'SAMPLE_ID',
    'MSI_STATUS',            ## MSI-H, MSS or Indeterminate
    'MSI_STATUS_BASIS',      ## 'caller' (status reported by caller) or 'derived' (by PCGR from the score)
    'MSI_CALLER',
    'CALLER_VERSION',
    'MSI_SCORE',
    'MSI_SCORE_TYPE',
    'MSI_SCORE_THRESHOLD_USED',  ## threshold on the score for MSI-H (by the caller, or by PCGR if status is derived)
    'N_LOCI_EVALUATED',
    'N_LOCI_UNSTABLE',
]

## accepted spellings of the MSI status (compared after upper-casing and mapping '_'/' ' to '-')
_STATUS_SYNONYMS = {
    MSI_STATUS_H: {'MSI-H', 'MSIH', 'MSI-HIGH', 'MSI', 'UNSTABLE'},
    ## MSI-L is regarded as microsatellite stable (not MSI-H), as in clinical practice
    MSI_STATUS_MSS: {'MSS', 'STABLE', 'MSI-L', 'MSIL', 'MSI-LOW'},
    MSI_STATUS_INDETERMINATE: {'INDETERMINATE', 'UNDETERMINED', 'EQUIVOCAL', 'INCONCLUSIVE', 'UNKNOWN'},
}
_LOW_STATUS_SPELLINGS = {'MSI-L', 'MSIL', 'MSI-LOW'}


class MsiInputError(ValueError):
    """Raised when an external MSI result file cannot be read or has invalid content."""


def _is_missing(value) -> bool:
    return value is None or str(value).strip() in ('', '.', 'NA', 'N/A', 'NaN', 'nan')


def _to_float(value, name: str) -> Optional[float]:
    if _is_missing(value):
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        raise MsiInputError(f"'{name}' must be numeric (found: '{value}')")


def _to_count(value, name: str) -> Optional[int]:
    number = _to_float(value, name)
    if number is None:
        return None
    if number < 0 or number != int(number):
        raise MsiInputError(f"'{name}' must be a non-negative integer (found: '{value}')")
    return int(number)


def _normalise_status(value, logger: Optional[logging.Logger] = None) -> str:
    key = str(value).strip().upper().replace('_', '-').replace(' ', '-')
    for status, spellings in _STATUS_SYNONYMS.items():
        if key in spellings:
            if key in _LOW_STATUS_SPELLINGS and logger is not None:
                logger.info(f"MSI status '{value}' is regarded as '{MSI_STATUS_MSS}' (not MSI-H)")
            return status
    raise MsiInputError(
        f"Unrecognised MSI status ('{value}') - must be one of: "
        f"{MSI_STATUS_H}, {MSI_STATUS_MSS}, {MSI_STATUS_INDETERMINATE} "
        "(or a common synonym, e.g. 'Unstable'/'Stable')")


def _optional_text(value) -> Optional[str]:
    return None if _is_missing(value) else str(value).strip()


def _read_lines(fname: str) -> list:
    if not os.path.isfile(fname) or os.stat(fname).st_size == 0:
        raise MsiInputError(f"MSI file ({fname}) does not exist or is empty")
    with open(fname, mode='rt') as fh:
        return [line.rstrip('\r\n') for line in fh if line.strip() != '']


def _split_tab(line: str) -> list:
    return [field.strip() for field in line.split('\t')]


def _hint_for_wrong_file(header: list, expected_caller: str) -> str:
    """Suggest the likely cause when a file does not match the expected caller."""
    first = header[0] if header else ''
    if first.startswith('Average Metric Value'):
        return " This looks like a MANTIS status file ('--input_msi_caller mantis')."
    if first == 'Total_Number_of_Sites' and expected_caller != 'msisensor-pro':
        return " This looks like an MSIsensor-pro summary file ('--input_msi_caller msisensor-pro')."
    if first == 'chromosome' and 'location' in header:
        return (" This looks like a locus-level MSIsensor-pro file (_all/_unstable) - "
                "please provide the summary file (header: " + ', '.join(MSISENSOR_PRO_SUMMARY_HEADER) + ").")
    if first == 'MSI_STATUS' or 'MSI_STATUS' in header:
        return " This looks like a file in the generic format ('--input_msi_caller generic')."
    return ""


def parse_msi_generic(fname: str, logger: Optional[logging.Logger] = None) -> dict:
    """
    Read an MSI result in the generic format: a tab-separated file with a header
    and exactly one data row (one sample). Columns:

      MSI_STATUS        required - MSI-H, MSS or Indeterminate (common synonyms such as
                        'Unstable'/'Stable' are accepted, case-insensitive)
      MSI_CALLER        required - name of the caller/method
      CALLER_VERSION    optional
      MSI_SCORE         optional - numeric score as reported by the caller
      MSI_SCORE_TYPE    optional - description of the score, e.g. 'percent unstable sites'
      MSI_SCORE_THRESHOLD_USED
                        optional - threshold on the score applied by the caller for MSI-H
      N_LOCI_EVALUATED  optional - number of microsatellite loci evaluated
      N_LOCI_UNSTABLE   optional - number of unstable loci

    Additional columns are ignored. Missing optional values may be left empty
    or set to '.' / 'NA'.

    Raises MsiInputError if the file is invalid.
    """
    lines = _read_lines(fname)
    header = _split_tab(lines[0])
    missing = [c for c in ('MSI_STATUS', 'MSI_CALLER') if c not in header]
    if missing:
        raise MsiInputError(
            f"Missing required column(s) in MSI file ({os.path.basename(fname)}): {', '.join(missing)}."
            + _hint_for_wrong_file(header, 'generic'))
    if len(set(header)) != len(header):
        raise MsiInputError(f"Duplicated column names in MSI file ({os.path.basename(fname)})")
    if len(lines) < 2:
        raise MsiInputError(f"MSI file ({os.path.basename(fname)}) has a header but no data row")
    if len(lines) > 2:
        raise MsiInputError(
            f"MSI file ({os.path.basename(fname)}) must contain a single data row (one sample), "
            f"found {len(lines) - 1}")

    fields = _split_tab(lines[1])
    if len(fields) != len(header):
        raise MsiInputError(
            f"MSI file ({os.path.basename(fname)}): data row has {len(fields)} field(s), "
            f"expected {len(header)} (as in the header)")
    row = dict(zip(header, fields))

    if _is_missing(row['MSI_CALLER']):
        raise MsiInputError("'MSI_CALLER' cannot be empty")
    if _is_missing(row['MSI_STATUS']):
        raise MsiInputError("'MSI_STATUS' cannot be empty")

    result = {
        'MSI_STATUS': _normalise_status(row['MSI_STATUS'], logger),
        'MSI_STATUS_BASIS': 'caller',
        'MSI_CALLER': row['MSI_CALLER'].strip(),
        'CALLER_VERSION': _optional_text(row.get('CALLER_VERSION')),
        'MSI_SCORE': _to_float(row.get('MSI_SCORE'), 'MSI_SCORE'),
        'MSI_SCORE_TYPE': _optional_text(row.get('MSI_SCORE_TYPE')),
        'MSI_SCORE_THRESHOLD_USED': _to_float(
            row.get('MSI_SCORE_THRESHOLD_USED'), 'MSI_SCORE_THRESHOLD_USED'),
        'N_LOCI_EVALUATED': _to_count(row.get('N_LOCI_EVALUATED'), 'N_LOCI_EVALUATED'),
        'N_LOCI_UNSTABLE': _to_count(row.get('N_LOCI_UNSTABLE'), 'N_LOCI_UNSTABLE'),
    }
    if result['N_LOCI_EVALUATED'] is not None and result['N_LOCI_UNSTABLE'] is not None \
            and result['N_LOCI_UNSTABLE'] > result['N_LOCI_EVALUATED']:
        raise MsiInputError(
            f"'N_LOCI_UNSTABLE' ({result['N_LOCI_UNSTABLE']}) cannot exceed "
            f"'N_LOCI_EVALUATED' ({result['N_LOCI_EVALUATED']})")
    return result


def parse_msisensor_pro(
        fname: str,
        msih_threshold: float = MSISENSOR_PRO_MSIH_CUTOFF,
        logger: Optional[logging.Logger] = None) -> dict:
    """
    Read the native summary file of MSIsensor-pro (header: Total_Number_of_Sites,
    Number_of_Unstable_Sites, %; one data row). The file has no sample identifier
    and no MSI status. The status is derived from the fraction of unstable sites:

      - no evaluated sites               -> Indeterminate
      - unstable sites >= threshold (%)  -> MSI-H
      - otherwise                        -> MSS

    The score is the percentage of unstable sites as reported by the tool. With the
    small number of sites typical for panel data, the percentage changes in coarse
    steps, so the counts are kept alongside the score.

    Raises MsiInputError if the file is not an MSIsensor-pro summary file.
    """
    lines = _read_lines(fname)
    header = _split_tab(lines[0])
    if header != MSISENSOR_PRO_SUMMARY_HEADER:
        raise MsiInputError(
            f"MSI file ({os.path.basename(fname)}) is not an MSIsensor-pro summary file "
            f"(expected header: {', '.join(MSISENSOR_PRO_SUMMARY_HEADER)})."
            + _hint_for_wrong_file(header, 'msisensor-pro'))
    if len(lines) < 2:
        raise MsiInputError(f"MSIsensor-pro summary file ({os.path.basename(fname)}) has a header but no data row")
    if len(lines) > 2:
        raise MsiInputError(
            f"MSIsensor-pro summary file ({os.path.basename(fname)}) must contain a single data row "
            f"(one sample), found {len(lines) - 1}")

    fields = _split_tab(lines[1])
    if len(fields) != 3:
        raise MsiInputError(
            f"MSIsensor-pro summary file ({os.path.basename(fname)}): data row has {len(fields)} field(s), expected 3")
    n_sites = _to_count(fields[0], 'Total_Number_of_Sites')
    n_unstable = _to_count(fields[1], 'Number_of_Unstable_Sites')
    pct_reported = _to_float(fields[2], '%')
    if n_sites is None or n_unstable is None or pct_reported is None:
        raise MsiInputError(f"MSIsensor-pro summary file ({os.path.basename(fname)}) has missing values")
    if n_unstable > n_sites:
        raise MsiInputError(
            f"MSIsensor-pro summary: number of unstable sites ({n_unstable}) exceeds "
            f"total number of sites ({n_sites})")

    result = {
        'MSI_STATUS': MSI_STATUS_INDETERMINATE,
        'MSI_STATUS_BASIS': 'derived',
        'MSI_CALLER': 'MSIsensor-pro',
        'CALLER_VERSION': None,
        'MSI_SCORE': None,
        'MSI_SCORE_TYPE': 'percent unstable sites',
        'MSI_SCORE_THRESHOLD_USED': float(msih_threshold),
        'N_LOCI_EVALUATED': n_sites,
        'N_LOCI_UNSTABLE': n_unstable,
    }
    if n_sites == 0:
        if logger is not None:
            logger.warning("MSIsensor-pro summary reports no evaluated sites - MSI status cannot be determined")
        return result

    result['MSI_SCORE'] = pct_reported
    if abs(pct_reported - 100.0 * n_unstable / n_sites) > 0.1:
        if logger is not None:
            logger.warning(
                f"MSIsensor-pro summary: reported percentage ({pct_reported}) differs from "
                f"{n_unstable}/{n_sites} unstable sites - using the reported percentage as score")
    ## compare using integer counts (not the rounded percentage) to avoid rounding artefacts at the threshold
    result['MSI_STATUS'] = MSI_STATUS_H if n_unstable * 100 >= msih_threshold * n_sites else MSI_STATUS_MSS
    return result


def read_external_msi(
        fname: str,
        caller: str = 'generic',
        logger: Optional[logging.Logger] = None,
        quiet: bool = False,
        msih_threshold: Optional[float] = None) -> dict:
    """
    Read and validate an external MSI result file from the given caller
    ('generic' or 'msisensor-pro'). Terminates with an error message (as for other
    PCGR input validation) if the file is invalid. Returns a dict with the
    normalised fields (see NORMALISED_COLUMNS, except SAMPLE_ID). With quiet = True,
    only errors are reported (used for early validation of the input file; the result
    is logged when the file is read again in the MSI section of the workflow).

    msih_threshold: user-defined threshold on the MSI score for MSI-H, only for callers
    that report a score but no status (see MSI_CALLER_DEFAULT_THRESHOLDS; default if None).
    """
    if logger is None:
        logger = logging.getLogger("pcgr-msi")
    if caller not in MSI_CALLERS:
        error_message(f"Unsupported MSI caller ('{caller}') - must be one of: {', '.join(MSI_CALLERS)}", logger)
    threshold_error = check_msi_threshold(caller, msih_threshold)
    if threshold_error is not None:
        error_message(threshold_error, logger)
    try:
        if caller == 'msisensor-pro':
            result = parse_msisensor_pro(
                fname, msih_threshold=MSI_CALLER_DEFAULT_THRESHOLDS[caller] if msih_threshold is None else msih_threshold,
                logger=None if quiet else logger)
        else:
            result = parse_msi_generic(fname, logger=None if quiet else logger)
    except MsiInputError as e:
        error_message(str(e), logger)

    if quiet:
        return result

    n_info = ''
    if result['N_LOCI_EVALUATED'] is not None:
        n_info = f", {result['N_LOCI_UNSTABLE']}/{result['N_LOCI_EVALUATED']} unstable loci"
    score_info = f", score = {result['MSI_SCORE']}" if result['MSI_SCORE'] is not None else ''
    basis = 'derived by PCGR' if result['MSI_STATUS_BASIS'] == 'derived' else 'reported by caller'
    logger.info(
        f"External MSI result ({result['MSI_CALLER']}): {result['MSI_STATUS']} ({basis}){score_info}{n_info}")
    if result['MSI_STATUS_BASIS'] == 'derived' and result['MSI_SCORE_THRESHOLD_USED'] is not None:
        source = 'default' if msih_threshold is None else 'user-defined'
        logger.info(
            f"MSI status derived with a {source} threshold on the MSI score: "
            f"MSI-H if score >= {result['MSI_SCORE_THRESHOLD_USED']:g} ({result['MSI_SCORE_TYPE']})")
    return result


def check_msi_threshold(caller: str, msih_threshold: Optional[float]) -> Optional[str]:
    """
    Check a user-defined threshold on the MSI score (--input_msi_threshold). Returns an
    error message if the threshold cannot be used with the caller, or is out of range
    (percentage of unstable sites: > 0 and <= 100), otherwise None.
    """
    if msih_threshold is None:
        return None
    if caller not in MSI_CALLER_DEFAULT_THRESHOLDS:
        return (
            f"Option '--input_msi_threshold' cannot be used with the '{caller}' MSI caller - it applies "
            f"only to callers that report a score but no MSI status "
            f"({', '.join(MSI_CALLER_DEFAULT_THRESHOLDS)}). For other callers, the MSI status "
            f"reported by the caller is used")
    if not 0 < msih_threshold <= 100:
        return (
            f"Option '--input_msi_threshold' ({msih_threshold:g}) must be a percentage of unstable "
            f"sites: larger than 0 and not larger than 100")
    return None


def write_external_msi(result: dict, sample_id: str, fname: str) -> None:
    """Write a normalised external MSI result (one row) as a tab-separated file; missing values as '.'."""
    row = {'SAMPLE_ID': sample_id}
    row.update(result)
    df = pd.DataFrame([{c: row.get(c) for c in NORMALISED_COLUMNS}], columns=NORMALISED_COLUMNS)
    pd_to_csv(df.astype(object).where(df.notna(), None), fname, compression=None)


if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(
        description="Read an external MSI result (generic format or MSIsensor-pro summary) and print the normalised result")
    parser.add_argument("msi_file", help="MSI summary file")
    parser.add_argument("--caller", choices=MSI_CALLERS, default="generic")
    parser.add_argument("--sample_id", default="SAMPLE")
    parser.add_argument("--threshold", type=float, default=None,
                        help="threshold on the MSI score for MSI-H (MSIsensor-pro: percent unstable sites; default 20)")
    args = parser.parse_args()
    res = read_external_msi(args.msi_file, args.caller, msih_threshold=args.threshold)
    res['SAMPLE_ID'] = args.sample_id
    for col in NORMALISED_COLUMNS:
        print(f"{col}\t{res.get(col) if res.get(col) is not None else '.'}")
