#!/usr/bin/env python
"""
Shared helpers for reading results of external assays/tools (MSI, HRD), given
as small tab-separated summary files with a header and a single data row.

Functions that validate content take the error class to raise (a subclass of
ExternalInputError, e.g. MsiInputError, HrdInputError), so that each module
reports errors with its own type.
"""

import os
from typing import Callable, Optional

import pandas as pd

from pcgr.utils import pd_to_csv


class ExternalInputError(ValueError):
    """Raised when an external result file cannot be read or has invalid content."""


def is_missing(value) -> bool:
    return value is None or str(value).strip() in ('', '.', 'NA', 'N/A', 'NaN', 'nan')


def optional_text(value) -> Optional[str]:
    return None if is_missing(value) else str(value).strip()


def to_float(value, name: str, error: type = ExternalInputError) -> Optional[float]:
    if is_missing(value):
        return None
    try:
        number = float(str(value).strip())
    except ValueError:
        raise error(f"'{name}' must be numeric (found: '{value}')")
    if number != number or number in (float('inf'), float('-inf')):
        raise error(f"'{name}' must be a finite number (found: '{value}')")
    return number


def to_count(value, name: str, error: type = ExternalInputError) -> Optional[int]:
    number = to_float(value, name, error)
    if number is None:
        return None
    if number < 0 or number != int(number):
        raise error(f"'{name}' must be a non-negative integer (found: '{value}')")
    return int(number)


def status_key(value) -> str:
    """Status value as compared with accepted spellings (upper-cased, '_'/' ' mapped to '-')."""
    return str(value).strip().upper().replace('_', '-').replace(' ', '-')


def match_synonym(value, synonyms: dict) -> Optional[str]:
    """
    Map a status value to its canonical form, using a dict of canonical status ->
    set of accepted spellings (see status_key()). Returns None if the value is not
    recognised.
    """
    key = status_key(value)
    for status, spellings in synonyms.items():
        if key in spellings:
            return status
    return None


def read_lines(fname: str, label: str, error: type = ExternalInputError) -> list:
    """Non-empty lines of a file (label, e.g. 'MSI', is used in error messages)."""
    if not os.path.isfile(fname) or os.stat(fname).st_size == 0:
        raise error(f"{label} file ({fname}) does not exist or is empty")
    with open(fname, mode='rt') as fh:
        return [line.rstrip('\r\n') for line in fh if line.strip() != '']


def split_tab(line: str) -> list:
    return [field.strip() for field in line.split('\t')]


def read_single_row_tsv(
        fname: str,
        label: str,
        required_cols: list,
        error: type = ExternalInputError,
        missing_hint: Optional[Callable[[list], str]] = None) -> dict:
    """
    Read a tab-separated file with a header and exactly one data row, and return
    the row as a dict (column -> value). Checks for required columns, duplicated
    column names, the number of data rows and the number of fields.

    missing_hint: optional function of the header, returning text appended to the
    error message on missing required columns.
    """
    base = os.path.basename(fname)
    lines = read_lines(fname, label, error)
    header = split_tab(lines[0])
    missing = [c for c in required_cols if c not in header]
    if missing:
        raise error(
            f"Missing required column(s) in {label} file ({base}): {', '.join(missing)}"
            + (missing_hint(header) if missing_hint is not None else ''))
    if len(set(header)) != len(header):
        raise error(f"Duplicated column names in {label} file ({base})")
    if len(lines) < 2:
        raise error(f"{label} file ({base}) has a header but no data row")
    if len(lines) > 2:
        raise error(
            f"{label} file ({base}) must contain a single data row (one sample), found {len(lines) - 1}")
    fields = split_tab(lines[1])
    if len(fields) != len(header):
        raise error(
            f"{label} file ({base}): data row has {len(fields)} field(s), "
            f"expected {len(header)} (as in the header)")
    return dict(zip(header, fields))


def write_external_result(result: dict, sample_id: str, fname: str, columns: list) -> None:
    """Write a normalised external result (one row) as a tab-separated file; missing values as '.'."""
    row = {'SAMPLE_ID': sample_id}
    row.update(result)
    df = pd.DataFrame([{c: row.get(c) for c in columns}], columns=columns)
    pd_to_csv(df.astype(object).where(df.notna(), None), fname, compression=None)
