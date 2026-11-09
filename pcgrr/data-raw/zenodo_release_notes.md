# Zenodo release steps

Publishes the example reports as a new version of the PCGR record on Zenodo,
using `data-raw/zenodo_release.py`. The README, site and output docs use the
concept DOI (`10.5281/zenodo.4022111`), so they don't need changing for a
release.

The record holds one zip, `v<version>.zip`, with a folder per sample. The
script builds it from a flat output folder, using the samples you name with
`--samples`, and builds it the same way every time, so `prepare` and `publish`
can compare checksums.

## Tokens

Set these in your shell. The two sites have separate accounts, and both tokens
need the `deposit:write` and `deposit:actions` scopes.

| Use         | Variable               | Create it at         |
|-------------|------------------------|----------------------|
| Testing     | `ZENODO_SANDBOX_TOKEN` | `sandbox.zenodo.org` |
| Real upload | `ZENODO_TOKEN`         | `zenodo.org`         |

## Steps

Replace `2.4.0` with the version being released. Run from `pcgrr/`.
`<dir>` is the Apptainer test output, e.g.
`workflow_test/output_data/v2.4.0/APPTAINER`, and `<samples>` the sample IDs
to publish, e.g. `WES-T001-APPTAINER-COAD-TN WES-T001-APPTAINER-BRCA-TO`.

1. Bump to the release version. The script checks `pcgrr/DESCRIPTION` and
   `pyproject.toml`.

2. Run the Apptainer test with the release code. The script refuses reports
   whose `.xlsx` or `conf.yaml` don't record this version. The test input
   VCFs carry a `##liftOverChainFile=/Users/...` header that ends up in the
   output VCFs; remove it from the inputs, or pass `--allow-local-paths`.

3. Run the local checks. This also builds the zip in a temporary folder and
   prints the description, generated from each sample's `conf.yaml`:

   ```bash
   python3 data-raw/zenodo_release.py check --version 2.4.0 --report-dir <dir> --samples <samples>
   ```

4. Test on the sandbox. The script hasn't run against the real API yet, so do
   this at least once before the first real release. `--new-record` is for a
   sandbox with no copy of the record.

   ```bash
   python3 data-raw/zenodo_release.py prepare --sandbox --new-record --version 2.4.0 --report-dir <dir> --samples <samples>
   python3 data-raw/zenodo_release.py publish --sandbox --draft-id <id> --version 2.4.0 --report-dir <dir> --samples <samples>
   ```

5. Create the real draft. This prints a review link and the draft ID, and
   nothing is public yet.

   ```bash
   python3 data-raw/zenodo_release.py prepare --version 2.4.0 --report-dir <dir> --samples <samples>
   ```

6. Review the draft, then publish with the same `--samples` (and `--exclude`,
   `--assembly`) as in step 5. You'll be asked to type `publish 2.4.0`.
   Published files can't be replaced, only superseded by a new version.

   ```bash
   python3 data-raw/zenodo_release.py publish --draft-id <id> --version 2.4.0 --report-dir <dir> --samples <samples>
   ```

## Options

- `--assembly grch37`: publish the GRCh37 reports (default `grch38`).
- `--exclude SUFFIX`: leave out files with this suffix, e.g. `--exclude log`
  (repeatable).
- `--allow-local-paths`: publish even if files contain paths like `/Users/<name>`.

Samples with a germline section must use artificial CPSR input (`-ART-`), and
the script refuses any other.
