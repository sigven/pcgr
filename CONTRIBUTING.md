# Contributing to PCGR

Thanks for your interest in contributing to the Personal Cancer Genome Reporter
(PCGR). Contributions of all kinds are welcome — bug reports, documentation
improvements, feature suggestions, and code.

This document explains how to report problems and how to propose changes.

## Code of conduct

This project follows a [Code of Conduct](CODE_OF_CONDUCT.md). By participating,
you are expected to uphold it. Please report unacceptable behaviour to the
maintainers.

## Ways to contribute

- **Report a bug** — open an issue (see below for what to include)
- **Request a feature** — open an issue describing the use case
- **Ask a question** — use GitHub Discussions or open an issue
- **Improve documentation** — corrections and clarifications are very welcome
- **Contribute code** — see *Development workflow* below

## Reporting bugs

Most issues we receive can only be diagnosed with the full run context, so
please include **all** of the following. Issues without this information will
usually need a follow-up before we can help.

1. **PCGR version** (`pcgr --version`)
2. **Reference data bundle version** and the genome assembly used
   (`grch37` or `grch38`)
3. **Installation method** — Conda, Docker, Singularity/Apptainer —
   including the image tag or environment specification
4. **Operating system** and, if relevant, the compute environment
   (laptop, HPC/cluster, cloud)
5. **The exact command you ran**, in full, with all arguments
6. **The complete error message and log output** (please paste as text in a
   fenced code block rather than as a screenshot)
7. **Input characteristics** — variant caller used, approximate number of
   variants, tumour type, whether CNA/expression/fusion input was supplied
8. **A minimal example input** that reproduces the problem, if you are able to
   share one

Please do not attach patient-identifiable data to a public issue. If a problem
can only be reproduced with sensitive data, say so in the issue and we will
find another way to investigate.

## Requesting features

When proposing a feature, describe the scientific or clinical use case rather
than only the implementation. PCGR aims to produce clinically interpretable
output aligned with established guidelines (AMP/ASCO/CAP, ClinGen/CGC/VICC), so
proposals that affect variant classification or tier assignment should reference
the relevant guideline or evidence source.

## Development workflow

We use a simple two-branch model:

- **`main`** — reflects the current released version. Do not open pull requests
  against `main`.
- **`dev`** — the active development branch. **All pull requests should target
  `dev`.**

To contribute code:

1. Fork the repository and create a branch from `dev`
   (e.g. `feature/short-description` or `fix/short-description`)
2. Make your changes, keeping the pull request focused on a single concern
3. Update documentation and the changelog where relevant
4. Open a pull request against `dev`, describing **what** the change does and
   **why** it is needed

For anything substantial — new annotation sources, changes to classification or
tiering logic, or architectural changes — please open an issue to discuss the
approach before writing code. This saves time on both sides.

## Code organisation

PCGR spans two languages, and it helps to know where a change belongs:

- **Python** — the pipeline, variant annotation, and much of the interpretation
  logic (including oncogenicity classification and biomarker matching)
- **R** (the [`pcgrr`](https://github.com/sigven/pcgrr) package) — report
  generation and visualisation via Quarto, together with parts of the clinical
  interpretation logic, notably AMP/ASCO/CAP tier assignment

Note that this split is not absolute: some classification logic currently lives
in `pcgrr` alongside the reporting code. If you are unsure where a change
belongs, open an issue and ask before starting work.

## Coding conventions

- Follow the style of the surrounding code
- Python: follow PEP 8 where practical; add docstrings to new functions
- R: follow the existing style in `pcgrr`
- Keep commits reasonably self-contained with descriptive messages
- Do **not** commit reference data bundles, large binary files, or test data
  containing patient information

## Testing and validation

Automated test coverage is currently limited and is actively being expanded. In
the meantime, we ask contributors to describe how a change was validated:

- State the command(s) you ran to exercise the change
- For changes affecting annotation, classification or tiering, show the effect
  on output — for example, the relevant rows or fields before and after
- Confirm that an end-to-end run completes on example data
- Where tests exist, please run them and add new ones covering your change

Continuous integration runs automatically on pull requests. Please make sure it
passes before requesting review.

## Releases and versioning

PCGR follows semantic versioning. Releases are cut from `dev` into `main` by the
maintainers, accompanied by an updated changelog. Reference data bundles are
versioned separately and tied to specific releases — see the documentation for
the compatible bundle for each version.

## Licence

PCGR is released under the MIT licence. By contributing, you agree that your
contributions will be licensed under the same terms.

## Getting in touch

- **Bugs and feature requests** — GitHub Issues
- **Questions and general discussion** — GitHub Discussions (to come)
- **Documentation** — https://sigven.github.io/pcgr/

Thanks for your effort to improve PCGR!
